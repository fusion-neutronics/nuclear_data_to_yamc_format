"""Test photon data roundtrip: HDF5 -> OpenMC object -> .yamc.arrow/ -> readback."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest
import openmc.data

from nuclear_data_to_yamc_format import (
    export_photon_to_arrow, read_photon_from_arrow, verify_photon,
)


class TestPhotonRoundtrip:
    """Test exporting and re-reading photon data preserves all values."""

    def test_basic_element(self, photon_h5_path):
        """Test roundtrip with a basic element."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path)

            # Check files exist
            assert (arrow_path / "version.json").exists()
            assert (arrow_path / "element.arrow").exists()

            # Read back
            arrow_data = read_photon_from_arrow(arrow_path)

            # Verify version.json
            assert "version" in arrow_data
            assert arrow_data["version"]["format_version"] == 1

            # Verify metadata
            elem = arrow_data["element"]
            assert elem["name"] == data.name
            assert elem["Z"] == data.atomic_number

            # Full verification
            assert verify_photon(data, arrow_path)

    def test_cross_sections(self, photon_h5_path):
        """Test that photon cross sections match."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            elem = arrow_data["element"]

            # Build union grid
            union_grid = np.array([])
            for rx in data:
                union_grid = np.union1d(union_grid, rx.xs.x)

            np.testing.assert_array_equal(
                np.array(elem["energy"]), union_grid)

    def test_log_space_data(self, photon_h5_path):
        """Test that log-space energy and XS are present and correct."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            elem = arrow_data["element"]

            # Check ln_energy exists and matches log of energy
            assert "ln_energy" in elem
            energy = np.array(elem["energy"])
            ln_energy = np.array(elem["ln_energy"])
            expected_ln = np.log(energy)
            np.testing.assert_allclose(ln_energy, expected_ln, rtol=1e-14)

    def test_subshells(self, photon_h5_path):
        """Test that subshell data is preserved with log-space XS."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            if arrow_data["subshells"]:
                for sub in arrow_data["subshells"]:
                    assert sub["designator"] is not None
                    assert len(sub["xs"]) > 0
                    # Verify ln_xs is present
                    assert "ln_xs" in sub
                    if sub["ln_xs"] is not None and len(sub["ln_xs"]) > 0:
                        xs = np.array(sub["xs"])
                        ln_xs = np.array(sub["ln_xs"])
                        # Log of XS values (with safe handling of zeros)
                        mask = xs > 0
                        if np.any(mask):
                            np.testing.assert_allclose(
                                ln_xs[mask], np.log(xs[mask]),
                                rtol=1e-14)

    def test_compton_profiles(self, photon_h5_path):
        """Test that Compton profile data is preserved with CDFs."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        if not data.compton_profiles:
            pytest.skip("No Compton profile data")

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            assert "compton" in arrow_data
            cmp = arrow_data["compton"]
            profile = data.compton_profiles

            np.testing.assert_array_equal(
                np.array(cmp["num_electrons"]),
                np.array(profile['num_electrons']))
            np.testing.assert_array_equal(
                np.array(cmp["binding_energy"]),
                np.array(profile['binding_energy']))
            np.testing.assert_array_equal(
                np.array(cmp["pz"]),
                np.array(profile['J'][0].x))

            # Verify CDFs exist and are valid
            assert "J_cdf_data" in cmp
            if cmp["J_cdf_data"] is not None and len(cmp["J_cdf_data"]) > 0:
                cdf = np.array(cmp["J_cdf_data"]).reshape(cmp["J_cdf_shape"])
                # CDFs should end at 1.0
                for s in range(cdf.shape[0]):
                    np.testing.assert_allclose(cdf[s, -1], 1.0, atol=1e-10)
                    # Should be monotonically non-decreasing
                    assert np.all(cdf[s, 1:] >= cdf[s, :-1] - 1e-15)

    def test_bremsstrahlung(self, photon_h5_path):
        """Test that bremsstrahlung data is preserved."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        if not data.bremsstrahlung:
            pytest.skip("No bremsstrahlung data")

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            assert "bremsstrahlung" in arrow_data
            abrem = arrow_data["bremsstrahlung"]
            brem = data.bremsstrahlung

            assert abrem["I"] == float(brem['I'])
            np.testing.assert_array_equal(
                np.array(abrem["electron_energy"]),
                np.array(brem['electron_energy']))

    def test_version_json(self, photon_h5_path):
        """Test that version.json contains expected fields."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.yamc.arrow"
            export_photon_to_arrow(data, arrow_path, library="test-lib")

            version = json.loads((arrow_path / "version.json").read_text())
            assert version["format_version"] == 1
            assert version["library"] == "test-lib"
            assert "converter_version" in version
            assert "created_utc" in version
