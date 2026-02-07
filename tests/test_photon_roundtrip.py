"""Test photon data roundtrip: HDF5 → OpenMC object → Arrow → readback."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import openmc.data

from nuclear_data_to_yamc_format import export_photon_to_arrow, read_photon_from_arrow, verify_photon


class TestPhotonRoundtrip:
    """Test exporting and re-reading photon data preserves all values."""

    def test_basic_element(self, photon_h5_path):
        """Test roundtrip with a basic element."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.arrow"
            export_photon_to_arrow(data, arrow_path)

            # Check files exist
            assert (arrow_path / "element.feather").exists()

            # Read back
            arrow_data = read_photon_from_arrow(arrow_path)

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
            arrow_path = Path(tmpdir) / f"{data.name}.photon.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            elem = arrow_data["element"]

            # Build union grid
            union_grid = np.array([])
            for rx in data:
                union_grid = np.union1d(union_grid, rx.xs.x)

            np.testing.assert_array_equal(
                np.array(elem["energy"]), union_grid)

    def test_subshells(self, photon_h5_path):
        """Test that subshell data is preserved."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            if arrow_data["subshells"]:
                for sub in arrow_data["subshells"]:
                    assert sub["designator"] is not None
                    assert len(sub["xs"]) > 0

    def test_compton_profiles(self, photon_h5_path):
        """Test that Compton profile data is preserved."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        if not data.compton_profiles:
            pytest.skip("No Compton profile data")

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.arrow"
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

    def test_bremsstrahlung(self, photon_h5_path):
        """Test that bremsstrahlung data is preserved."""
        data = openmc.data.IncidentPhoton.from_hdf5(photon_h5_path)

        if not data.bremsstrahlung:
            pytest.skip("No bremsstrahlung data")

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.photon.arrow"
            export_photon_to_arrow(data, arrow_path)
            arrow_data = read_photon_from_arrow(arrow_path)

            assert "bremsstrahlung" in arrow_data
            abrem = arrow_data["bremsstrahlung"]
            brem = data.bremsstrahlung

            assert abrem["I"] == float(brem['I'])
            np.testing.assert_array_equal(
                np.array(abrem["electron_energy"]),
                np.array(brem['electron_energy']))
