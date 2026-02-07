"""Test neutron data roundtrip: HDF5 → OpenMC object → Arrow → readback."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import openmc.data

from nuclear_data_to_yamc_format import export_neutron_to_arrow, read_neutron_from_arrow, verify_neutron


class TestNeutronRoundtrip:
    """Test exporting and re-reading neutron data preserves all values."""

    def test_simple_nuclide(self, neutron_h5_path):
        """Test roundtrip with a simple nuclide (e.g., Li6, H1)."""
        data = openmc.data.IncidentNeutron.from_hdf5(neutron_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.arrow"
            export_neutron_to_arrow(data, arrow_path)

            # Check files exist
            assert (arrow_path / "nuclide.feather").exists()
            assert (arrow_path / "reactions.feather").exists()

            # Read back
            arrow_data = read_neutron_from_arrow(arrow_path)

            # Verify metadata
            nuc = arrow_data["nuclide"]
            assert nuc["name"] == data.name
            assert nuc["Z"] == data.atomic_number
            assert nuc["A"] == data.mass_number
            assert nuc["metastable"] == data.metastable
            assert nuc["atomic_weight_ratio"] == data.atomic_weight_ratio

            # Verify temperatures
            assert nuc["temperatures"] == list(data.temperatures)

            # Verify energy grids
            for i, t in enumerate(nuc["energy_temperatures"]):
                np.testing.assert_array_equal(
                    np.array(nuc["energy_values"][i]),
                    np.array(data.energy[t])
                )

            # Verify reactions exist
            assert len(arrow_data["reactions"]) > 0

            # Full verification
            assert verify_neutron(data, arrow_path)

    def test_fissile_nuclide(self, fissile_h5_path):
        """Test roundtrip with a fissile nuclide (U235, URR, total_nu)."""
        data = openmc.data.IncidentNeutron.from_hdf5(fissile_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.arrow"
            export_neutron_to_arrow(data, arrow_path)

            # Check files exist
            assert (arrow_path / "nuclide.feather").exists()
            assert (arrow_path / "reactions.feather").exists()

            # Check URR if present
            if data.urr:
                assert (arrow_path / "urr.feather").exists()
                arrow_data = read_neutron_from_arrow(arrow_path)
                assert "urr" in arrow_data
                assert len(arrow_data["urr"]) == len(data.urr)
                for i, (temp, urr) in enumerate(data.urr.items()):
                    au = arrow_data["urr"][i]
                    assert au["temperature"] == temp
                    np.testing.assert_array_equal(
                        np.array(au["energy"]),
                        np.array(urr.energy)
                    )

            # Check total_nu if present
            has_derived = any(
                len(rx.derived_products) > 0
                for rx in data.reactions.values()
            )
            if has_derived:
                assert (arrow_path / "total_nu.feather").exists()

            # Full verification
            assert verify_neutron(data, arrow_path)

    def test_cross_section_values(self, neutron_h5_path):
        """Test that cross section values are bit-identical."""
        data = openmc.data.IncidentNeutron.from_hdf5(neutron_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.arrow"
            export_neutron_to_arrow(data, arrow_path)
            arrow_data = read_neutron_from_arrow(arrow_path)

            for rx_orig, rx_arrow in zip(
                data.reactions.values(),
                arrow_data["reactions"]
            ):
                for t_idx, T in enumerate(rx_arrow["xs_temperatures"]):
                    if rx_orig.xs[T] is not None:
                        expected = np.array(rx_orig.xs[T].y)
                        actual = np.array(rx_arrow["xs_values"][t_idx])
                        np.testing.assert_array_equal(
                            expected, actual,
                            err_msg=f"XS mismatch for MT={rx_orig.mt} T={T}")

    def test_product_yields(self, neutron_h5_path):
        """Test that product yield data is preserved."""
        data = openmc.data.IncidentNeutron.from_hdf5(neutron_h5_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            arrow_path = Path(tmpdir) / f"{data.name}.arrow"
            export_neutron_to_arrow(data, arrow_path)
            arrow_data = read_neutron_from_arrow(arrow_path)

            prod_idx = 0
            for rx in data.reactions.values():
                for p in rx.products:
                    if prod_idx >= len(arrow_data["products"]):
                        break
                    ap = arrow_data["products"][prod_idx]
                    assert ap["particle"] == str(p.particle)
                    assert ap["emission_mode"] == str(p.emission_mode)
                    assert ap["decay_rate"] == float(p.decay_rate)
                    prod_idx += 1
