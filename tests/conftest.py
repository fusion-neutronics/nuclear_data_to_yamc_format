"""Shared test fixtures for nuclear_data_to_yamc_format tests."""

import os
from pathlib import Path

import pytest

# Known locations with HDF5 nuclear data
_SEARCH_DIRS = [
    Path("/home/jon/yamc-org/yamc/tests"),
    Path("/home/jon/yamc/tests"),
    Path("/home/jon/cross_section_data_fendl_3.2c/fendl-3.2c-hdf5"),
    Path("/home/jon/tendl-2021-hdf5/tendl-2021-hdf5"),
    Path("/home/jon/nuclear_data/fendl-3.2c-hdf5/photon"),
    Path("/home/jon/nuclear_data/endf-b8.0-hdf5/photon"),
]


def _find_hdf5_file(name):
    """Search common locations for an HDF5 nuclear data file."""
    # Check environment variable first
    if os.environ.get("OPENMC_CROSS_SECTIONS"):
        xs_xml = Path(os.environ["OPENMC_CROSS_SECTIONS"])
        if xs_xml.exists():
            lib_root = xs_xml.parent
            for sub in ["neutron", "photon", ""]:
                candidate = lib_root / sub / name
                if candidate.exists():
                    return candidate

    # Check known directories
    for d in _SEARCH_DIRS:
        candidate = d / name
        if candidate.exists():
            return candidate

    return None


@pytest.fixture
def neutron_h5_path():
    """Path to a neutron HDF5 file for testing."""
    for name in ["Li6.h5", "H1.h5", "Fe56.h5"]:
        p = _find_hdf5_file(name)
        if p is not None:
            return p
    pytest.skip("No neutron HDF5 file found; set OPENMC_CROSS_SECTIONS")


@pytest.fixture
def fissile_h5_path():
    """Path to a fissile nuclide HDF5 for testing (U235)."""
    for name in ["U235.h5", "U238.h5", "Pu239.h5"]:
        p = _find_hdf5_file(name)
        if p is not None:
            return p
    pytest.skip("No fissile HDF5 file found; set OPENMC_CROSS_SECTIONS")


@pytest.fixture
def photon_h5_path():
    """Path to a photon HDF5 file for testing."""
    for name in ["Fe.h5", "Li.h5", "Be.h5", "H.h5"]:
        p = _find_hdf5_file(name)
        if p is not None:
            return p
    pytest.skip("No photon HDF5 file found; set OPENMC_CROSS_SECTIONS")
