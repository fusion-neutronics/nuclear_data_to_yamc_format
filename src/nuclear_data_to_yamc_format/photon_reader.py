"""Read Arrow IPC photon data back into Python dicts for verification."""

from pathlib import Path

import numpy as np
import pyarrow.feather as pf


def _to_list(val):
    """Convert pyarrow scalar/list to Python list or value."""
    if val is None:
        return None
    v = val.as_py() if hasattr(val, 'as_py') else val
    return v


def _read_single_row(table):
    """Read a single-row table into a dict."""
    assert len(table) == 1
    row = {}
    for name in table.schema.names:
        row[name] = _to_list(table.column(name)[0])
    return row


def read_photon_from_arrow(path):
    """Read an Arrow photon directory into a dict of dicts/lists.

    Parameters
    ----------
    path : str or Path
        Path to the Arrow directory (e.g., "H.photon.arrow").

    Returns
    -------
    dict
        Dictionary with keys: "element", "subshells", "compton" (optional),
        "bremsstrahlung" (optional).
    """
    path = Path(path)
    result = {}

    # element
    element_table = pf.read_table(path / "element.feather")
    result["element"] = _read_single_row(element_table)
    result["element"]["_metadata"] = element_table.schema.metadata

    # subshells
    if (path / "subshells.feather").exists():
        subshells_table = pf.read_table(path / "subshells.feather")
        result["subshells"] = []
        for i in range(len(subshells_table)):
            row = {}
            for name in subshells_table.schema.names:
                row[name] = _to_list(subshells_table.column(name)[i])
            result["subshells"].append(row)
    else:
        result["subshells"] = []

    # compton (optional)
    if (path / "compton.feather").exists():
        compton_table = pf.read_table(path / "compton.feather")
        result["compton"] = _read_single_row(compton_table)

    # bremsstrahlung (optional)
    if (path / "bremsstrahlung.feather").exists():
        brem_table = pf.read_table(path / "bremsstrahlung.feather")
        result["bremsstrahlung"] = _read_single_row(brem_table)

    return result
