"""Read Arrow IPC neutron data back into Python dicts for verification."""

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


def read_neutron_from_arrow(path):
    """Read an Arrow neutron directory into a dict of dicts/lists.

    Parameters
    ----------
    path : str or Path
        Path to the Arrow directory (e.g., "Li6.arrow").

    Returns
    -------
    dict
        Dictionary with keys: "nuclide", "reactions", "products",
        "distributions", "urr" (optional), "total_nu" (optional).
    """
    path = Path(path)
    result = {}

    # nuclide
    nuclide_table = pf.read_table(path / "nuclide.feather")
    result["nuclide"] = _read_single_row(nuclide_table)
    result["nuclide"]["_metadata"] = nuclide_table.schema.metadata

    # reactions
    if (path / "reactions.feather").exists():
        reactions_table = pf.read_table(path / "reactions.feather")
        result["reactions"] = []
        for i in range(len(reactions_table)):
            row = {}
            for name in reactions_table.schema.names:
                row[name] = _to_list(reactions_table.column(name)[i])
            result["reactions"].append(row)
    else:
        result["reactions"] = []

    # products
    if (path / "products.feather").exists():
        products_table = pf.read_table(path / "products.feather")
        result["products"] = []
        for i in range(len(products_table)):
            row = {}
            for name in products_table.schema.names:
                row[name] = _to_list(products_table.column(name)[i])
            result["products"].append(row)
    else:
        result["products"] = []

    # distributions
    if (path / "distributions.feather").exists():
        distributions_table = pf.read_table(path / "distributions.feather")
        result["distributions"] = []
        for i in range(len(distributions_table)):
            row = {}
            for name in distributions_table.schema.names:
                row[name] = _to_list(distributions_table.column(name)[i])
            result["distributions"].append(row)
    else:
        result["distributions"] = []

    # urr (optional)
    if (path / "urr.feather").exists():
        urr_table = pf.read_table(path / "urr.feather")
        result["urr"] = []
        for i in range(len(urr_table)):
            row = {}
            for name in urr_table.schema.names:
                row[name] = _to_list(urr_table.column(name)[i])
            result["urr"].append(row)

    # total_nu (optional)
    if (path / "total_nu.feather").exists():
        total_nu_table = pf.read_table(path / "total_nu.feather")
        result["total_nu"] = _read_single_row(total_nu_table)

    return result
