"""Export an OpenMC IncidentPhoton object to a simulation-ready .arrow/ directory.

Stores energy and cross sections in both linear and log space for fast
interpolation.  Pre-computes Compton profile CDFs via trapezoidal integration.
"""

import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc

from openmc.data.function import Tabulated1D

from .schemas import (
    ELEMENT_SCHEMA,
    SUBSHELLS_SCHEMA,
    COMPTON_SCHEMA,
    BREMSSTRAHLUNG_SCHEMA,
)

# Map MT numbers to (name, key) — mirrors _REACTION_NAME in photon.py
_MT_KEY_MAP = {
    502: ("coherent", "coherent"),
    504: ("incoherent", "incoherent"),
    515: ("pair_production_electron", "pair_production_electron"),
    517: ("pair_production_nuclear", "pair_production_nuclear"),
    522: ("photoelectric", "photoelectric"),
    525: ("heating", "heating"),
}


def _safe_log(arr):
    """Compute log, replacing zeros/negatives with a very small value."""
    arr = np.asarray(arr, dtype=np.float64)
    safe = np.where(arr > 0, arr, 1e-300)
    return np.log(safe)


def _compute_compton_cdfs(J_arr):
    """Compute CDFs from Compton profile J values via trapezoidal integration.

    Parameters
    ----------
    J_arr : numpy.ndarray
        Shape (n_shells, n_pz).

    Returns
    -------
    numpy.ndarray
        CDF array of same shape, normalized so last value is 1.0.
    """
    n_shells, n_pz = J_arr.shape
    cdf = np.zeros_like(J_arr)
    for s in range(n_shells):
        # Trapezoidal integration (cumulative)
        for i in range(1, n_pz):
            cdf[s, i] = cdf[s, i-1] + 0.5 * (J_arr[s, i-1] + J_arr[s, i])
        # Normalize
        if cdf[s, -1] > 0:
            cdf[s, :] /= cdf[s, -1]
    return cdf


def _write_arrow_ipc(table, filepath):
    """Write a PyArrow table as Arrow IPC (not Feather)."""
    with pa.OSFile(str(filepath), 'wb') as f:
        writer = ipc.new_file(f, table.schema)
        writer.write_table(table)
        writer.close()


def export_photon_to_arrow(data, path, *, library=""):
    """Export an IncidentPhoton object to a simulation-ready .arrow/ directory.

    Parameters
    ----------
    data : openmc.data.IncidentPhoton
        The incident photon data to export.
    path : str or Path
        Directory path to write Arrow files to (e.g., "Fe.arrow").
    library : str, optional
        Library name (e.g., "endfb-8.0", "fendl-3.2c").
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    Z = data.atomic_number

    # ------------------------------------------------------------------
    # 0. version.json
    # ------------------------------------------------------------------
    from . import __version__
    version_info = {
        "format_version": 1,
        "library": library,
        "converter_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (path / "version.json").write_text(json.dumps(version_info, indent=2))

    # Build union energy grid (same as export_to_hdf5)
    union_grid = np.array([])
    for rx in data:
        union_grid = np.union1d(union_grid, rx.xs.x)

    ln_energy = _safe_log(union_grid)

    # ------------------------------------------------------------------
    # Collect cross sections and form factors from reactions
    # ------------------------------------------------------------------
    xs_data = {}  # key -> xs array on union grid
    coherent_rx = None
    incoherent_rx = None
    subshell_rows = []

    from openmc.data.photon import _REACTION_NAME

    for mt, rx in data.reactions.items():
        name, key = _REACTION_NAME[mt]

        if mt in (502, 504, 515, 517, 522, 525):
            if mt >= 534 and mt <= 572:
                pass  # handled below
            else:
                xs_data[key] = np.asarray(rx.xs(union_grid), dtype=np.float64)

            if mt == 502:
                coherent_rx = rx
            elif mt == 504:
                incoherent_rx = rx

        elif mt >= 534 and mt <= 572:
            threshold = rx.xs.x[0]
            idx = int(np.searchsorted(union_grid, threshold, side='right') - 1)
            photoionization = np.asarray(rx.xs(union_grid[idx:]), dtype=np.float64)
            ln_photoionization = _safe_log(photoionization)

            binding_energy = 0.0
            num_electrons = 0.0
            transitions_data = None
            transitions_shape = None

            if data.atomic_relaxation is not None:
                if key in data.atomic_relaxation.subshells:
                    ar = data.atomic_relaxation
                    binding_energy = float(ar.binding_energy[key])
                    num_electrons = float(ar.num_electrons[key])
                    if key in ar.transitions:
                        from openmc.data.photon import _SUBSHELLS
                        import pandas as pd
                        with pd.option_context('future.no_silent_downcasting', True):
                            df = ar.transitions[key].replace(
                                _SUBSHELLS, range(len(_SUBSHELLS)))
                        t_arr = df.values.astype(float)
                        transitions_data = t_arr.ravel(order='C').tolist()
                        transitions_shape = list(t_arr.shape)

            subshell_rows.append({
                "designator": key,
                "binding_energy": binding_energy,
                "num_electrons": num_electrons,
                "xs": photoionization.tolist(),
                "ln_xs": ln_photoionization.tolist(),
                "threshold_idx": idx,
                "transitions_data": transitions_data,
                "transitions_shape": transitions_shape,
            })

    # ------------------------------------------------------------------
    # element.arrow
    # ------------------------------------------------------------------
    element_row = {
        "name": data.name,
        "Z": int(Z),
        "energy": union_grid.tolist(),
        "ln_energy": ln_energy.tolist(),
    }

    _xs_key_map = {
        "coherent": "coherent_xs",
        "incoherent": "incoherent_xs",
        "photoelectric": "photoelectric_xs",
        "pair_production_nuclear": "pair_production_nuclear_xs",
        "pair_production_electron": "pair_production_electron_xs",
        "heating": "heating_xs",
    }
    _ln_xs_key_map = {
        "coherent": "ln_coherent_xs",
        "incoherent": "ln_incoherent_xs",
        "photoelectric": "ln_photoelectric_xs",
    }
    for key, col_name in _xs_key_map.items():
        if key in xs_data:
            element_row[col_name] = xs_data[key].tolist()
        else:
            element_row[col_name] = []
    for key, col_name in _ln_xs_key_map.items():
        if key in xs_data:
            element_row[col_name] = _safe_log(xs_data[key]).tolist()
        else:
            element_row[col_name] = []

    # Form factors for coherent scattering
    if coherent_rx is not None and coherent_rx.scattering_factor is not None:
        ff = coherent_rx.scattering_factor
        ff_copy = deepcopy(ff)
        ff_copy.x = ff_copy.x * ff_copy.x
        ff_copy.y = ff_copy.y * ff_copy.y / Z**2
        int_ff = Tabulated1D(ff_copy.x, ff_copy.integral())
        element_row["coherent_int_ff_x"] = np.asarray(int_ff.x, dtype=np.float64).tolist()
        element_row["coherent_int_ff_y"] = np.asarray(int_ff.y, dtype=np.float64).tolist()
        element_row["coherent_ff_x"] = np.asarray(ff.x, dtype=np.float64).tolist()
        element_row["coherent_ff_y"] = np.asarray(ff.y, dtype=np.float64).tolist()
    else:
        element_row["coherent_int_ff_x"] = []
        element_row["coherent_int_ff_y"] = []
        element_row["coherent_ff_x"] = []
        element_row["coherent_ff_y"] = []

    if coherent_rx is not None and coherent_rx.anomalous_real is not None:
        element_row["coherent_anomalous_real_x"] = np.asarray(
            coherent_rx.anomalous_real.x, dtype=np.float64).tolist()
        element_row["coherent_anomalous_real_y"] = np.asarray(
            coherent_rx.anomalous_real.y, dtype=np.float64).tolist()
    else:
        element_row["coherent_anomalous_real_x"] = []
        element_row["coherent_anomalous_real_y"] = []

    if coherent_rx is not None and coherent_rx.anomalous_imag is not None:
        element_row["coherent_anomalous_imag_x"] = np.asarray(
            coherent_rx.anomalous_imag.x, dtype=np.float64).tolist()
        element_row["coherent_anomalous_imag_y"] = np.asarray(
            coherent_rx.anomalous_imag.y, dtype=np.float64).tolist()
    else:
        element_row["coherent_anomalous_imag_x"] = []
        element_row["coherent_anomalous_imag_y"] = []

    if incoherent_rx is not None and incoherent_rx.scattering_factor is not None:
        element_row["incoherent_ff_x"] = np.asarray(
            incoherent_rx.scattering_factor.x, dtype=np.float64).tolist()
        element_row["incoherent_ff_y"] = np.asarray(
            incoherent_rx.scattering_factor.y, dtype=np.float64).tolist()
    else:
        element_row["incoherent_ff_x"] = []
        element_row["incoherent_ff_y"] = []

    element_table = pa.table(
        {col: [element_row[col]] for col in ELEMENT_SCHEMA.names},
        schema=ELEMENT_SCHEMA,
    )
    _write_arrow_ipc(element_table, path / "element.arrow")

    # ------------------------------------------------------------------
    # subshells.arrow
    # ------------------------------------------------------------------
    if subshell_rows:
        subshells_table = pa.table(
            {col: [r[col] for r in subshell_rows]
             for col in SUBSHELLS_SCHEMA.names},
            schema=SUBSHELLS_SCHEMA,
        )
        _write_arrow_ipc(subshells_table, path / "subshells.arrow")

    # ------------------------------------------------------------------
    # compton.arrow
    # ------------------------------------------------------------------
    if data.compton_profiles:
        profile = data.compton_profiles
        J_arr = np.array([Jk.y for Jk in profile['J']], dtype=np.float64)

        # Pre-compute CDFs
        J_cdf = _compute_compton_cdfs(J_arr)

        compton_row = {
            "num_electrons": np.asarray(profile['num_electrons'], dtype=np.float64).tolist(),
            "binding_energy": np.asarray(profile['binding_energy'], dtype=np.float64).tolist(),
            "pz": np.asarray(profile['J'][0].x, dtype=np.float64).tolist(),
            "J_data": J_arr.ravel(order='C').tolist(),
            "J_shape": list(J_arr.shape),
            "J_cdf_data": J_cdf.ravel(order='C').tolist(),
            "J_cdf_shape": list(J_cdf.shape),
        }
        compton_table = pa.table(
            {col: [compton_row[col]] for col in COMPTON_SCHEMA.names},
            schema=COMPTON_SCHEMA,
        )
        _write_arrow_ipc(compton_table, path / "compton.arrow")

    # ------------------------------------------------------------------
    # bremsstrahlung.arrow
    # ------------------------------------------------------------------
    if data.bremsstrahlung:
        brem = data.bremsstrahlung
        dcs_arr = np.asarray(brem['dcs'], dtype=np.float64)
        brem_row = {
            "I": float(brem['I']),
            "electron_energy": np.asarray(brem['electron_energy'], dtype=np.float64).tolist(),
            "photon_energy": np.asarray(brem['photon_energy'], dtype=np.float64).tolist(),
            "num_electrons": np.asarray(brem['num_electrons'], dtype=np.float64).tolist(),
            "ionization_energy": np.asarray(brem['ionization_energy'], dtype=np.float64).tolist(),
            "dcs_data": dcs_arr.ravel(order='C').tolist(),
            "dcs_shape": list(dcs_arr.shape),
        }
        brem_table = pa.table(
            {col: [brem_row[col]] for col in BREMSSTRAHLUNG_SCHEMA.names},
            schema=BREMSSTRAHLUNG_SCHEMA,
        )
        _write_arrow_ipc(brem_table, path / "bremsstrahlung.arrow")
