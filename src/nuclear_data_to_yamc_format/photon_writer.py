"""Export an OpenMC IncidentPhoton object to an Arrow IPC directory."""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.feather as pf

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


def export_photon_to_arrow(data, path):
    """Export an IncidentPhoton object to an Arrow IPC directory.

    Parameters
    ----------
    data : openmc.data.IncidentPhoton
        The incident photon data to export.
    path : str or Path
        Directory path to write Arrow files to (e.g., "H.photon.arrow").
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    Z = data.atomic_number

    # Build union energy grid (same as export_to_hdf5)
    union_grid = np.array([])
    for rx in data:
        union_grid = np.union1d(union_grid, rx.xs.x)

    # ------------------------------------------------------------------
    # Collect cross sections and form factors from reactions
    # ------------------------------------------------------------------
    xs_data = {}  # key -> xs array on union grid
    coherent_rx = None
    incoherent_rx = None
    subshell_rows = []
    designators = []

    # Map from photon.py _REACTION_NAME
    from openmc.data.photon import _REACTION_NAME

    for mt, rx in data.reactions.items():
        name, key = _REACTION_NAME[mt]

        if mt in (502, 504, 515, 517, 522, 525):
            # Main cross section on union grid
            if mt >= 534 and mt <= 572:
                pass  # handled below
            else:
                xs_data[key] = np.asarray(rx.xs(union_grid), dtype=np.float64)

            if mt == 502:
                coherent_rx = rx
            elif mt == 504:
                incoherent_rx = rx

        elif mt >= 534 and mt <= 572:
            designators.append(key)

            # Determine threshold
            threshold = rx.xs.x[0]
            idx = int(np.searchsorted(union_grid, threshold, side='right') - 1)
            photoionization = np.asarray(rx.xs(union_grid[idx:]), dtype=np.float64)

            # Subshell atomic relaxation data
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
                "threshold_idx": idx,
                "transitions_data": transitions_data,
                "transitions_shape": transitions_shape,
            })

    # ------------------------------------------------------------------
    # element.feather
    # ------------------------------------------------------------------
    element_row = {
        "name": data.name,
        "Z": int(Z),
        "energy": union_grid.tolist(),
    }

    # Main cross sections
    _xs_key_map = {
        "coherent": "coherent_xs",
        "incoherent": "incoherent_xs",
        "photoelectric": "photoelectric_xs",
        "pair_production_nuclear": "pair_production_nuclear_xs",
        "pair_production_electron": "pair_production_electron_xs",
        "heating": "heating_xs",
    }
    for key, col_name in _xs_key_map.items():
        if key in xs_data:
            element_row[col_name] = xs_data[key].tolist()
        else:
            element_row[col_name] = []

    # Form factors for coherent scattering
    if coherent_rx is not None and coherent_rx.scattering_factor is not None:
        ff = coherent_rx.scattering_factor
        # Integrated form factor (same computation as photon.py to_hdf5)
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

    # Incoherent scattering factor
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
    pf.write_feather(element_table, path / "element.feather")

    # ------------------------------------------------------------------
    # subshells.feather
    # ------------------------------------------------------------------
    if subshell_rows:
        subshells_table = pa.table(
            {col: [r[col] for r in subshell_rows]
             for col in SUBSHELLS_SCHEMA.names},
            schema=SUBSHELLS_SCHEMA,
        )
        pf.write_feather(subshells_table, path / "subshells.feather")

    # ------------------------------------------------------------------
    # compton.feather
    # ------------------------------------------------------------------
    if data.compton_profiles:
        profile = data.compton_profiles
        J_arr = np.array([Jk.y for Jk in profile['J']], dtype=np.float64)
        compton_row = {
            "num_electrons": np.asarray(profile['num_electrons'], dtype=np.float64).tolist(),
            "binding_energy": np.asarray(profile['binding_energy'], dtype=np.float64).tolist(),
            "pz": np.asarray(profile['J'][0].x, dtype=np.float64).tolist(),
            "J_data": J_arr.ravel(order='C').tolist(),
            "J_shape": list(J_arr.shape),
        }
        compton_table = pa.table(
            {col: [compton_row[col]] for col in COMPTON_SCHEMA.names},
            schema=COMPTON_SCHEMA,
        )
        pf.write_feather(compton_table, path / "compton.feather")

    # ------------------------------------------------------------------
    # bremsstrahlung.feather
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
        pf.write_feather(brem_table, path / "bremsstrahlung.feather")
