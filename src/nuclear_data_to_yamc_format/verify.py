"""Compare OpenMC objects against Arrow readback for verification."""

from pathlib import Path

import numpy as np

from openmc.data.function import Tabulated1D, Polynomial
from openmc.data.reaction import REACTION_NAME

from .neutron_reader import read_neutron_from_arrow
from .photon_reader import read_photon_from_arrow


def _arrays_equal(a, b, label=""):
    """Compare two arrays or lists for exact equality."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        print(f"  MISMATCH {label}: shape {a.shape} != {b.shape}")
        return False
    if not np.array_equal(a, b):
        diff_mask = a != b
        n_diff = diff_mask.sum()
        max_diff = np.max(np.abs(a[diff_mask] - b[diff_mask])) if n_diff > 0 else 0
        print(f"  MISMATCH {label}: {n_diff} differing values, max diff = {max_diff}")
        return False
    return True


def _check_scalar(expected, actual, label=""):
    """Check scalar equality."""
    if expected != actual:
        print(f"  MISMATCH {label}: {expected!r} != {actual!r}")
        return False
    return True


def verify_neutron(data, arrow_path):
    """Verify that Arrow export matches the OpenMC IncidentNeutron object.

    Parameters
    ----------
    data : openmc.data.IncidentNeutron
        Original data object.
    arrow_path : str or Path
        Path to the Arrow directory.

    Returns
    -------
    bool
        True if all data matches exactly.
    """
    arrow_path = Path(arrow_path)
    arrow_data = read_neutron_from_arrow(arrow_path)
    ok = True

    # --- Nuclide metadata ---
    nuc = arrow_data["nuclide"]
    ok &= _check_scalar(data.name, nuc["name"], "name")
    ok &= _check_scalar(int(data.atomic_number), nuc["Z"], "Z")
    ok &= _check_scalar(int(data.mass_number), nuc["A"], "A")
    ok &= _check_scalar(int(data.metastable), nuc["metastable"], "metastable")
    ok &= _check_scalar(
        float(data.atomic_weight_ratio), nuc["atomic_weight_ratio"], "AWR")

    # Temperatures
    ok &= _check_scalar(
        list(data.temperatures), nuc["temperatures"], "temperatures")

    # kTs
    expected_kts = [float(data.kTs[i]) for i in range(len(data.temperatures))]
    ok &= _arrays_equal(expected_kts, nuc["kTs"], "kTs")

    # Energy grids
    energy_temps = nuc["energy_temperatures"]
    energy_vals = nuc["energy_values"]
    for i, t in enumerate(energy_temps):
        ok &= _arrays_equal(
            data.energy[t], energy_vals[i], f"energy[{t}]")

    # --- Reactions ---
    expected_rxs = list(data.reactions.values())

    arrow_rxs = arrow_data["reactions"]
    ok &= _check_scalar(len(expected_rxs), len(arrow_rxs), "n_reactions")

    for rx, arx in zip(expected_rxs, arrow_rxs):
        mt = rx.mt
        ok &= _check_scalar(mt, arx["mt"], f"rx[{mt}].mt")
        ok &= _check_scalar(float(rx.q_value), arx["Q_value"], f"rx[{mt}].Q_value")
        ok &= _check_scalar(
            bool(rx.center_of_mass), arx["center_of_mass"],
            f"rx[{mt}].center_of_mass")

        # Cross sections
        for t_idx, T in enumerate(arx["xs_temperatures"]):
            if rx.xs[T] is not None:
                ok &= _arrays_equal(
                    rx.xs[T].y, arx["xs_values"][t_idx],
                    f"rx[{mt}].xs[{T}]")
                expected_thresh = getattr(rx.xs[T], '_threshold_idx', 0)
                ok &= _check_scalar(
                    int(expected_thresh), arx["xs_threshold_idx"][t_idx],
                    f"rx[{mt}].xs_threshold[{T}]")

    # --- Products ---
    arrow_prods = arrow_data["products"]
    prod_idx = 0
    for rx in expected_rxs:
        for p_idx, p in enumerate(rx.products):
            if prod_idx >= len(arrow_prods):
                print(f"  MISMATCH: ran out of arrow products at rx {rx.mt} p {p_idx}")
                ok = False
                continue
            ap = arrow_prods[prod_idx]
            ok &= _check_scalar(
                str(p.particle), ap["particle"],
                f"product[{rx.mt},{p_idx}].particle")
            ok &= _check_scalar(
                str(p.emission_mode), ap["emission_mode"],
                f"product[{rx.mt},{p_idx}].emission_mode")
            ok &= _check_scalar(
                float(p.decay_rate), ap["decay_rate"],
                f"product[{rx.mt},{p_idx}].decay_rate")

            # Verify yield
            if isinstance(p.yield_, Tabulated1D):
                expected_data = np.concatenate([p.yield_.x, p.yield_.y])
                ok &= _arrays_equal(
                    expected_data, ap["yield_data"],
                    f"product[{rx.mt},{p_idx}].yield_data")
            elif isinstance(p.yield_, Polynomial):
                ok &= _arrays_equal(
                    p.yield_.coef, ap["yield_data"],
                    f"product[{rx.mt},{p_idx}].yield_data")

            prod_idx += 1

    # --- Distributions (spot check types) ---
    arrow_dists = arrow_data["distributions"]
    dist_idx = 0
    for rx in expected_rxs:
        for p_idx, p in enumerate(rx.products):
            for d_idx, d in enumerate(p.distribution):
                if dist_idx >= len(arrow_dists):
                    print(f"  MISMATCH: ran out of arrow distributions")
                    ok = False
                    continue
                ad = arrow_dists[dist_idx]
                from openmc.data.uncorrelated import UncorrelatedAngleEnergy
                from openmc.data.correlated import CorrelatedAngleEnergy
                from openmc.data.kalbach_mann import KalbachMann
                from openmc.data.nbody import NBodyPhaseSpace

                if isinstance(d, UncorrelatedAngleEnergy):
                    ok &= _check_scalar(
                        "uncorrelated", ad["type"],
                        f"dist[{rx.mt},{p_idx},{d_idx}].type")
                elif isinstance(d, CorrelatedAngleEnergy):
                    ok &= _check_scalar(
                        "correlated", ad["type"],
                        f"dist[{rx.mt},{p_idx},{d_idx}].type")
                elif isinstance(d, KalbachMann):
                    ok &= _check_scalar(
                        "kalbach-mann", ad["type"],
                        f"dist[{rx.mt},{p_idx},{d_idx}].type")
                elif isinstance(d, NBodyPhaseSpace):
                    ok &= _check_scalar(
                        "nbody", ad["type"],
                        f"dist[{rx.mt},{p_idx},{d_idx}].type")
                    ok &= _check_scalar(
                        int(d.n_particles), ad["nbody_n"],
                        f"dist[{rx.mt},{p_idx},{d_idx}].nbody_n")
                    ok &= _check_scalar(
                        float(d.total_mass), ad["nbody_total_mass"],
                        f"dist[{rx.mt},{p_idx},{d_idx}].nbody_total_mass")
                dist_idx += 1

    # --- URR ---
    if data.urr:
        if "urr" not in arrow_data:
            print("  MISMATCH: URR data missing from Arrow")
            ok = False
        else:
            arrow_urrs = arrow_data["urr"]
            for i, (temp, urr) in enumerate(data.urr.items()):
                au = arrow_urrs[i]
                ok &= _check_scalar(temp, au["temperature"], f"urr[{i}].temperature")
                ok &= _arrays_equal(urr.energy, au["energy"], f"urr[{i}].energy")
                table_flat = np.asarray(urr.table).ravel(order='C')
                ok &= _arrays_equal(
                    table_flat, au["table_data"], f"urr[{i}].table_data")
                ok &= _check_scalar(
                    int(urr.interpolation), au["interpolation"],
                    f"urr[{i}].interpolation")

    # --- total_nu ---
    total_nu_written = False
    for rx in data.reactions.values():
        if len(rx.derived_products) > 0 and not total_nu_written:
            dp = rx.derived_products[0]
            total_nu_written = True
            if "total_nu" not in arrow_data:
                print("  MISMATCH: total_nu missing from Arrow")
                ok = False
            else:
                tnu = arrow_data["total_nu"]
                ok &= _check_scalar(
                    str(dp.particle), tnu["particle"], "total_nu.particle")
                if isinstance(dp.yield_, Tabulated1D):
                    expected_data = np.concatenate([dp.yield_.x, dp.yield_.y])
                    ok &= _arrays_equal(
                        expected_data, tnu["yield_data"], "total_nu.yield_data")
                elif isinstance(dp.yield_, Polynomial):
                    ok &= _arrays_equal(
                        dp.yield_.coef, tnu["yield_data"], "total_nu.yield_data")

    if ok:
        print("Neutron verification PASSED")
    else:
        print("Neutron verification FAILED")
    return ok


def verify_photon(data, arrow_path):
    """Verify that Arrow export matches the OpenMC IncidentPhoton object.

    Parameters
    ----------
    data : openmc.data.IncidentPhoton
        Original data object.
    arrow_path : str or Path
        Path to the Arrow directory.

    Returns
    -------
    bool
        True if all data matches exactly.
    """
    arrow_path = Path(arrow_path)
    arrow_data = read_photon_from_arrow(arrow_path)
    ok = True

    elem = arrow_data["element"]
    ok &= _check_scalar(data.name, elem["name"], "name")
    ok &= _check_scalar(int(data.atomic_number), elem["Z"], "Z")

    # Union energy grid
    union_grid = np.array([])
    for rx in data:
        union_grid = np.union1d(union_grid, rx.xs.x)
    ok &= _arrays_equal(union_grid, elem["energy"], "energy")

    # Main cross sections
    from openmc.data.photon import _REACTION_NAME
    for mt, rx in data.reactions.items():
        name, key = _REACTION_NAME[mt]
        if mt in (502, 504, 515, 517, 522, 525):
            expected_xs = np.asarray(rx.xs(union_grid), dtype=np.float64)
            col_map = {
                "coherent": "coherent_xs",
                "incoherent": "incoherent_xs",
                "photoelectric": "photoelectric_xs",
                "pair_production_nuclear": "pair_production_nuclear_xs",
                "pair_production_electron": "pair_production_electron_xs",
                "heating": "heating_xs",
            }
            if key in col_map:
                ok &= _arrays_equal(
                    expected_xs, elem[col_map[key]], f"xs[{key}]")

    # Subshells
    if arrow_data["subshells"]:
        for i, sub in enumerate(arrow_data["subshells"]):
            desig = sub["designator"]
            # Find corresponding reaction
            for mt, rx in data.reactions.items():
                _, rkey = _REACTION_NAME[mt]
                if rkey == desig and mt >= 534 and mt <= 572:
                    threshold = rx.xs.x[0]
                    idx = int(np.searchsorted(union_grid, threshold, side='right') - 1)
                    expected_xs = np.asarray(
                        rx.xs(union_grid[idx:]), dtype=np.float64)
                    ok &= _arrays_equal(
                        expected_xs, sub["xs"], f"subshell[{desig}].xs")
                    ok &= _check_scalar(
                        idx, sub["threshold_idx"],
                        f"subshell[{desig}].threshold_idx")
                    break

    # Compton profiles
    if data.compton_profiles and "compton" in arrow_data:
        profile = data.compton_profiles
        cmp = arrow_data["compton"]
        ok &= _arrays_equal(
            profile['num_electrons'], cmp["num_electrons"], "compton.num_electrons")
        ok &= _arrays_equal(
            profile['binding_energy'], cmp["binding_energy"], "compton.binding_energy")
        ok &= _arrays_equal(
            profile['J'][0].x, cmp["pz"], "compton.pz")
        J_arr = np.array([Jk.y for Jk in profile['J']])
        ok &= _arrays_equal(
            J_arr.ravel(order='C'), cmp["J_data"], "compton.J_data")

    # Bremsstrahlung
    if data.bremsstrahlung and "bremsstrahlung" in arrow_data:
        brem = data.bremsstrahlung
        abrem = arrow_data["bremsstrahlung"]
        ok &= _check_scalar(float(brem['I']), abrem["I"], "brem.I")
        ok &= _arrays_equal(
            brem['electron_energy'], abrem["electron_energy"],
            "brem.electron_energy")
        ok &= _arrays_equal(
            brem['photon_energy'], abrem["photon_energy"],
            "brem.photon_energy")
        dcs_arr = np.asarray(brem['dcs'])
        ok &= _arrays_equal(
            dcs_arr.ravel(order='C'), abrem["dcs_data"], "brem.dcs_data")

    if ok:
        print("Photon verification PASSED")
    else:
        print("Photon verification FAILED")
    return ok
