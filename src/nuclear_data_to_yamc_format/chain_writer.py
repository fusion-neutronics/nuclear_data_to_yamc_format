"""Export an OpenMC depletion Chain object to a simulation-ready .chain.arrow/ directory.

Writes a directory of Arrow IPC tables describing the depletion network:
nuclides, decay modes, reactions, decay photon/electron sources, and
fission product yields.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc

from openmc.stats.univariate import Discrete, Tabular, Mixture

from .schemas import (
    CHAIN_NUCLIDES_SCHEMA,
    CHAIN_DECAYS_SCHEMA,
    CHAIN_REACTIONS_SCHEMA,
    CHAIN_SOURCES_SCHEMA,
    CHAIN_FISSION_YIELDS_SCHEMA,
)


def _write_arrow_ipc(table, filepath):
    with pa.OSFile(str(filepath), 'wb') as f:
        writer = ipc.new_file(f, table.schema)
        writer.write_table(table)
        writer.close()


def _source_rows(nuclide_name, particle, dist):
    """Yield (type, energies, intensities) tuples from a source distribution.

    Mixture distributions are flattened into one row per component, with the
    component probability multiplied into its intensities.
    """
    if isinstance(dist, Discrete):
        yield (
            "discrete",
            np.asarray(dist.x, dtype=np.float64).tolist(),
            np.asarray(dist.p, dtype=np.float64).tolist(),
        )
    elif isinstance(dist, Tabular):
        yield (
            "tabular",
            np.asarray(dist.x, dtype=np.float64).tolist(),
            np.asarray(dist.p, dtype=np.float64).tolist(),
        )
    elif isinstance(dist, Mixture):
        for prob, sub in zip(dist.probability, dist.distribution):
            for sub_type, energies, intensities in _source_rows(
                nuclide_name, particle, sub
            ):
                scaled = (np.asarray(intensities, dtype=np.float64) * float(prob)).tolist()
                yield (sub_type, energies, scaled)
    else:
        raise NotImplementedError(
            f"{nuclide_name} source ({particle}): unsupported distribution "
            f"type {type(dist).__name__}"
        )


def export_chain_to_arrow(chain, path, *, library=""):
    """Export an openmc.deplete.Chain to a .chain.arrow/ directory.

    Parameters
    ----------
    chain : openmc.deplete.Chain
        The depletion chain to export.
    path : str or Path
        Output directory (e.g., "chain_endf_b8.0.chain.arrow").
    library : str, optional
        Library name (e.g., "endfb-8.0").
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    from . import __version__
    version_info = {
        "format_version": 1,
        "library": library,
        "converter_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (path / "version.json").write_text(json.dumps(version_info, indent=2))

    nuclide_rows = []
    decay_rows = []
    reaction_rows = []
    source_rows = []
    fy_rows = []

    for nuc in chain.nuclides:
        fy_parent = getattr(nuc, "_fpy", None)
        has_fy = nuc.yield_data is not None
        n_sources = sum(1 for _ in nuc.sources) if nuc.sources else 0

        nuclide_rows.append({
            "name": nuc.name,
            "half_life": float(nuc.half_life) if nuc.half_life is not None else None,
            "decay_energy": float(nuc.decay_energy),
            "n_decay_modes": int(len(nuc.decay_modes)),
            "n_reactions": int(len(nuc.reactions)),
            "n_sources": int(n_sources),
            "has_fission_yields": bool(has_fy),
            "fission_yield_parent": fy_parent,
        })

        for d in nuc.decay_modes:
            decay_rows.append({
                "nuclide": nuc.name,
                "type": d.type,
                "target": d.target,
                "branching_ratio": float(d.branching_ratio),
            })

        for r in nuc.reactions:
            reaction_rows.append({
                "nuclide": nuc.name,
                "type": r.type,
                "target": r.target,
                "Q": float(r.Q),
                "branching_ratio": float(r.branching_ratio),
            })

        for particle, dist in nuc.sources.items():
            for src_type, energies, intensities in _source_rows(nuc.name, particle, dist):
                source_rows.append({
                    "nuclide": nuc.name,
                    "particle": particle,
                    "type": src_type,
                    "energies": energies,
                    "intensities": intensities,
                })

        if has_fy and fy_parent is None:
            for energy, fy in nuc.yield_data.items():
                products = list(fy.products)
                yields = np.asarray(fy.yields, dtype=np.float64).tolist()
                fy_rows.append({
                    "nuclide": nuc.name,
                    "energy": float(energy),
                    "products": products,
                    "yields": yields,
                })

    def _table(rows, schema):
        return pa.table(
            {col: [r[col] for r in rows] for col in schema.names},
            schema=schema,
        )

    nuclides_table = _table(nuclide_rows, CHAIN_NUCLIDES_SCHEMA)
    _write_arrow_ipc(nuclides_table, path / "nuclides.arrow")

    if decay_rows:
        _write_arrow_ipc(_table(decay_rows, CHAIN_DECAYS_SCHEMA), path / "decays.arrow")

    if reaction_rows:
        _write_arrow_ipc(_table(reaction_rows, CHAIN_REACTIONS_SCHEMA), path / "reactions.arrow")

    if source_rows:
        _write_arrow_ipc(_table(source_rows, CHAIN_SOURCES_SCHEMA), path / "sources.arrow")

    if fy_rows:
        _write_arrow_ipc(_table(fy_rows, CHAIN_FISSION_YIELDS_SCHEMA), path / "fission_yields.arrow")
