# Getting started

## Prerequisites

* Python 3.9+
* [OpenMC Python API](https://docs.openmc.org) (must be installed separately)
* NJOY is required when converting from ENDF format (`source_format="endf"`)

## Installation

```bash
pip install nuclear_data_to_yamc_format
```

Or install from source:

```bash
git clone https://github.com/yamc-org/nuclear_data_to_yamc_format.git
cd nuclear_data_to_yamc_format
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

OpenMC must be installed separately (it is not on PyPI):

```bash
pip install --extra-index-url https://shimwell.github.io/wheels openmc
```

## Quick example

Convert a single ACE neutron file:

```python
from nuclear_data_to_yamc_format import convert_neutron

convert_neutron(input_path="92235.710nc", output_dir="output/", library="endfb-8.0")
```

This creates `output/U235.yamc.arrow/` containing:

```text
output/U235.yamc.arrow/
├── version.json
├── nuclide.arrow
├── reactions.arrow       ← includes synthesized MTs 1, 3, 4, 27, 101
├── products.arrow
├── distributions.arrow
├── fast_xs.arrow         ← FastXSGrid lookup tables
├── urr.arrow             (if unresolved resonance data exists)
└── total_nu.arrow        (if fission data exists)
```

Each `.arrow` file is an Arrow IPC file that can be memory-mapped for
zero-copy reads.  The `version.json` file records the format version,
source library, and converter version.

## What's in the output?

The `.yamc.arrow/` format is **simulation-ready**.  Unlike raw Arrow exports,
it includes all the post-processing that yamc would otherwise do at load time:

Hierarchical MT synthesis
: MTs 1 (total), 3 (non-elastic), 4 (inelastic), 27 (absorption), and
  101 (disappearance) are pre-computed from constituent reactions.

FastXSGrid lookup tables
: A log-binned index (8000 bins) over the energy grid enables O(1) cross-section
  retrieval during particle transport.  The 4-column XS array stores total,
  absorption, scattering, and fission cross sections.

Log-space photon data
: Photon energy grids and cross sections are stored in both linear and
  log space.  Compton profile CDFs are pre-computed via trapezoidal
  integration.

## Running the tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

Tests require HDF5 nuclear data files.  Set `OPENMC_CROSS_SECTIONS` to point
at your cross-section library, or place `.h5` files in a known search path.
