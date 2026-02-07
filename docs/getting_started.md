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
git clone https://github.com/your-org/nuclear_data_to_yamc_format.git
cd nuclear_data_to_yamc_format
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Quick example

Convert a single ACE neutron file:

```python
from nuclear_data_to_yamc_format import convert_neutron

convert_neutron("92235.710nc", "output/")
```

This creates `output/U235.arrow/` containing:

```text
output/U235.arrow/
├── nuclide.feather
├── reactions.feather
├── products.feather
├── distributions.feather
├── urr.feather          (if unresolved resonance data exists)
└── total_nu.feather     (if fission data exists)
```

Each `.feather` file is an Arrow IPC file that can be memory-mapped for
zero-copy reads.
