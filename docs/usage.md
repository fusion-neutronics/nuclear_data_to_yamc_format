# Usage

## Python API

The two main functions accept nuclear data source files and write
simulation-ready Arrow output:

```python
from nuclear_data_to_yamc_format import convert_neutron, convert_photon
```

### Neutron data

From an ACE file:

```python
convert_neutron(input_path="92235.710nc", output_dir="output/")
# creates output/U235.arrow/
```

From an ENDF file (NJOY is invoked automatically):

```python
convert_neutron(
    input_path="n-092_U_235.endf",
    output_dir="output/",
    source_format="endf",
    temperatures=[293.6, 600.0, 900.0],
)
```

With a library name (stored in version.json):

```python
convert_neutron(input_path="92235.710nc", output_dir="output/", library="endfb-8.0")
```

### Photon data

Photon files are always ENDF.  If a separate atomic relaxation file is
available, pass it with `atom_path`:

```python
convert_photon(
    input_path="photoat-026_Fe_000.endf",
    output_dir="output/",
    atom_path="atom-026_Fe_000.endf",
    library="endfb-8.0",
)
# creates output/Fe.arrow/
```

## Command-line script

A convenience script is included for single-file conversion:

```bash
# Neutron from ACE
python scripts/convert_single_file.py neutron 92235.710nc -o output/

# Neutron from ACE with library name
python scripts/convert_single_file.py neutron 92235.710nc -o output/ --library endfb-8.0

# Neutron from ENDF
python scripts/convert_single_file.py neutron n-092_U_235.endf -f endf -o output/

# Photon from ENDF with atomic relaxation
python scripts/convert_single_file.py photon photoat-026_Fe_000.endf \
    --atom atom-026_Fe_000.endf -o output/ --library endfb-8.0
```

### Bulk conversion scripts

Two scripts are provided for converting entire libraries:

`scripts/convert_endf.py`
: Downloads and converts the ENDF/B-VIII.0 (or VII.1) library.  Mirrors
  the `openmc_data` `generate_endf.py` pipeline.  Automatically sets
  the library name ("endfb-8.0" or "endfb-7.1").

`scripts/convert_fendl.py`
: Downloads and converts the FENDL library.  Mirrors the `openmc_data`
  `convert_fendl.py` pipeline.  Automatically sets the library name
  (e.g., "fendl-3.2c").

Both scripts accept `--download` / `--no-download` and `--extract` /
`--no-extract` flags so you can skip steps if the source files are already
present.

## Reading Arrow files back

For verification or inspection, you can read the Arrow files back into
Python dicts:

```python
from nuclear_data_to_yamc_format import read_neutron_from_arrow, read_photon_from_arrow

data = read_neutron_from_arrow(path="output/U235.arrow")
print(data["version"])                   # format and converter metadata
print(data["nuclide"]["name"])           # "U235"
print(len(data["reactions"]))            # number of reactions (incl. synthesized)
print(data["reactions"][0]["mt"])         # first reaction MT number

# FastXSGrid data
for fxs in data["fast_xs"]:
    print(fxs["temperature"], len(fxs["log_grid_index"]))  # 8001 entries

photon = read_photon_from_arrow(path="output/Fe.arrow")
print(photon["element"]["Z"])            # 26
print(len(photon["element"]["ln_energy"]))  # log-space energy grid
```

## Verification

To confirm that the Arrow output matches OpenMC data within tolerance,
use the verification functions:

```python
import openmc.data
from nuclear_data_to_yamc_format import export_neutron_to_arrow, verify_neutron

data = openmc.data.IncidentNeutron.from_ace("92235.710nc")
export_neutron_to_arrow(data=data, path="U235.arrow")
assert verify_neutron(data=data, arrow_path="U235.arrow")
```

Verification uses `np.allclose` with `rtol=1e-12` for floating-point
comparisons and additionally checks:

- Synthesized MTs (1, 3, 4, 27, 101) are present and marked redundant
- FastXSGrid: `total == absorption + scattering + fission`
- `log_grid_index` is monotonically non-decreasing with 8001 entries
- `version.json` is present with valid `format_version`
