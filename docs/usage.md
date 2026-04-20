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

For ENDF files containing multiple elements (e.g. FENDL bundles):

```python
from nuclear_data_to_yamc_format import convert_photon_endf

convert_photon_endf(
    input_path="fendl-photoat.endf",
    output_dir="output/",
    library="fendl-3.2c",
)
# creates output/H.arrow/, output/He.arrow/, ...
```

## Command-line tools

Installing the package registers four console entry points.

### `convert-single-file`

Convert a single ACE or ENDF file:

```bash
# Neutron from ACE
convert-single-file neutron 92235.710nc -o output/

# Neutron from ACE with library name
convert-single-file neutron 92235.710nc -o output/ --library endfb-8.0

# Neutron from ENDF
convert-single-file neutron n-092_U_235.endf -f endf -o output/

# Photon from ENDF with atomic relaxation
convert-single-file photon photoat-026_Fe_000.endf \
    --atom atom-026_Fe_000.endf -o output/ --library endfb-8.0
```

### Bulk conversion

`convert-endf`
: Converts the ENDF/B-VIII.1 (or VIII.0, VII.1) library using NJOY for Doppler
  broadening.  Automatically downloads ENDF files from NNDC if not found
  locally (searches `./endfb-{release}-endf/` then
  `~/nuclear_data/endfb-{release}-endf/`).  Processes neutrons in parallel.

  ```bash
  # Default: ENDF/B-VIII.1, 6 temperatures, output to ~/nuclear_data/endf-b8.1-arrow/
  convert-endf

  # ENDF/B-VII.1 with custom temperatures
  convert-endf -r vii.1 --temperatures 293.6 600.0
  ```

`convert-fendl`
: Converts the FENDL library (ACE neutrons + ENDF photons).
  Automatically downloads from IAEA if not found locally (searches
  `./fendl-{release}-ace/` and `./fendl-{release}-endf/` then
  `~/nuclear_data/fendl-{release}-{ace,endf}/`).

  ```bash
  # Default: FENDL 3.2c, output to ./fendl-3.2c-arrow/
  convert-fendl

  # FENDL 3.1d, custom output
  convert-fendl -r 3.1d -d /path/to/output
  ```

`convert-tendl`
: Converts the TENDL library (neutron only). Pipeline: download ENDF → NJOY →
  OpenMC → Arrow.

  ```bash
  # Default: TENDL 2025
  convert-tendl

  # TENDL 2023, filter to specific nuclides
  convert-tendl -r 2023 --nuclides Fe56 U235
  ```

All bulk commands accept `--cleanup` to remove source files after conversion
and `--nuclides` to filter to a subset of isotopes.

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
