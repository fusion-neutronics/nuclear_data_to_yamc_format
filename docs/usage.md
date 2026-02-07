# Usage

## Python API

The two main functions accept nuclear data source files and write Arrow
output:

```python
from nuclear_data_to_yamc_format import convert_neutron, convert_photon
```

### Neutron data

From an ACE file:

```python
convert_neutron("92235.710nc", "output/")
# creates output/U235.arrow/
```

From an ENDF file (NJOY is invoked automatically):

```python
convert_neutron(
    "n-092_U_235.endf",
    "output/",
    source_format="endf",
    temperatures=[293.6, 600.0, 900.0],
)
```

### Photon data

Photon files are always ENDF.  If a separate atomic relaxation file is
available, pass it with `atom_path`:

```python
convert_photon(
    "photoat-026_Fe_000.endf",
    "output/",
    atom_path="atom-026_Fe_000.endf",
)
# creates output/Fe.photon.arrow/
```

## Command-line script

A convenience script is included for single-file conversion:

```bash
# Neutron from ACE
python scripts/convert_single_file.py neutron 92235.710nc -o output/

# Neutron from ENDF
python scripts/convert_single_file.py neutron n-092_U_235.endf -f endf -o output/

# Photon from ENDF with atomic relaxation
python scripts/convert_single_file.py photon photoat-026_Fe_000.endf \
    --atom atom-026_Fe_000.endf -o output/
```

### Bulk conversion scripts

Two scripts are provided for converting entire libraries:

`scripts/convert_endf.py`
: Downloads and converts the ENDF/B-VIII.0 (or VII.1) library.  Mirrors
  the `openmc_data` `generate_endf.py` pipeline.

`scripts/convert_fendl.py`
: Downloads and converts the FENDL library.  Mirrors the `openmc_data`
  `convert_fendl.py` pipeline.

Both scripts accept `--download` / `--no-download` and `--extract` /
`--no-extract` flags so you can skip steps if the source files are already
present.

## Reading Arrow files back

For verification or inspection, you can read the Arrow files back into
Python dicts:

```python
from nuclear_data_to_yamc_format import read_neutron_from_arrow, read_photon_from_arrow

data = read_neutron_from_arrow("output/U235.arrow")
print(data["nuclide"]["name"])        # "U235"
print(len(data["reactions"]))         # number of reactions
print(data["reactions"][0]["mt"])      # first reaction MT number

photon = read_photon_from_arrow("output/Fe.photon.arrow")
print(photon["element"]["Z"])         # 26
```

## Verification

To confirm that the Arrow output is bit-identical to what OpenMC would
produce, use the verification functions:

```python
import openmc.data
from nuclear_data_to_yamc_format import export_neutron_to_arrow, verify_neutron

data = openmc.data.IncidentNeutron.from_ace("92235.710nc")
export_neutron_to_arrow(data, "U235.arrow")
assert verify_neutron(data, "U235.arrow")
```
