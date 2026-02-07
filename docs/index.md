# nuclear_data_to_yamc_format

**nuclear_data_to_yamc_format** converts nuclear data files (ACE or ENDF format) into
[Apache Arrow IPC](https://arrow.apache.org/docs/format/Columnar.html#ipc-file-format)
files.  It replaces the HDF5 export step of the standard OpenMC data pipeline
with a set of mmap-friendly Feather v2 files that can be loaded with zero-copy
reads.

```text
ENDF ──► NJOY ──► ACE ──► OpenMC objects ──► Arrow IPC
                              │
                              └── (previously) ──► HDF5
```

```python
from nuclear_data_to_yamc_format import convert_neutron, convert_photon

# Neutron from ACE
convert_neutron("92235.710nc", "output/")

# Neutron from ENDF (runs NJOY internally)
convert_neutron("n-092_U_235.endf", "output/", source_format="endf")

# Photon from ENDF
convert_photon("photoat-026_Fe_000.endf", "output/",
               atom_path="atom-026_Fe_000.endf")
```

```{toctree}
:maxdepth: 2
:caption: Contents

getting_started
usage
format
api
```
