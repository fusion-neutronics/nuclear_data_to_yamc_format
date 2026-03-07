# nuclear_data_to_yamc_format

**nuclear_data_to_yamc_format** converts nuclear data files (ACE or ENDF) into
simulation-ready [Apache Arrow IPC](https://arrow.apache.org/docs/format/Columnar.html#ipc-file-format)
directories (`.arrow/`).  The output includes pre-computed hierarchical MT
cross sections, FastXSGrid lookup tables, and log-space photon data — everything
yamc needs to start a simulation with near-zero load time.

```text
ACE  ──────────────────► OpenMC objects ──► .arrow/
                              │
ENDF ──► NJOY ──────────►    │    └── synthesize MTs (1,3,4,27,101)
                                  └── build FastXSGrid (8000 log bins)
                                  └── log-space photon XS + Compton CDFs
```

## Quick start

```python
from nuclear_data_to_yamc_format import convert_neutron, convert_photon

# Neutron from ACE
convert_neutron(input_path="92235.710nc", output_dir="output/", library="endfb-8.0")

# Neutron from ENDF (runs NJOY internally)
convert_neutron(
    input_path="n-092_U_235.endf",
    output_dir="output/",
    source_format="endf",
    library="endfb-8.0",
)

# Photon from ENDF
convert_photon(
    input_path="photoat-026_Fe_000.endf",
    output_dir="output/",
    atom_path="atom-026_Fe_000.endf",
    library="endfb-8.0",
)
```

## Output structure

### Neutron: `Fe56.arrow/`

```text
Fe56.arrow/
├── version.json          ← format version, library name, timestamps
├── nuclide.arrow         ← metadata, energy grids
├── reactions.arrow       ← all MTs incl. synthesized 1,3,4,27,101
├── products.arrow        ← reaction products
├── distributions.arrow   ← angle-energy distributions
├── fast_xs.arrow         ← FastXSGrid lookup tables per temperature
├── urr.arrow             (optional: unresolved resonance data)
└── total_nu.arrow        (optional: fission neutron yield)
```

### Photon: `Fe.arrow/`

```text
Fe.arrow/
├── version.json
├── element.arrow         ← metadata, energy + ln(energy), XS + ln(XS)
├── subshells.arrow       ← photoionization per subshell with ln(XS)
├── compton.arrow         (optional: profiles + pre-computed CDFs)
└── bremsstrahlung.arrow  (optional)
```

```{toctree}
:maxdepth: 2
:caption: Contents

getting_started
usage
format
synthesis
api
```
