# Nuclear Data_to YAMC Format

**nuclear_data_to_yamc_format** converts nuclear data files (ACE or ENDF) into
simulation-ready [Apache Arrow IPC](https://arrow.apache.org/docs/format/Columnar.html#ipc-file-format)
directories (`.arrow/`).  The output is tailored for GPU-ready Monte Carlo
transport: columnar, fixed-width, zero-copy-mappable layouts that stream
directly into device buffers.  It includes pre-computed hierarchical MT
cross sections, FastXSGrid lookup tables, and log-space photon data — everything YAMC needs to start a simulation with near-zero load time.

Depletion and activation workflows are covered by a matching
`transmutation_{library}.arrow/` format: one directory of Arrow tables holding the
full transmutation network (decay data, decay-product sources, transmutation
reactions, and fission product yields).

```{mermaid}
flowchart LR
    ACE[ACE] --> OMC_N[OpenMC IncidentNeutron]
    ENDF_N[ENDF neutron] --> NJOY[NJOY] --> OMC_N
    OMC_N --> NARROW["neutron .arrow/"]
    NARROW --- N1(["synthesize MTs · FastXSGrid"])

    ENDF_P["ENDF photon + atomic relax"] --> OMC_P[OpenMC IncidentPhoton]
    OMC_P --> PARROW["photon .arrow/"]
    PARROW --- N2(["log-space photon XS · Compton CDFs"])

    ENDF_D[ENDF decay] --> CHAIN[OpenMC Chain]
    ENDF_Y[ENDF NFY] --> CHAIN
    ENDF_N --> CHAIN
    XML[Chain XML] --> CHAIN
    JSON[branch ratios JSON] -.-> CHAIN
    CHAIN --> TARROW["transmutation_{library}.arrow/"]
    TARROW --- N3(["nuclides · decays · reactions · decay γ/e⁻ sources · fission yields"])

    style N1 fill:#f5f5f5,stroke:#ccc
    style N2 fill:#f5f5f5,stroke:#ccc
    style N3 fill:#f5f5f5,stroke:#ccc
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

### Transmutation network: `transmutation_{library}.arrow/`

```text
transmutation_endf_b8.1.arrow/
├── version.json
├── nuclides.arrow        ← one row per nuclide; index for the other tables
├── decays.arrow          (optional: decay modes + branching ratios)
├── reactions.arrow       (optional: transmutation reactions)
├── sources.arrow         (optional: decay photon/electron spectra)
└── fission_yields.arrow  (optional: only for fissioning parents)
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
