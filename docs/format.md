# Arrow file format

Each nuclide or element is stored as a directory of Arrow IPC (Feather v2)
files.  Every `.feather` file is a self-contained table that can be
independently memory-mapped.

## Neutron data: `{name}.arrow/`

```text
Li6.arrow/
├── nuclide.feather
├── reactions.feather
├── products.feather
├── distributions.feather
├── urr.feather            (optional)
└── total_nu.feather       (optional)
```

### nuclide.feather

One row of nuclide-level metadata and energy grids.

| Column | Type | Description |
|---|---|---|
| `name` | `utf8` | Nuclide name in GNDS format (e.g. `"Li6"`) |
| `Z` | `int32` | Atomic number |
| `A` | `int32` | Mass number |
| `metastable` | `int32` | Metastable state (0 = ground) |
| `atomic_weight_ratio` | `float64` | Atomic weight ratio |
| `temperatures` | `list<utf8>` | Temperature labels (`["294K", "600K", ...]`) |
| `kTs` | `list<float64>` | kT values in eV |
| `energy_temperatures` | `list<utf8>` | Temperature keys for energy grids |
| `energy_values` | `list<list<float64>>` | Energy grids per temperature |

File-level metadata: `filetype=data_neutron`, `version=3.0`

### reactions.feather

One row per reaction.

| Column | Type | Description |
|---|---|---|
| `mt` | `int32` | ENDF MT reaction number |
| `label` | `utf8` | Reaction name |
| `Q_value` | `float64` | Q-value in eV |
| `center_of_mass` | `bool` | Center-of-mass frame flag |
| `redundant` | `bool` | Redundant reaction flag |
| `xs_temperatures` | `list<utf8>` | Temperature keys |
| `xs_values` | `list<list<float64>>` | Cross-section arrays per temperature |
| `xs_threshold_idx` | `list<int32>` | Threshold index per temperature |
| `n_products` | `int32` | Number of products for this reaction |

### products.feather

One row per product across all reactions.

| Column | Type | Description |
|---|---|---|
| `reaction_mt` | `int32` | Parent reaction MT |
| `product_idx` | `int32` | Product index within the reaction |
| `particle` | `utf8` | Particle type (`"neutron"`, `"photon"`, etc.) |
| `emission_mode` | `utf8` | `"prompt"`, `"delayed"`, or `"total"` |
| `decay_rate` | `float64` | Decay rate (0.0 if prompt) |
| `n_distribution` | `int32` | Number of angle-energy distributions |
| `yield_type` | `utf8` | `"Tabulated1D"` or `"Polynomial"` |
| `yield_data` | `list<float64>` | Yield values (flat) |
| `yield_shape` | `list<int32>` | Shape of the yield array |
| `yield_breakpoints` | `list<int32>` | Interpolation breakpoints |
| `yield_interpolation` | `list<int32>` | Interpolation codes |

### distributions.feather

One row per angle-energy distribution.  Most columns are nullable — only the
columns relevant to each distribution type are populated.

**Common columns:**

- `reaction_mt` (`int32`), `product_idx` (`int32`), `dist_idx` (`int32`)
- `type` (`utf8`): `"uncorrelated"`, `"correlated"`, `"kalbach-mann"`, or `"nbody"`
- `applicability_*`: Applicability tabulation (nullable)

**Uncorrelated angle columns** (nullable):
`angle_energies`, `angle_mu_data`, `angle_mu_offsets`,
`angle_mu_interpolation`

**Uncorrelated energy columns** (nullable):
`energy_dist_type`, `energy_dist_energies`, `energy_dist_interpolation`,
`energy_dist_data`, `energy_dist_offsets`, `energy_dist_out_interp`,
`energy_dist_n_discrete`, `energy_param_x`, `energy_param_y`,
`energy_param2_x`, `energy_param2_y`, `energy_restriction_u`,
`energy_threshold`, `energy_mass_ratio`, `energy_primary_flag`,
`energy_atomic_weight_ratio`, `energy_discrete_energy`

**Correlated columns** (nullable):
`corr_energies`, `corr_breakpoints`, `corr_interpolation`,
`corr_eout_data`, `corr_eout_offsets`, `corr_eout_interp`,
`corr_eout_n_discrete`, `corr_mu_data`, `corr_mu_offsets`,
`corr_mu_interp`

**Kalbach-Mann columns** (nullable):
`km_energies`, `km_breakpoints`, `km_interpolation`,
`km_data`, `km_offsets`, `km_interp`, `km_n_discrete`

**N-body columns** (nullable):
`nbody_n`, `nbody_total_mass`, `nbody_atomic_weight_ratio`,
`nbody_q_value`

### urr.feather

One row per temperature.  Only present if unresolved resonance probability
tables exist.

| Column | Type | Description |
|---|---|---|
| `temperature` | `utf8` | e.g. `"294K"` |
| `energy` | `list<float64>` | URR energy grid |
| `table_data` | `list<float64>` | Flattened probability table (C order) |
| `table_shape` | `list<int32>` | `[n_energy, 6, n_band]` |
| `interpolation` | `int32` | 2 (lin-lin) or 5 (log-log) |
| `inelastic` | `int32` | Inelastic scattering flag |
| `absorption` | `int32` | Absorption flag |
| `multiply_smooth` | `bool` | Multiply-by-smooth-XS flag |

### total_nu.feather

One row.  Only present for fissile nuclides with total nu data.

| Column | Type | Description |
|---|---|---|
| `particle` | `utf8` | `"neutron"` |
| `emission_mode` | `utf8` | `"prompt"`, `"delayed"`, or `"total"` |
| `decay_rate` | `float64` | Decay rate |
| `yield_type` | `utf8` | `"Tabulated1D"` or `"Polynomial"` |
| `yield_data` | `list<float64>` | Yield values |
| `yield_shape` | `list<int32>` | Shape |
| `yield_breakpoints` | `list<int32>` | Interpolation breakpoints |
| `yield_interpolation` | `list<int32>` | Interpolation codes |


## Photon data: `{element}.photon.arrow/`

```text
Fe.photon.arrow/
├── element.feather
├── subshells.feather
├── compton.feather
└── bremsstrahlung.feather
```

### element.feather

One row of element-level data: metadata, union energy grid, main cross
sections, and form factors.

| Column | Type | Description |
|---|---|---|
| `name` | `utf8` | Element symbol |
| `Z` | `int32` | Atomic number |
| `energy` | `list<float64>` | Union energy grid (eV) |
| `coherent_xs` | `list<float64>` | Coherent scattering cross section |
| `incoherent_xs` | `list<float64>` | Incoherent scattering cross section |
| `photoelectric_xs` | `list<float64>` | Photoelectric cross section |
| `pair_production_nuclear_xs` | `list<float64>` | Nuclear pair production cross section |
| `pair_production_electron_xs` | `list<float64>` | Electron pair production cross section |
| `heating_xs` | `list<float64>` | Heating cross section |
| `coherent_int_ff_x` / `_y` | `list<float64>` | Integrated coherent form factor |
| `coherent_ff_x` / `_y` | `list<float64>` | Coherent form factor |
| `coherent_anomalous_real_x` / `_y` | `list<float64>` | Real anomalous scattering factor |
| `coherent_anomalous_imag_x` / `_y` | `list<float64>` | Imaginary anomalous scattering factor |
| `incoherent_ff_x` / `_y` | `list<float64>` | Incoherent scattering function |

File-level metadata: `filetype=data_photon`, `version=3.0`

### subshells.feather

One row per subshell.

| Column | Type | Description |
|---|---|---|
| `designator` | `utf8` | Subshell name (`"K"`, `"L1"`, etc.) |
| `binding_energy` | `float64` | Binding energy (eV) |
| `num_electrons` | `float64` | Number of electrons |
| `xs` | `list<float64>` | Photoionization cross section |
| `threshold_idx` | `int32` | Index into the union energy grid |
| `transitions_data` | `list<float64>` | Flattened transition matrix (nullable) |
| `transitions_shape` | `list<int32>` | Shape of the transition matrix (nullable) |

### compton.feather

One row.  Present when Compton profile data exists.

| Column | Type | Description |
|---|---|---|
| `num_electrons` | `list<float64>` | Electron occupancy per shell |
| `binding_energy` | `list<float64>` | Binding energies (eV) |
| `pz` | `list<float64>` | Electron momentum grid |
| `J_data` | `list<float64>` | Flattened Compton profiles (C order) |
| `J_shape` | `list<int32>` | `[n_shells, n_momentum]` |

### bremsstrahlung.feather

One row.  Present when bremsstrahlung data exists.

| Column | Type | Description |
|---|---|---|
| `I` | `float64` | Mean excitation energy (eV) |
| `electron_energy` | `list<float64>` | Electron kinetic energies (eV) |
| `photon_energy` | `list<float64>` | Reduced photon energies |
| `num_electrons` | `list<float64>` | Subshell occupancies |
| `ionization_energy` | `list<float64>` | Ionization potentials (eV) |
| `dcs_data` | `list<float64>` | Flattened scaled DCS array (C order) |
| `dcs_shape` | `list<int32>` | `[n_electron_energies, n_photon_energies]` |
