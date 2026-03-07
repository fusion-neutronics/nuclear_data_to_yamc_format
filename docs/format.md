# Arrow file format

Each nuclide or element is stored as a directory of Arrow IPC files.  Every
`.arrow` file is a self-contained table that can be independently
memory-mapped.  The format is simulation-ready: synthesized MTs, FastXSGrid
lookup tables, and log-space cross sections are pre-computed.

## Neutron data: `{name}.arrow/`

```text
Li6.arrow/
├── version.json
├── nuclide.arrow
├── reactions.arrow
├── products.arrow
├── distributions.arrow
├── fast_xs.arrow
├── urr.arrow            (optional)
└── total_nu.arrow       (optional)
```

### version.json

Format metadata written at conversion time.

```json
{
  "format_version": 1,
  "library": "fendl-3.2c",
  "converter_version": "0.1.0",
  "created_utc": "2026-02-08T14:30:00Z"
}
```

The reader checks `format_version` and rejects incompatible versions.

### nuclide.arrow

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

File-level metadata: `filetype=data_neutron`, `version=4.0`

### reactions.arrow

One row per reaction, including synthesized MTs (1, 3, 4, 27, 101).

| Column | Type | Description |
|---|---|---|
| `mt` | `int32` | ENDF MT reaction number |
| `label` | `utf8` | Reaction name |
| `Q_value` | `float64` | Q-value in eV |
| `center_of_mass` | `bool` | Center-of-mass frame flag |
| `redundant` | `bool` | Redundant reaction flag (True for synthesized MTs) |
| `xs_temperatures` | `list<utf8>` | Temperature keys |
| `xs_values` | `list<list<float64>>` | Cross-section arrays per temperature |
| `xs_threshold_idx` | `list<int32>` | Threshold index per temperature |
| `n_products` | `int32` | Number of products for this reaction |

### products.arrow

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

### distributions.arrow

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

### fast_xs.arrow

One row per temperature.  Contains the FastXSGrid lookup table for O(1)
cross-section retrieval during transport.

| Column | Type | Description |
|---|---|---|
| `temperature` | `utf8` | e.g. `"294K"` |
| `log_e_min` | `float64` | log of minimum energy |
| `inv_log_delta` | `float64` | Inverse of log bin width |
| `log_grid_index` | `list<int32>` | Grid index per log bin (8001 entries) |
| `xs` | `list<float64>` | Flattened (n_energy, 4) — [total, absorption, scatter, fission] |
| `xs_shape` | `list<int32>` | `[n_energy, 4]` |
| `energy` | `list<float64>` | Full energy grid |
| `scatter_mt_numbers` | `list<int32>` | MT numbers for scattering channels |
| `scatter_mt_xs` | `list<float64>` | Flattened (n_energy, n_scatter) |
| `scatter_mt_shape` | `list<int32>` | Shape of scatter_mt_xs |
| `fission_mt_numbers` | `list<int32>` | MT numbers for fission channels |
| `fission_mt_xs` | `list<float64>` | Flattened (n_energy, n_fission) |
| `fission_mt_shape` | `list<int32>` | Shape of fission_mt_xs |
| `has_partial_fission` | `bool` | True if partial fission MTs present |
| `xs_ngamma` | `list<float64>` | (n,gamma) capture XS |
| `photon_prod` | `list<float64>` | Photon production XS |

### urr.arrow

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

### total_nu.arrow

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


## Photon data: `{element}.arrow/`

```text
Fe.arrow/
├── version.json
├── element.arrow
├── subshells.arrow
├── compton.arrow
└── bremsstrahlung.arrow
```

### version.json

Same structure as neutron version.json.

### element.arrow

One row of element-level data: metadata, union energy grid, main cross
sections (both linear and log-space), and form factors.

| Column | Type | Description |
|---|---|---|
| `name` | `utf8` | Element symbol |
| `Z` | `int32` | Atomic number |
| `energy` | `list<float64>` | Union energy grid (eV) |
| `ln_energy` | `list<float64>` | log of energy grid |
| `coherent_xs` | `list<float64>` | Coherent scattering cross section |
| `ln_coherent_xs` | `list<float64>` | log of coherent XS |
| `incoherent_xs` | `list<float64>` | Incoherent scattering cross section |
| `ln_incoherent_xs` | `list<float64>` | log of incoherent XS |
| `photoelectric_xs` | `list<float64>` | Photoelectric cross section |
| `ln_photoelectric_xs` | `list<float64>` | log of photoelectric XS |
| `pair_production_nuclear_xs` | `list<float64>` | Nuclear pair production cross section |
| `pair_production_electron_xs` | `list<float64>` | Electron pair production cross section |
| `heating_xs` | `list<float64>` | Heating cross section |
| `coherent_int_ff_x` / `_y` | `list<float64>` | Integrated coherent form factor |
| `coherent_ff_x` / `_y` | `list<float64>` | Coherent form factor |
| `coherent_anomalous_real_x` / `_y` | `list<float64>` | Real anomalous scattering factor |
| `coherent_anomalous_imag_x` / `_y` | `list<float64>` | Imaginary anomalous scattering factor |
| `incoherent_ff_x` / `_y` | `list<float64>` | Incoherent scattering function |

File-level metadata: `filetype=data_photon`, `version=4.0`

### subshells.arrow

One row per subshell.

| Column | Type | Description |
|---|---|---|
| `designator` | `utf8` | Subshell name (`"K"`, `"L1"`, etc.) |
| `binding_energy` | `float64` | Binding energy (eV) |
| `num_electrons` | `float64` | Number of electrons |
| `xs` | `list<float64>` | Photoionization cross section |
| `ln_xs` | `list<float64>` | log of photoionization XS |
| `threshold_idx` | `int32` | Index into the union energy grid |
| `transitions_data` | `list<float64>` | Flattened transition matrix (nullable) |
| `transitions_shape` | `list<int32>` | Shape of the transition matrix (nullable) |

### compton.arrow

One row.  Present when Compton profile data exists.  Includes pre-computed
CDFs from trapezoidal integration.

| Column | Type | Description |
|---|---|---|
| `num_electrons` | `list<float64>` | Electron occupancy per shell |
| `binding_energy` | `list<float64>` | Binding energies (eV) |
| `pz` | `list<float64>` | Electron momentum grid |
| `J_data` | `list<float64>` | Flattened Compton profiles (C order) |
| `J_shape` | `list<int32>` | `[n_shells, n_momentum]` |
| `J_cdf_data` | `list<float64>` | Flattened Compton profile CDFs (C order) |
| `J_cdf_shape` | `list<int32>` | `[n_shells, n_momentum]` |

### bremsstrahlung.arrow

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
