# Plan: Simulation-ready Arrow format for yamc

## Why

yamc currently reads nuclear data exclusively from OpenMC-format HDF5 files.
Each time a simulation starts, every nuclide goes through: HDF5 parsing →
CDF normalisation → hierarchical MT synthesis → FastXSGrid construction.  This
work is repeated identically every run.

A simulation-ready Arrow format moves that work to a one-time precompute step.
At simulation time, yamc mmap's the Arrow files and populates its structs with
near-zero load time.

**What this gives us:**
- Faster startup (seconds → sub-second for large libraries)
- MPI shared-node memory (mmap pages shared across ranks via OS page cache)
- A path to eventually dropping the HDF5 C library dependency
- Standard tooling (Arrow files are inspectable with pandas, polars, DuckDB)

**What this does NOT change:**
- Simulation physics — the in-memory `Nuclide` struct is identical
- Runtime memory footprint — same data lives in RAM either way
- Simulation speed — the hot path is unchanged

---

## Will simulation results be identical to the HDF5 route?

**Yes.  Bitwise identical.  Here is why.**

The precompute binary (`yamc-precompute`) calls yamc's own HDF5 reader, which
runs the full loading path:

1. `read_nuclide_from_hdf5()` — parses groups/datasets, reads CDF arrays
   directly from the file, calls `.normalize()` (divides `p` and `c` by the
   last CDF value — a simple `x / scalar` division)
2. `synthesize_hierarchical_mts()` — sums constituent cross sections on the
   full energy grid to produce MT 1, 3, 4, 27, 101
3. `init_fast_xs()` — builds the log-binned O(1) lookup tables

The precompute binary then serialises the resulting `Nuclide` struct to Arrow
IPC.  Arrow writes IEEE 754 floats as raw bytes — no rounding, no
transformation.  Reading them back produces the exact same bit pattern.

**CDFs are not independently recomputed.**  The CDF arrays originate from the
source HDF5 file (OpenMC writes them).  yamc normalises them with a simple
division.  The one exception is Compton profile CDFs (`photon_hdf5.rs`),
which are computed by trapezoidal integration — but the precompute binary
runs that same code, so the output is identical.

If the precompute were reimplemented in Python, floating-point ordering
differences in the MT synthesis and CDF integration could cause last-bit
divergence.  That is why the plan uses a Rust binary that reuses yamc's own
functions — no logic is reimplemented, no results can diverge.

---

## Dual-format verification mode

During development and before the Arrow format is trusted in production, yamc
should support loading the **same nuclide from both formats and comparing them
field-by-field**.  This catches bugs in the Arrow reader/writer that tests
against known values might miss.

### How it works

Add a `--verify-arrow` flag (or equivalent config option) that, when enabled,
makes the loading path do this for every nuclide:

```
1. Load from .yamc.arrow/ via Arrow reader  →  arrow_nuclide
2. Load from .h5 via HDF5 reader            →  hdf5_nuclide
3. Compare every field, panic on mismatch
4. Use arrow_nuclide for simulation (it passed)
```

This is slow (loads everything twice) but gives absolute confidence.  It runs
during CI and during the initial validation phase, then gets turned off for
production runs.

### What gets compared

```rust
/// Compare two Nuclides field-by-field. Panics on any mismatch.
pub fn assert_nuclides_identical(a: &Nuclide, b: &Nuclide) {
    assert_eq!(a.name, b.name);
    assert_eq!(a.atomic_number, b.atomic_number);
    assert_eq!(a.mass_number, b.mass_number);
    assert_eq!(a.atomic_weight_ratio, b.atomic_weight_ratio);
    assert_eq!(a.loaded_temperatures, b.loaded_temperatures);
    assert_eq!(a.fissionable, b.fissionable);

    // Energy grids
    let a_energy = a.energy.as_ref().unwrap();
    let b_energy = b.energy.as_ref().unwrap();
    for temp in &a.loaded_temperatures {
        assert_eq!(a_energy[temp], b_energy[temp],
            "energy grid mismatch for temp {temp}");
    }

    // Reactions (per temperature, per MT)
    for (temp_idx, temp) in a.loaded_temperatures.iter().enumerate() {
        let a_rxns = &a.reactions[temp_idx];
        let b_rxns = &b.reactions[temp_idx];
        assert_eq!(a_rxns.len(), b_rxns.len(),
            "reaction count mismatch for temp {temp}");
        for (mt, a_r) in a_rxns {
            let b_r = b_rxns.get(mt).unwrap_or_else(||
                panic!("MT {mt} missing in second nuclide at temp {temp}"));
            assert_eq!(a_r.cross_section, b_r.cross_section,
                "XS mismatch for MT {mt} at temp {temp}");
            assert_eq!(a_r.threshold_idx, b_r.threshold_idx);
            assert_eq!(a_r.energy, b_r.energy);
            assert_eq!(a_r.q_value, b_r.q_value);
            assert_eq!(a_r.scatter_in_cm, b_r.scatter_in_cm);
            assert_eq!(a_r.redundant, b_r.redundant);
            assert_eq!(a_r.products.len(), b_r.products.len(),
                "product count mismatch for MT {mt}");
            // ... compare products, distributions, yields ...
        }
    }

    // FastXSGrid
    for i in 0..a.fast_xs.len() {
        assert_eq!(a.fast_xs[i].xs, b.fast_xs[i].xs,
            "fast_xs.xs mismatch at temp_idx {i}");
        assert_eq!(a.fast_xs[i].energy, b.fast_xs[i].energy);
        assert_eq!(a.fast_xs[i].log_grid_index, b.fast_xs[i].log_grid_index);
        assert_eq!(a.fast_xs[i].log_e_min, b.fast_xs[i].log_e_min);
        assert_eq!(a.fast_xs[i].inv_log_delta, b.fast_xs[i].inv_log_delta);
    }

    // URR data
    for i in 0..a.urr_data.len() {
        // ... compare UrrData fields ...
    }
}
```

This function lives in the main crate (not behind a feature gate) so it can
be used in integration tests regardless of build configuration.

### Verification levels

| Level | What it does | When to use |
|---|---|---|
| `--verify-arrow` | Load both formats for every nuclide, compare, panic on mismatch | Development, CI, first use of a new Arrow library |
| Integration test | Compare a handful of known nuclides (Li6, U235, U238, Fe56) | Every CI run |
| Simulation-level test | Run same problem with HDF5 and Arrow data, compare tallies | Release validation |
| None (default) | Load Arrow only, trust the data | Production runs |

---

## Pipeline overview

```
Current (HDF5 only):
  ENDF → ACE → OpenMC → .h5 → yamc load (parse + normalise + synthesize + fast_xs) → simulate

New (simulation-ready Arrow):
  ENDF → ACE → OpenMC → .h5 → yamc-precompute (Rust binary) → .yamc.arrow/
                                    │
                                    └── runs yamc's own HDF5 reader +
                                        synthesize_hierarchical_mts() +
                                        init_fast_xs(), then serialises
                                        the fully-populated Nuclide to Arrow

  At simulation time:
  .yamc.arrow/ → yamc mmap + populate structs → simulate
```

The `.yamc.arrow/` directory contains the fully-loaded, post-processed state
of the `Nuclide` struct — exactly what sits in memory after the HDF5 reader
finishes.

We use `.yamc.arrow` (not plain `.arrow`) to distinguish the simulation-ready
format from the raw data format that `nuclear_data_to_yamc_format` produces.

---

## Format versioning

The `.yamc.arrow/` format is produced by `yamc-precompute`, which runs yamc's
own loading code to populate the data.  If that code changes (new fields, bug
fixes to synthesis logic, different CDF handling), old Arrow files become stale.
Silent use of stale data would give wrong simulation results with no error.

Every `.yamc.arrow/` directory contains a `version.json` file written by
`yamc-precompute` at generation time.  yamc checks this file at load time and
refuses to proceed if the versions are incompatible.

### `version.json` contents

```json
{
  "format_version": 1,
  "library": "fendl-3.2c",
  "yamc_version": "0.3.0",
  "yamc_precompute_version": "0.1.0",
  "created_utc": "2026-02-08T14:30:00Z",
  "source_h5": "Li6.h5",
  "source_h5_sha256": "a1b2c3d4..."
}
```

| Field | Purpose |
|---|---|
| `format_version` | Integer schema version.  Bumped when the Arrow column layout changes (new columns, renamed columns, different encoding).  yamc checks `format_version` and rejects files with an unsupported version. |
| `library` | Nuclear data library name (e.g. `"fendl-3.2c"`, `"endfb-8.0"`, `"jeff-3.3"`).  Passed to `yamc-precompute` via `--library` flag.  Informational — lets users and yamc know which library the data came from.  yamc can optionally warn if a simulation mixes nuclides from different libraries. |
| `yamc_version` | The version of the `yamc` crate that generated the file (from `Cargo.toml`).  Informational — lets the user know which yamc produced the data. |
| `yamc_precompute_version` | The version of the `yamc-precompute` binary.  Informational. |
| `created_utc` | Timestamp.  Informational — helps users know how old the data is. |
| `source_h5` | Original HDF5 filename.  Informational — provenance tracking. |
| `source_h5_sha256` | SHA-256 of the source `.h5` file.  Allows detecting if the source data has changed since the Arrow files were generated. |

### Version checking at load time

```rust
// In nuclide_arrow.rs:
const SUPPORTED_FORMAT_VERSIONS: &[u32] = &[1];

fn check_version(dir: &Path) -> Result<(), Box<dyn Error>> {
    let version_path = dir.join("version.json");
    let version: VersionInfo = serde_json::from_reader(
        std::fs::File::open(&version_path)?
    )?;

    if !SUPPORTED_FORMAT_VERSIONS.contains(&version.format_version) {
        return Err(format!(
            "Unsupported .yamc.arrow format version {} (supported: {:?}). \
             Regenerate with: yamc-precompute {}",
            version.format_version,
            SUPPORTED_FORMAT_VERSIONS,
            version.source_h5,
        ).into());
    }
    Ok(())
}
```

The error message tells the user exactly what to do: rerun `yamc-precompute` on
the original `.h5` file.

### When to bump `format_version`

| Change | Bump? |
|---|---|
| Add a new optional column to a Arrow file | No — old files still work, the column is absent |
| Add a new required column | **Yes** — old files missing the column would fail |
| Rename or remove a column | **Yes** |
| Change the encoding of a column (e.g. AoS → SoA) | **Yes** |
| Bug fix in `synthesize_hierarchical_mts` or `init_fast_xs` | **Yes** — old data is wrong |
| New yamc version but no schema change | No — `yamc_version` field is informational |

This means `format_version` changes rarely.  Most yamc releases won't require
regenerating Arrow files.

---

## Simulation-ready Arrow directory layout

### Neutron: `Li6.yamc.arrow/`

```
Li6.yamc.arrow/
├── version.json             – format version, yamc version, provenance
├── nuclide.arrow            – 1 row: name, Z, A, AWR, temperatures, energy grids
├── reactions.arrow           – N rows per temperature: mt, Q_value, xs, threshold_idx,
│                               scatter_in_cm, redundant (includes synthesised MT 1/3/4/27/101)
├── products.arrow            – N rows: particle, emission_mode, yield
├── distributions.arrow       – N rows: type + all angle/energy data (CDFs already normalised)
├── fast_xs.arrow             – 1 row per temperature: log_grid_index, xs[total,abs,scat,fiss],
│                               energy, scatter_mt data, fission_mt data
├── urr.arrow                 – N rows per temperature (optional)
└── total_nu.arrow            – 1 row (optional, fissile only)
```

### Photon: `Fe.photon.yamc.arrow/`

```
Fe.photon.yamc.arrow/
├── version.json             – format version, yamc version, provenance
├── element.arrow             – 1 row: name, Z, ln(energy), ln(xs), form factors
├── subshells.arrow           – N rows: designator, binding_energy, ln(xs), transitions
├── compton.arrow             – 1 row: profiles with CDF (optional)
└── bremsstrahlung.arrow      – 1 row: DCS data (optional)
```

### What the simulation-ready format stores vs what HDF5 stores

| Data | HDF5 file | `.yamc.arrow/` | Notes |
|---|---|---|---|
| Reaction cross sections | Per-reaction, threshold-truncated | Same — per-reaction, threshold-truncated | Identical data |
| Hierarchical MTs (1,3,4,27,101) | Absent | Present as rows in reactions.arrow | Avoids re-synthesis at load |
| Per-reaction energy grids | Implicit (threshold_idx into top-level) | Explicit Vec per reaction | Avoids re-slicing at load |
| Angle/energy CDFs | Present (from OpenMC) | Present (post yamc `.normalize()`) | Only difference: divided by integral |
| FastXSGrid | Absent | Present in fast_xs.arrow | Avoids `init_fast_xs()` at load |
| Photon XS / energy | Raw eV / barns | ln(E) / ln(sigma) | Avoids log-transform at load |
| Compton profile CDFs | Absent | Present (computed by yamc's own code) | Avoids trapezoidal integration at load |

---

## Changes

### 1. Add `arrow` feature + deps to `Cargo.toml`

**File:** `yamc/Cargo.toml`

```toml
[dependencies]
arrow-ipc    = { version = "54", optional = true }
arrow-schema = { version = "54", optional = true }
arrow-array  = { version = "54", optional = true }
arrow-cast   = { version = "54", optional = true }

[features]
arrow = ["dep:arrow-ipc", "dep:arrow-schema", "dep:arrow-array", "dep:arrow-cast"]
```

### 2. Refactor shared post-processing out of `nuclide_hdf5.rs`

**File:** `yamc/src/nuclide_hdf5.rs` → extract into `yamc/src/nuclide.rs`

`synthesize_hierarchical_mts()` is currently private to `nuclide_hdf5.rs`
(line 573).  `init_fast_xs()` is already a public method on `Nuclide`
(line 820).

Make `synthesize_hierarchical_mts()` a public method on `Nuclide` so the
precompute binary can call it:

```rust
// Move from nuclide_hdf5.rs to nuclide.rs:
impl Nuclide {
    pub fn synthesize_hierarchical_mts(&mut self) { ... }
    // init_fast_xs() is already pub
}
```

The HDF5 reader calls these the same way — no behaviour change.

### 3. Create `nuclide_arrow.rs` — Arrow reader and writer

**File:** `yamc/src/nuclide_arrow.rs` (new, ~800-1000 lines)

Contains both reading and writing functions:

**Writer** (used by `yamc-precompute`):
```rust
pub fn write_nuclide_to_arrow(nuclide: &Nuclide, dir: &Path) -> Result<()>;
```

**Reader** (used at simulation time):
```rust
pub fn read_nuclide_from_arrow(
    path: &Path,
    temps_filter: Option<&HashSet<String>>,
) -> Result<Nuclide>;
```

The reader:
1. mmap's each `.arrow` file in the `.yamc.arrow/` directory
2. Reads Arrow RecordBatch columns into yamc's Rust structs
3. Re-links `Arc<Reaction>` pointers in FastXSGrid by MT number lookup
4. Returns a fully populated `Nuclide` — no synthesis or fast_xs building

#### Column mappings

**nuclide.arrow → Nuclide header:**

| Arrow column | yamc field | Notes |
|---|---|---|
| `name` | `Nuclide.name` | String direct |
| `Z` | `Nuclide.atomic_number` | i32 → u32 |
| `A` | `Nuclide.mass_number` | i32 → u32 |
| `atomic_weight_ratio` | `Nuclide.atomic_weight_ratio` | f64 direct |
| `temperatures` | `available_temperatures` / `loaded_temperatures` | `list<utf8>`, apply temps_filter |
| `energy_temperatures` + `energy_values` | `Nuclide.energy` | `HashMap<String, Vec<f64>>` |

Derived fields: `atomic_symbol` (strip digits from name), `neutron_number`
(A - Z), `fissionable` (scan reactions for fission MTs).

**reactions.arrow → Reaction structs:**

| Arrow column | yamc field | Notes |
|---|---|---|
| `temperature` | — | Group by this to build `Vec<HashMap<i32, Reaction>>` |
| `mt` | `Reaction.mt_number` | i32 direct |
| `Q_value` | `Reaction.q_value` | f64 direct |
| `cross_section` | `Reaction.cross_section` | `list<f64>` — already on full grid |
| `energy` | `Reaction.energy` | `list<f64>` — already threshold-sliced |
| `threshold_idx` | `Reaction.threshold_idx` | i32 → usize |
| `center_of_mass` | `Reaction.scatter_in_cm` | bool |
| `redundant` | `Reaction.redundant` | bool |

Includes pre-synthesised MT 1, 3, 4, 27, 101 rows — no synthesis needed.

**products.arrow → ReactionProduct structs:**

| Arrow column | yamc field | Notes |
|---|---|---|
| `reaction_mt` + `product_idx` | — | Link to parent reaction |
| `particle` | `ReactionProduct.particle` | String → `ParticleType` enum |
| `emission_mode` | `ReactionProduct.emission_mode` | String direct |
| `decay_rate` | `ReactionProduct.decay_rate` | f64 direct |
| `yield_type` + `yield_data` + ... | `ReactionProduct.product_yield` | Reconstruct `Yield` variant |

**distributions.arrow → AngleEnergyDistribution variants:**

CDFs have already been through yamc's `.normalize()` pass.  They originate
from the source data (OpenMC/HDF5) and are only divided by the last CDF value
— no independent recomputation.  The Arrow reader does not need to call
`.normalize()`.

Distribution type dispatch by `type` column:
- `"uncorrelated"` → `UncorrelatedAngleEnergy`
- `"correlated"` → `CorrelatedAngleEnergy`
- `"kalbach-mann"` → `KalbachMann`
- `"nbody"` → `NBodyPhaseSpace`
- `"evaporation"` → `Evaporation`

**fast_xs.arrow → FastXSGrid structs:**

| Arrow column | yamc field | Notes |
|---|---|---|
| `temperature` | — | Index alignment with `loaded_temperatures` |
| `log_e_min` | `FastXSGrid.log_e_min` | f64 direct |
| `inv_log_delta` | `FastXSGrid.inv_log_delta` | f64 direct |
| `log_grid_index` | `FastXSGrid.log_grid_index` | `list<u32>` → `Vec<usize>` (8192 entries) |
| `xs` | `FastXSGrid.xs` | `FixedSizeList<f64, 4>` → `Vec<[f64; 4]>` |
| `energy` | `FastXSGrid.energy` | `list<f64>` direct |
| `scatter_mt_numbers` | `FastXSGrid.scatter_mt_xs[].0` | `list<i32>` |
| `scatter_mt_xs` | `FastXSGrid.scatter_mt_xs[].1` | `list<list<f64>>` |
| `fission_mt_numbers` | `FastXSGrid.fission_mt_xs[].0` | `list<i32>` |
| `fission_mt_xs` | `FastXSGrid.fission_mt_xs[].1` | `list<list<f64>>` |
| `xs_ngamma` | `FastXSGrid.xs_ngamma` | `list<f64>` direct |

The `Arc<Reaction>` pointers in `scatter_mt_xs` and `fission_mt_xs` cannot be
serialised.  They are re-linked at load time by looking up the MT number in the
already-loaded reactions HashMap.  This is a cheap O(n) pass, not a
recomputation.

**urr.arrow → UrrData:** Flat `table_data` + `table_shape` reshaped into
`cdf_values` + `xs_values`.  Same as raw format.

**total_nu.arrow → FissionNuData:** `yield_data` with shape `[2, N]` →
energy + nu vectors.  Same as raw format.

### 4. Create `photon_arrow.rs` — photon Arrow reader and writer

**File:** `yamc/src/photon_arrow.rs` (new, ~400-500 lines)

Reads `.photon.yamc.arrow/` directories.  Cross-sections and energy grids are
already in log space.  Compton profile CDFs are already computed (by yamc's
own `photon_hdf5.rs` code during precompute).

| Arrow column | yamc field | Notes |
|---|---|---|
| `energy` | `PhotonInteraction.energy` | Already ln(E) — direct |
| `coherent_xs` etc. | `PhotonInteraction.coherent_xs` | Already ln(sigma) — direct |
| `coherent_int_ff_x/y` | `coherent_int_form_factor` | Build `Tabulated1D` |
| `J_data` + `J_cdf` | `profile_pdf` + `profile_cdf` | Already computed — direct |
| Subshell `xs` | `ElectronSubshell.cross_section` | Already ln(sigma) — direct |
| Subshell transitions | `ElectronSubshell.transitions` | Already CDF — direct |

### 5. Update `nuclide_loader.rs` — register Arrow format

**File:** `yamc/src/nuclide_loader.rs`

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NuclideFormat {
    #[default]
    Hdf5,
    #[cfg(feature = "arrow")]
    Arrow,
}

impl NuclideFormat {
    pub fn loader(&self) -> Box<dyn NuclideLoader> {
        match self {
            NuclideFormat::Hdf5 => Box::new(Hdf5Loader),
            #[cfg(feature = "arrow")]
            NuclideFormat::Arrow => Box::new(ArrowLoader),
        }
    }

    pub fn from_extension(path: &Path) -> Option<Self> {
        // Check directory name (Li6.yamc.arrow/ or Li6.yamc.arrow)
        let name = path.file_name()?.to_str()?;
        if name.ends_with(".yamc.arrow") {
            #[cfg(feature = "arrow")]
            return Some(NuclideFormat::Arrow);
            #[cfg(not(feature = "arrow"))]
            return None;
        }
        let ext = path.extension()?.to_str()?.to_lowercase();
        match ext.as_str() {
            "h5" | "hdf5" => Some(NuclideFormat::Hdf5),
            _ => None,
        }
    }
}
```

Both formats coexist.  The user controls which is used by the file path they
pass to `material.read_nuclide()`.  Passing `Li6.h5` uses HDF5; passing
`Li6.yamc.arrow` uses Arrow.  No global switch needed.

### 6. Register modules in `lib.rs`

**File:** `yamc/src/lib.rs`

```rust
#[cfg(feature = "arrow")]
pub mod nuclide_arrow;
#[cfg(feature = "arrow")]
pub mod photon_arrow;
```

### 7. Create `yamc-precompute` binary

**File:** `yamc/src/bin/yamc_precompute.rs` (new, ~200-300 lines)

```
Usage:
  yamc-precompute neutron Li6.h5 -o Li6.yamc.arrow/
  yamc-precompute photon  Fe.h5  -o Fe.photon.yamc.arrow/
  yamc-precompute all     /path/to/h5/library/ -o /path/to/arrow/library/
```

The binary:

1. Takes `.h5` file path(s) as input
2. Loads via `read_nuclide_from_hdf5()` — this runs the full loading path
   including CDF normalisation, hierarchical MT synthesis, and `init_fast_xs()`.
   **This is the same code that runs when yamc loads an HDF5 file for
   simulation, which is what guarantees bitwise identical results.**
3. Serialises the fully-populated `Nuclide` to `.yamc.arrow/` using
   `write_nuclide_to_arrow()`
4. Does the same for photon data

This binary depends on both `hdf5` and `arrow` features.  It is a build-time
tool, not needed at simulation time.

### 8. Optional: Add photon loader abstraction

**File:** `yamc/src/photon_loader.rs` (new, small)

Currently `photon_hdf5.rs` is called directly with no trait dispatch.
For symmetry with `NuclideLoader`:

```rust
pub trait PhotonLoader: Send + Sync {
    fn load(&self, path: &Path) -> Result<PhotonInteraction, Box<dyn Error>>;
}
```

---

## Files changed (summary)

| File | Action | Est. lines |
|---|---|---|
| `yamc/Cargo.toml` | Edit: add arrow deps + feature + binary | ~10 |
| `yamc/src/lib.rs` | Edit: add module declarations | ~4 |
| `yamc/src/nuclide_hdf5.rs` | Refactor: move `synthesize_hierarchical_mts` to nuclide.rs | ~5 |
| `yamc/src/nuclide.rs` | Edit: receive extracted function + `assert_nuclides_identical` | ~250 |
| `yamc/src/nuclide_loader.rs` | Edit: add Arrow variant + ArrowLoader | ~50 |
| `yamc/src/nuclide_arrow.rs` | **New**: Arrow reader + writer | ~800-1000 |
| `yamc/src/photon_arrow.rs` | **New**: photon Arrow reader + writer | ~400-500 |
| `yamc/src/bin/yamc_precompute.rs` | **New**: precompute CLI binary | ~200-300 |
| `yamc/src/photon_loader.rs` | **New** (optional): photon loader trait | ~60 |

---

## Implementation order

1. **Refactor `synthesize_hierarchical_mts`** — move to `nuclide.rs` as a pub
   method.  `init_fast_xs` is already pub.  Run existing tests.

2. **Add `assert_nuclides_identical`** — the comparison function.  Test it by
   loading the same HDF5 file twice and comparing (must pass trivially).

3. **Add `arrow` feature + deps** to `Cargo.toml`, module stubs in `lib.rs`.

4. **Implement Arrow writer** (`write_nuclide_to_arrow`,
   `write_photon_to_arrow`) — defines the on-disk schema.

5. **Implement `yamc-precompute` binary** — load HDF5, write Arrow.  Then
   immediately test: load the Arrow back, compare against HDF5 using
   `assert_nuclides_identical`.

6. **Implement Arrow reader** (`read_nuclide_from_arrow`,
   `read_photon_from_arrow`).

7. **Wire up `nuclide_loader.rs`** — add `Arrow` variant and `ArrowLoader`.

8. **Dual-format verification tests** — load same nuclide from both formats,
   compare everything.

9. **Simulation-level comparison** — run identical simulation with HDF5 and
   Arrow data, same RNG seed, compare tallies.

10. **Remove HDF5 dependency** — once verification passes (see below).

---

## Phase 2: Removing HDF5

Once the dual-format verification tests pass — bitwise field-by-field for a
representative set of nuclides AND simulation-level with identical particle
histories — the HDF5 reading code is no longer needed at simulation time.

### What gets removed

| File | Action |
|---|---|
| `yamc/src/nuclide_hdf5.rs` (~1500 lines) | Delete |
| `yamc/src/photon_hdf5.rs` (~800 lines) | Delete |
| `Cargo.toml`: `hdf5 = { package = "hdf5-metno", ... }` | Remove from `[dependencies]` |
| `nuclide_loader.rs`: `Hdf5Loader`, `NuclideFormat::Hdf5` | Remove |
| All `use hdf5::` imports | Remove |
| `Cargo.toml`: `exclude = ["tests/*.h5", ...]` | Remove |

### What this gains

- **~2300 lines of reader code deleted** — less to maintain
- **Drop the `hdf5-metno` crate + the HDF5 C library** — this is a heavy
  dependency.  It requires a C compiler, CMake, and static linking (or a system
  libhdf5).  It is the single most painful dependency to build, especially on
  CI, Windows, and cross-compilation targets.
- **Simpler build** — `cargo build` no longer needs to compile the HDF5 C
  library from source (which the `static` feature currently does)
- **Smaller binary** — the HDF5 C library adds significant size
- **Reduced attack surface** — the HDF5 C library has had multiple CVEs; the
  pure-Rust Arrow reader has no C code

### What stays

- **`yamc-precompute` binary** — this still needs HDF5 to read the source
  `.h5` files.  It becomes the only place HDF5 is used, and it runs once
  per library, not per simulation.  It can be a separate crate or a
  feature-gated binary (`cargo build --features hdf5,arrow --bin
  yamc_precompute`).
- **`synthesize_hierarchical_mts()` and `init_fast_xs()`** — these live in
  `nuclide.rs`, not in the HDF5 reader, so they survive the deletion.  They
  are still needed by the precompute binary.
- **All simulation code** — nothing in the physics, geometry, tallying, or
  transport code touches the file format.  The `Nuclide` struct is the
  boundary.

### When to do it

Not immediately.  The removal happens after:

1. All dual-format verification tests pass for a full nuclear data library
2. At least one real simulation has been validated (k-effective, tallies)
3. The Arrow format has been used in production for long enough to catch
   any edge cases (unusual nuclides, metastable states, URR edge cases)

A reasonable timeline: keep both formats for one release cycle, then remove
HDF5 in the following release.

### Migration path for users

```
Before:  user has .h5 library
         ↓
         yamc-precompute all /path/to/h5/ -o /path/to/arrow/
         ↓
After:   user points yamc at .yamc.arrow/ library
         (yamc no longer needs HDF5 at all)
```

Users who already have `.h5` libraries run `yamc-precompute` once.  New users
never touch HDF5 — the conversion pipeline (`nuclear_data_to_yamc_format`)
could be extended to produce `.yamc.arrow/` directly in the future, bypassing
HDF5 entirely.

---

## Testing strategy

### Level 1: Bitwise field-by-field comparison

```rust
#[test]
fn arrow_matches_hdf5_li6() {
    let h5 = read_nuclide_from_hdf5("Li6.h5", None).unwrap();
    let ar = read_nuclide_from_arrow("Li6.yamc.arrow/", None).unwrap();
    assert_nuclides_identical(&h5, &ar);
}
```

Run for a set of representative nuclides: Li6 (light, simple), Fe56 (medium),
U235 (fissile, many products), U238 (heavy, URR, many reactions), H1 (minimal).

### Level 2: Simulation-level comparison

Run the same simulation with both data sources, same RNG seed:

```rust
#[test]
fn simulation_identical_hdf5_vs_arrow() {
    let result_h5 = run_simulation(data_path="library_h5/", seed=42);
    let result_ar = run_simulation(data_path="library_arrow/", seed=42);
    assert_eq!(result_h5.k_effective, result_ar.k_effective);
    assert_eq!(result_h5.tallies, result_ar.tallies);
}
```

Because the `Nuclide` structs are bitwise identical, the same RNG seed must
produce identical particle histories — not just statistically consistent, but
exactly the same numbers.

### Level 3: Feature gate test

Compile without the `arrow` feature.  All existing HDF5-only tests must pass
unchanged.

### Level 4: Performance test

Measure load time for a full library (~400 nuclides):
- HDF5: parse + synthesise + build fast_xs
- Arrow: mmap + populate structs

---

## Size and memory considerations

### Disk size

The simulation-ready format is larger than HDF5 because it includes
`FastXSGrid` data that the HDF5 path computes at load time.

`FastXSGrid` contains per-reaction cross-sections re-interpolated onto the
full nuclide energy grid.  Every scattering/fission reaction gets a `Vec<f64>`
of the same length as the top-level grid — even threshold reactions where most
values are zeros.

Rough estimate for U238 (~150k energy points, ~30 reactions, 1 temperature):

| Component | Size |
|---|---|
| `fast_xs.xs` (interleaved `[total, abs, scatter, fission]`) | ~4.8 MB |
| `fast_xs.scatter_mt_xs` (~20 reactions * full grid) | ~24 MB |
| `fast_xs.fission_mt_xs` (~5 reactions * full grid) | ~6 MB |
| Reactions, products, distributions, nuclide header | ~1-2 MB |
| **Total per nuclide** | **~35 MB** |

Full library (~400 nuclides, 1 temperature): ~14 GB uncompressed vs ~2-3 GB
for HDF5.

**Recommended: LZ4 compression.**  Arrow IPC supports per-column compression.
The long runs of zeros in threshold reactions compress at ~10:1.  LZ4
decompresses at ~4 GB/s — negligible load-time cost.  Compressed library size:
~3-5 GB, comparable to HDF5.

### Runtime memory

**The runtime memory footprint is identical regardless of format.**  The
`FastXSGrid` data lives in heap memory for the entire simulation whether it
was parsed from HDF5 or mmap'd from Arrow.  The Arrow format does not increase
memory usage — it just makes the disk cost of that data visible.

A typical simulation loads 20-30 nuclides (materials in the problem), not the
full 400-nuclide library.  At ~35 MB each that's ~0.7-1 GB for nuclide data,
leaving plenty of room for mesh tallies.

Where Arrow mmap helps with memory: **MPI multi-rank on shared nodes**.  100
ranks loading from HDF5 each hold their own copy in heap memory (100 * 1 GB =
100 GB).  100 ranks mmap'ing the same Arrow files share physical pages via the
OS page cache (~1 GB shared + per-rank overhead).

---

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Schema drift: yamc changes struct layout, stale `.yamc.arrow` files silently give wrong results | `version.json` in every `.yamc.arrow/` directory contains `format_version` (integer).  Reader checks it at load time, rejects mismatches with a clear error message telling the user to rerun `yamc-precompute`.  See "Format versioning" section above. |
| `arrow-rs` crate slows builds | Feature-gated — only users who opt in pay the cost.  The `hdf5-metno` crate with static linking is arguably heavier. |
| `Arc<Reaction>` pointers can't be serialised | Store MT numbers + XS data.  Re-link `Arc` pointers at load time via reaction HashMap lookup (cheap O(n) pass). |
| Directory path handling (trailing slash or not) | `from_extension` checks the path name suffix, not the filesystem extension.  Works for both `Li6.yamc.arrow/` and `Li6.yamc.arrow`. |
| Precompute binary requires HDF5 at build time | It's a build-time tool in a separate binary.  After Phase 2, HDF5 is only needed for `yamc-precompute`, not for `yamc` itself. |
| User doesn't trust Arrow results | Dual-format verification mode loads both and compares.  Run it once per library to build confidence.  Phase 2 removal only happens after verification passes. |
| Removing HDF5 is irreversible | It's not — `nuclide_hdf5.rs` lives in git history and can be restored.  The `NuclideLoader` trait means re-adding HDF5 support later is a clean addition, not a refactor.  But the real safety net is `yamc-precompute`: as long as it can read `.h5` files, users are never stranded. |

---

## PR breakdown

### Safety strategy: why this can't leave us stuck

Every PR below is either a **pure refactor** (no new functionality, existing tests
must pass) or an **additive feature behind `#[cfg(feature = "arrow")]`**.  At no
point does an unmerged or failed PR break the existing HDF5 path.

If the entire Arrow effort were abandoned after any PR, the worst case is:
- PRs 1-2 remain as small improvements to the main crate (pub method, comparison utility)
- PRs 3-7 are behind a feature gate — `cargo build` without `--features arrow` compiles exactly the same binary as before

The `arrow` feature is **opt-in**.  Default builds never touch Arrow code.  This
means every PR can be merged to `main` and released without affecting any user
who hasn't explicitly opted in.

### Go/no-go checkpoints

Between certain PRs, verification must pass before proceeding.  These are marked
with a checkpoint icon below.  If a checkpoint fails, we stop and debug — the
codebase is in a working state at every point.

---

### PR 1 — Refactor: make `synthesize_hierarchical_mts` public
**Repo:** `yamc`
**Risk:** None — pure refactor

Move `synthesize_hierarchical_mts()` from `nuclide_hdf5.rs` (where it's private)
to `nuclide.rs` as `impl Nuclide { pub fn synthesize_hierarchical_mts(...) }`.
The HDF5 reader calls the new location.  `init_fast_xs()` is already pub.

**What to check:** All existing tests pass.  `cargo test` green.  No behaviour change.

**Files changed:**
- `src/nuclide.rs` — add the method
- `src/nuclide_hdf5.rs` — replace inline code with `nuclide.synthesize_hierarchical_mts()`

---

### PR 2 — Add `assert_nuclides_identical()` comparison function
**Repo:** `yamc`
**Risk:** None — adds a test utility, no production code changes

Add `pub fn assert_nuclides_identical(a: &Nuclide, b: &Nuclide)` to `nuclide.rs`.
Compares every field: metadata, energy grids, reactions (per temperature, per MT),
products, distributions, FastXSGrid, URR data, fission nu.

Add a trivial test that loads the same HDF5 file twice and compares — must pass.

**What to check:** `cargo test` green.  The identity test passes.

**Files changed:**
- `src/nuclide.rs` — add comparison function + test

---

### PR 3 — Add `arrow` feature gate + Arrow deps
**Repo:** `yamc`
**Risk:** None — adds optional dependencies, no code uses them yet

Add `arrow-ipc`, `arrow-schema`, `arrow-array`, `arrow-cast` as optional deps.
Add `arrow` feature flag.  Add empty `nuclide_arrow.rs` and `photon_arrow.rs`
module stubs behind `#[cfg(feature = "arrow")]`.

**What to check:** `cargo build` (no features) still works.
`cargo build --features arrow` compiles.

**Files changed:**
- `Cargo.toml` — add deps + feature
- `src/lib.rs` — add conditional module declarations
- `src/nuclide_arrow.rs` — empty stub
- `src/photon_arrow.rs` — empty stub

---

### PR 4 — Implement Arrow writer + `yamc-precompute` binary
**Repo:** `yamc`
**Risk:** Low — new code, doesn't touch existing paths

Implement `write_nuclide_to_arrow()` and `write_photon_to_arrow()` in
`nuclide_arrow.rs` / `photon_arrow.rs`.  Create the `yamc-precompute` binary
in `src/bin/yamc_precompute.rs`.

The binary loads HDF5 via yamc's existing reader, then serialises to
`.yamc.arrow/` directories.

**What to check:** Run `yamc-precompute` on a few nuclides (H1, Li6, Fe56, U235,
U238).  Inspect the output directories — they should contain the expected
`.arrow` files with reasonable sizes.

**Files changed:**
- `src/nuclide_arrow.rs` — writer implementation
- `src/photon_arrow.rs` — writer implementation
- `src/bin/yamc_precompute.rs` — new binary
- `Cargo.toml` — add `[[bin]]` entry

---

### CHECKPOINT: round-trip test before writing the reader

At this point we have a writer but no reader.  Before investing in the reader,
add a **write-then-read-back test** using Arrow's own API to verify the Arrow
files contain the right data.  This catches schema bugs early.

---

### PR 5 — Implement Arrow reader
**Repo:** `yamc`
**Risk:** Low — new code behind feature gate

Implement `read_nuclide_from_arrow()` and `read_photon_from_arrow()`.  The
reader mmap's Arrow files, reconstructs `Nuclide`/`PhotonInteraction` structs,
and re-links `Arc<Reaction>` pointers by MT number lookup.

**What to check:** Unit tests that write a `Nuclide` to Arrow and read it back,
then compare using `assert_nuclides_identical()`.

**Files changed:**
- `src/nuclide_arrow.rs` — reader implementation
- `src/photon_arrow.rs` — reader implementation

---

### CHECKPOINT: dual-format field-by-field verification

This is the critical gate.  For each representative nuclide:
```
HDF5 → Nuclide (h5)
HDF5 → yamc-precompute → .yamc.arrow/ → Nuclide (arrow)
assert_nuclides_identical(&h5, &arrow)  // must pass
```

If this fails, the Arrow writer or reader has a bug.  Fix it before proceeding.
**Do not merge PR 6 until this checkpoint passes for H1, Li6, Fe56, U235, U238.**

---

### PR 6 — Wire Arrow reader into `nuclide_loader.rs`
**Repo:** `yamc`
**Risk:** Low — extends existing dispatch, HDF5 path unchanged

Add `NuclideFormat::Arrow` variant.  Add `ArrowLoader` implementing
`NuclideLoader`.  Update `from_extension()` to recognise `.yamc.arrow`.

Users can now pass `.yamc.arrow` paths to `material.read_nuclide()`.

**What to check:** Integration test that loads a nuclide via the `NuclideLoader`
trait with both formats and compares.  All existing HDF5 tests still pass.

**Files changed:**
- `src/nuclide_loader.rs` — add Arrow variant + loader

---

### PR 7 — Dual-format verification mode + simulation-level tests
**Repo:** `yamc`
**Risk:** None — test infrastructure only

Add `--verify-arrow` flag support.  Add simulation-level comparison test: run the
same problem with HDF5 and Arrow data, same RNG seed, compare k-effective and
tallies.

**What to check:** The simulation test produces bitwise identical results.

**Files changed:**
- Integration test files
- Config/CLI support for `--verify-arrow` flag

---

### CHECKPOINT: simulation-level verification

Run a real simulation (criticality or fixed-source) with both data sources.
Same seed → same particle histories → same tallies.  If this passes, the Arrow
format is validated.

**This is the go/no-go for Phase 2.**

---

### PR 8 — (Phase 2) Remove HDF5 from simulation path
**Repo:** `yamc`
**Risk:** Medium — removes code, but behind the verification above

This only happens after PR 7's checkpoint passes.  Consider keeping both formats
for one release cycle before merging this.

- Delete `nuclide_hdf5.rs` (~1500 lines)
- Delete `photon_hdf5.rs` (~800 lines)
- Remove `hdf5-metno` from default deps (keep it for `yamc-precompute` binary
  behind a feature gate)
- Remove `NuclideFormat::Hdf5` variant from `nuclide_loader.rs`
- Make `arrow` feature non-optional (or the default)

**What to check:** `cargo build` works.  All tests pass.
`cargo build --features hdf5 --bin yamc_precompute` still compiles the
conversion tool.

**Files changed:**
- Delete `src/nuclide_hdf5.rs`
- Delete `src/photon_hdf5.rs`
- `Cargo.toml` — restructure features
- `src/nuclide_loader.rs` — remove HDF5 variant
- `src/lib.rs` — remove HDF5 module declarations

---

### PR 9 (optional, future) — Direct Arrow output from conversion pipeline
**Repo:** `nuclear_data_to_yamc_format`
**Risk:** Low — additive feature

Extend the Python conversion scripts to call `yamc-precompute` as a subprocess
(or reimplement the write in Python, accepting that it would need the exact same
post-processing logic).  This lets new users go from ENDF → `.yamc.arrow/`
without ever touching HDF5.

This is a future optimisation, not part of the core plan.

---

### Summary table

| PR | Repo | Description | Risk | Depends on |
|---|---|---|---|---|
| 1 | `yamc` | Refactor: pub `synthesize_hierarchical_mts` | None | — |
| 2 | `yamc` | Add `assert_nuclides_identical` | None | — |
| 3 | `yamc` | Add `arrow` feature + deps + stubs | None | — |
| 4 | `yamc` | Arrow writer + `yamc-precompute` binary | Low | 1, 3 |
| — | — | **Checkpoint: round-trip Arrow validation** | — | 4 |
| 5 | `yamc` | Arrow reader | Low | 3, 4 |
| — | — | **Checkpoint: field-by-field HDF5 vs Arrow** | — | 2, 5 |
| 6 | `yamc` | Wire into `nuclide_loader.rs` | Low | 5 |
| 7 | `yamc` | Verification mode + simulation tests | None | 2, 6 |
| — | — | **Checkpoint: simulation-level verification** | — | 7 |
| 8 | `yamc` | Phase 2: remove HDF5 from simulation path | Medium | 7 ✅ |
| 9 | `nuclear_data_to_yamc_format` | Direct Arrow output (future) | Low | 4 |

PRs 1, 2, and 3 have no dependencies on each other and can be developed and
merged in parallel.

---

## GPU relevance

### Does Arrow help port yamc to GPUs?

**Yes — Arrow is a significantly better starting point for GPU work than HDF5.**

The core problem for GPU Monte Carlo is getting nuclear data into GPU memory
efficiently.  With HDF5, the path is:

```
.h5 → CPU HDF5 C library → parse → build structs → copy each field to GPU buffers
```

Every struct field, every reaction's cross-section array, every distribution
table must be individually marshalled into a GPU-compatible layout.  HDF5 has no
concept of GPU memory and provides no tools for this.

With the simulation-ready Arrow format, the path becomes:

```
.yamc.arrow/ → mmap Arrow files → Arrow buffers (contiguous, aligned) → GPU upload
```

Arrow buffers are already contiguous byte arrays with known alignment.  This is
exactly what GPU APIs (CUDA, Vulkan, WebGPU) want: a flat chunk of memory with a
known layout that can be bulk-transferred via DMA or `cudaMemcpy`.

### What makes Arrow GPU-friendly

**Contiguous memory layout.**  Each Arrow column (e.g. all cross-section values
for all reactions at a given temperature) is a single contiguous allocation.
GPU bulk transfers work on contiguous ranges — no gathering from scattered heap
allocations.

**Fixed-width types with known alignment.**  Arrow's `Float64Array`,
`FixedSizeList<f64, 4>`, `UInt32Array` etc. have well-defined byte layouts that
map directly to GPU buffer formats.  The `FastXSGrid.xs` column
(`FixedSizeList<f64, 4>`) could be uploaded as-is and accessed as `double4` or
`dvec4` on the GPU.

**Zero-copy potential.**  On systems with unified memory (Apple Silicon, some
AMD APUs) or with GPU-direct storage (NVIDIA GPUDirect), Arrow's mmap'd buffers
could potentially be made visible to the GPU without any CPU-side copy at all.
HDF5 has no path to this — its data must be parsed before it can be used.

**Ecosystem.**  NVIDIA's RAPIDS (cuDF) and Apache Arrow GPU libraries already
handle Arrow ↔ GPU transfers.  If yamc ever needed to do cross-section lookups
in bulk on a GPU, the Arrow buffers could be passed directly to these libraries.

### What Arrow doesn't solve for GPU — and what we can do about it

yamc's transport loop is **history-based**: each particle is tracked through its
entire life (collisions, secondaries, death) before the next particle starts.
This is natural for CPUs but causes **warp divergence** on GPUs — threads in the
same warp take different branches (absorption vs scattering vs fission) and
execute at the speed of the slowest path.

The standard GPU fix is **event-based transport**: batch many particles together
and process them by event type.  All "sample reaction type" events run as one
kernel, all "sample scattering reaction" events as another.  This is an
architectural change to yamc's transport loop, independent of the data format.

But the data format *can* make the GPU kernels faster or slower.  Here's what
the current `FastXSGrid` layout does well and where it can be improved.

#### What already works well

**The log-bin energy lookup is GPU-friendly.**  `lookup_grid_index()` does:
`log(E)` → integer bin → `log_grid_index[bin]` → narrow search (1-2
iterations).  Every thread in a warp does the same number of operations — no
divergence.  The `log_grid_index` array (8192 entries, ~64 KB) fits comfortably
in GPU L2 cache.

**The `xs` array is compact.**  `Vec<[f64; 4]>` stores `[total, absorption,
scattering, fission]` interleaved at each energy point.  Two adjacent energy
points (needed for interpolation) are 64 bytes apart — one cache line on most
GPUs.  This is reasonable for random-access lookups.

#### What can be improved for GPU

**1. SoA layout for `FastXSGrid.xs`**

The current `FixedSizeList<f64, 4>` (AoS) stores `[total, abs, scat, fiss]`
together.  In event-based GPU transport, a warp of 32 threads all doing "sample
reaction type" for the same nuclide each need `total_xs` at their own energy
index.  With AoS, these reads have stride 32 bytes — not ideal for coalesced
GPU memory access.

Storing four separate columns (`xs_total`, `xs_absorption`, `xs_scattering`,
`xs_fission`) as plain `Float64Array`s gives stride 8 bytes — perfectly
coalesced.  All 32 threads reading `xs_total` at nearby energy indices hit
sequential memory addresses.

The Arrow format can store **both layouts** at negligible cost:

```
fast_xs.arrow:
  xs          – FixedSizeList<f64, 4>   (AoS, for CPU path)
  xs_total    – List<f64>               (SoA, for GPU path)
  xs_abs      – List<f64>
  xs_scat     – List<f64>
  xs_fiss     – List<f64>
```

The CPU reader uses `xs` (current behaviour, no change).  A future GPU loader
uploads `xs_total` etc. as separate GPU buffers.  The extra disk cost is ~2x for
this one table — ~5 MB per nuclide, well within the LZ4-compressed budget.

**2. Dense 2D matrix for `scatter_mt_xs`**

This is the biggest GPU pain point.  `sample_scatter_reaction()` iterates over
10-50 scattering reactions, doing linear interpolation for each, accumulating
cross-sections, then sampling with a random number.  On a GPU, this is a
**variable-length loop with early exit** — the worst case for warp efficiency.

The fix: pack `scatter_mt_xs` into a **dense 2D matrix**
`[n_energy_points × n_scatter_reactions]` stored as a single contiguous
`Float64Array` with a shape column.  This is exactly what a GPU texture or SSBO
(Shader Storage Buffer Object) wants.

```
fast_xs.arrow:
  scatter_mt_matrix  – Float64Array (flat, row-major: n_energy × n_scatter)
  scatter_mt_shape   – [n_energy, n_scatter]
  scatter_mt_ids     – List<i32> (MT numbers, length = n_scatter)
```

On the GPU, reaction sampling becomes a single read of one row of the matrix
(all reaction XS at the interpolated energy) followed by a prefix sum and
binary search — no variable-length loop, no early exit, deterministic execution.

**3. Precomputed scatter/fission CDFs per energy point**

Take the dense matrix one step further: precompute the cumulative XS at each
energy point.  Store a `[n_energy × n_scatter_reactions]` CDF matrix where
entry `[i, j]` = sum of XS for reactions 0..j at energy point i.

Sampling then becomes: look up energy index, read one row of the CDF matrix,
binary search for `xi * total_scatter_xs`.  Binary search on a fixed-size row
(10-50 entries) is O(log N) with no warp divergence — every thread does the
same number of comparisons.

This CDF matrix can be precomputed by `yamc-precompute` alongside the existing
data.  It's derived from `scatter_mt_xs` (a cumulative sum along the reaction
axis) so it adds no new physics logic.

**4. GPU texture mapping for XS interpolation**

GPU texture units do hardware-accelerated linear interpolation for free.  If
cross-section data is stored as a 1D texture with energy as the texture
coordinate, the GPU does the `xs0 + f * (xs1 - xs0)` interpolation in hardware
— zero ALU cost.

The Arrow format enables this because the XS arrays are contiguous float
buffers that map directly to texture memory.  With HDF5, you'd have to parse
the data into a separate buffer first.

The log-bin index would become a texture coordinate transform:
`tex_coord = (log(E) - log_e_min) * inv_log_delta / n_energy_points`.  The
GPU's texture unit handles the rest.

#### Summary of GPU-ready additions to the Arrow schema

| Addition | Where | Disk cost | GPU benefit |
|---|---|---|---|
| SoA XS columns | `fast_xs.arrow` | ~5 MB/nuclide | Coalesced warp memory access |
| Dense scatter matrix | `fast_xs.arrow` | ~0 (reformat of existing data) | Eliminates variable-length loop |
| Scatter CDF matrix | `fast_xs.arrow` | ~same as scatter matrix | O(log N) sampling, no divergence |
| Fission CDF matrix | `fast_xs.arrow` | small (1-5 reactions) | Same benefit for fission sampling |

These are **optional additional columns** in the same Arrow files.  The CPU
reader ignores them.  They cost minimal extra disk space (mostly zeros that
compress well with LZ4).  They can be added in a follow-up PR after the core
Arrow format is working — no need to delay the initial implementation.

### Comparison: HDF5 vs Arrow for GPU porting

| Aspect | HDF5 | Arrow |
|---|---|---|
| Data loading | C library, CPU only, no GPU awareness | Contiguous buffers, mmap-able, GPU-transfer-ready |
| Bulk upload to GPU | Must copy field-by-field from scattered heap allocations | Can upload whole columns as contiguous byte ranges |
| SoA / AoS flexibility | Would require writing a second HDF5 dataset | Extra columns in the same Arrow file |
| Dense 2D XS matrices | Not stored — must build on CPU, then upload | Stored contiguously, upload directly |
| Zero-copy GPU access | Not possible — must parse first | Possible on unified memory architectures |
| GPU texture mapping | Need to parse then copy to texture buffer | Contiguous float arrays map directly to textures |
| GPU ecosystem support | None | RAPIDS cuDF, Arrow GPU libraries, Vulkano interop |
| Pre-computed lookup tables | Not stored — must build on CPU | Stored in `.yamc.arrow/`, upload directly |
| Cross-platform GPU | N/A | Arrow byte layout is platform-independent — same file works for CUDA, Vulkan, WebGPU |
| WebGPU (browser) | HDF5 has no WASM story | Arrow IPC can be fetched + parsed in WASM, buffers passed to WebGPU |

### Bottom line

Arrow does not magically make yamc run on GPUs.  The transport loop needs to
become event-based, the physics kernels need GPU implementations, and geometry
navigation is its own challenge.  But Arrow lets us **shape the data for GPU
access patterns at precompute time** rather than at simulation time:

- SoA columns give coalesced memory access for GPU warps
- Dense 2D matrices and precomputed CDFs eliminate the variable-length reaction sampling loop — the single worst GPU performance killer in the current code
- Contiguous aligned buffers can be uploaded (or zero-copied) to GPU memory without transformation
- GPU texture units can do XS interpolation in hardware, for free, if the data is in contiguous float arrays — which Arrow gives us
- The WebGPU path benefits doubly: Arrow IPC can be fetched over HTTP and parsed in WASM, HDF5 has no browser story
- These GPU-oriented columns are **optional extras** in the same files — they don't complicate the CPU path and can be added incrementally
