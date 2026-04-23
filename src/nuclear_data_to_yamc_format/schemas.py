"""Arrow schema definitions for simulation-ready nuclear data tables."""

import pyarrow as pa


# ---------------------------------------------------------------------------
# Neutron schemas
# ---------------------------------------------------------------------------

NUCLIDE_SCHEMA = pa.schema(
    [
        pa.field("name", pa.utf8()),
        pa.field("Z", pa.int32()),
        pa.field("A", pa.int32()),
        pa.field("metastable", pa.int32()),
        pa.field("atomic_weight_ratio", pa.float64()),
        pa.field("temperatures", pa.list_(pa.utf8())),
        pa.field("kTs", pa.list_(pa.float64())),
        pa.field("energy_temperatures", pa.list_(pa.utf8())),
        pa.field("energy_values", pa.list_(pa.list_(pa.float64()))),
        # Pre-computed at build time: true iff any fission MT (18, 19, 20, 21, 38)
        # is present. Nullable for back-compat with older files.
        pa.field("fissionable", pa.bool_(), nullable=True),
    ],
    metadata={b"filetype": b"data_neutron", b"version": b"4.0"},
)

REACTIONS_SCHEMA = pa.schema(
    [
        pa.field("mt", pa.int32()),
        pa.field("label", pa.utf8()),
        pa.field("Q_value", pa.float64()),
        pa.field("center_of_mass", pa.bool_()),
        pa.field("redundant", pa.bool_()),
        pa.field("xs_temperatures", pa.list_(pa.utf8())),
        pa.field("xs_values", pa.list_(pa.list_(pa.float64()))),
        pa.field("xs_threshold_idx", pa.list_(pa.int32())),
        pa.field("n_products", pa.int32()),
    ]
)

PRODUCTS_SCHEMA = pa.schema(
    [
        pa.field("reaction_mt", pa.int32()),
        pa.field("product_idx", pa.int32()),
        pa.field("particle", pa.utf8()),
        pa.field("emission_mode", pa.utf8()),
        pa.field("decay_rate", pa.float64()),
        pa.field("n_distribution", pa.int32()),
        pa.field("yield_type", pa.utf8()),
        pa.field("yield_data", pa.list_(pa.float64())),
        pa.field("yield_shape", pa.list_(pa.int32())),
        pa.field("yield_breakpoints", pa.list_(pa.int32())),
        pa.field("yield_interpolation", pa.list_(pa.int32())),
    ]
)

DISTRIBUTIONS_SCHEMA = pa.schema(
    [
        pa.field("reaction_mt", pa.int32()),
        pa.field("product_idx", pa.int32()),
        pa.field("dist_idx", pa.int32()),
        pa.field("type", pa.utf8()),
        # Applicability (nullable)
        pa.field("applicability_data", pa.list_(pa.float64())),
        pa.field("applicability_shape", pa.list_(pa.int32())),
        pa.field("applicability_breakpoints", pa.list_(pa.int32())),
        pa.field("applicability_interpolation", pa.list_(pa.int32())),
        # --- Angle distribution (nullable) ---
        pa.field("angle_energies", pa.list_(pa.float64())),
        pa.field("angle_mu_data", pa.list_(pa.float64())),
        pa.field("angle_mu_offsets", pa.list_(pa.int32())),
        pa.field("angle_mu_interpolation", pa.list_(pa.int32())),
        # --- Energy distribution (nullable) ---
        pa.field("energy_dist_type", pa.utf8()),
        pa.field("energy_dist_energies", pa.list_(pa.float64())),
        pa.field("energy_dist_interpolation", pa.list_(pa.int32())),
        pa.field("energy_dist_data", pa.list_(pa.float64())),
        pa.field("energy_dist_offsets", pa.list_(pa.int32())),
        pa.field("energy_dist_out_interp", pa.list_(pa.int32())),
        pa.field("energy_dist_n_discrete", pa.list_(pa.int32())),
        pa.field("energy_param_x", pa.list_(pa.float64())),
        pa.field("energy_param_y", pa.list_(pa.float64())),
        pa.field("energy_param2_x", pa.list_(pa.float64())),
        pa.field("energy_param2_y", pa.list_(pa.float64())),
        pa.field("energy_restriction_u", pa.float64()),
        pa.field("energy_threshold", pa.float64()),
        pa.field("energy_mass_ratio", pa.float64()),
        pa.field("energy_primary_flag", pa.int32()),
        pa.field("energy_atomic_weight_ratio", pa.float64()),
        pa.field("energy_discrete_energy", pa.float64()),
        # --- Correlated columns (nullable) ---
        pa.field("corr_energies", pa.list_(pa.float64())),
        pa.field("corr_breakpoints", pa.list_(pa.int32())),
        pa.field("corr_interpolation", pa.list_(pa.int32())),
        pa.field("corr_eout_data", pa.list_(pa.float64())),
        pa.field("corr_eout_offsets", pa.list_(pa.int32())),
        pa.field("corr_eout_interp", pa.list_(pa.int32())),
        pa.field("corr_eout_n_discrete", pa.list_(pa.int32())),
        pa.field("corr_mu_data", pa.list_(pa.float64())),
        pa.field("corr_mu_offsets", pa.list_(pa.int32())),
        pa.field("corr_mu_interp", pa.list_(pa.int32())),
        # --- Kalbach-Mann columns (nullable) ---
        pa.field("km_energies", pa.list_(pa.float64())),
        pa.field("km_breakpoints", pa.list_(pa.int32())),
        pa.field("km_interpolation", pa.list_(pa.int32())),
        pa.field("km_data", pa.list_(pa.float64())),
        pa.field("km_offsets", pa.list_(pa.int32())),
        pa.field("km_interp", pa.list_(pa.int32())),
        pa.field("km_n_discrete", pa.list_(pa.int32())),
        # --- N-body columns (nullable) ---
        pa.field("nbody_n", pa.int32()),
        pa.field("nbody_total_mass", pa.float64()),
        pa.field("nbody_atomic_weight_ratio", pa.float64()),
        pa.field("nbody_q_value", pa.float64()),
    ]
)

URR_SCHEMA = pa.schema(
    [
        pa.field("temperature", pa.utf8()),
        pa.field("energy", pa.list_(pa.float64())),
        pa.field("table_data", pa.list_(pa.float64())),
        pa.field("table_shape", pa.list_(pa.int32())),
        pa.field("interpolation", pa.int32()),
        pa.field("inelastic", pa.int32()),
        pa.field("absorption", pa.int32()),
        pa.field("multiply_smooth", pa.bool_()),
    ]
)

TOTAL_NU_SCHEMA = pa.schema(
    [
        pa.field("particle", pa.utf8()),
        pa.field("emission_mode", pa.utf8()),
        pa.field("decay_rate", pa.float64()),
        pa.field("yield_type", pa.utf8()),
        pa.field("yield_data", pa.list_(pa.float64())),
        pa.field("yield_shape", pa.list_(pa.int32())),
        pa.field("yield_breakpoints", pa.list_(pa.int32())),
        pa.field("yield_interpolation", pa.list_(pa.int32())),
    ]
)

# ---------------------------------------------------------------------------
# FastXSGrid schema — one row per temperature
# ---------------------------------------------------------------------------

#: Fixed size of the dense MT → index lookup tables. Large enough for every
#: ENDF-defined MT number in practice (levels run up to ~890).
MT_LOOKUP_SIZE = 1024

FAST_XS_SCHEMA = pa.schema(
    [
        pa.field("temperature", pa.utf8()),
        pa.field("log_e_min", pa.float64()),
        pa.field("inv_log_delta", pa.float64()),
        pa.field("log_grid_index", pa.list_(pa.int32())),
        pa.field("xs", pa.list_(pa.float64())),
        pa.field("xs_shape", pa.list_(pa.int32())),
        pa.field("energy", pa.list_(pa.float64())),
        pa.field("scatter_mt_numbers", pa.list_(pa.int32())),
        pa.field("scatter_mt_xs", pa.list_(pa.float64())),
        pa.field("scatter_mt_shape", pa.list_(pa.int32())),
        pa.field("fission_mt_numbers", pa.list_(pa.int32())),
        pa.field("fission_mt_xs", pa.list_(pa.float64())),
        pa.field("fission_mt_shape", pa.list_(pa.int32())),
        pa.field("has_partial_fission", pa.bool_()),
        pa.field("xs_ngamma", pa.list_(pa.float64())),
        pa.field("photon_prod", pa.list_(pa.float64())),
        # Pre-computed scalars (derived from the shape arrays above).
        # Explicit for GPU buffer sizing and faster load. All nullable for
        # back-compat with older files.
        pa.field("n_energies", pa.int32(), nullable=True),
        pa.field("n_scatter_mts", pa.int32(), nullable=True),
        pa.field("n_fission_mts", pa.int32(), nullable=True),
        # Dense MT → index lookup tables of length MT_LOOKUP_SIZE.
        # scatter_mt_to_idx[mt] = position in `scatter_mt_numbers`, or -1 if
        # absent. Removes the linear scan that would otherwise run per MT
        # query (elastic, (n,2n), etc.). Same shape for fission.
        pa.field("scatter_mt_to_idx", pa.list_(pa.int32()), nullable=True),
        pa.field("fission_mt_to_idx", pa.list_(pa.int32()), nullable=True),
    ]
)

# ---------------------------------------------------------------------------
# Photon schemas
# ---------------------------------------------------------------------------

ELEMENT_SCHEMA = pa.schema(
    [
        pa.field("name", pa.utf8()),
        pa.field("Z", pa.int32()),
        pa.field("energy", pa.list_(pa.float64())),
        pa.field("ln_energy", pa.list_(pa.float64())),
        pa.field("coherent_xs", pa.list_(pa.float64())),
        pa.field("ln_coherent_xs", pa.list_(pa.float64())),
        pa.field("incoherent_xs", pa.list_(pa.float64())),
        pa.field("ln_incoherent_xs", pa.list_(pa.float64())),
        pa.field("photoelectric_xs", pa.list_(pa.float64())),
        pa.field("ln_photoelectric_xs", pa.list_(pa.float64())),
        pa.field("pair_production_nuclear_xs", pa.list_(pa.float64())),
        pa.field("pair_production_electron_xs", pa.list_(pa.float64())),
        pa.field("heating_xs", pa.list_(pa.float64())),
        pa.field("coherent_int_ff_x", pa.list_(pa.float64())),
        pa.field("coherent_int_ff_y", pa.list_(pa.float64())),
        pa.field("coherent_ff_x", pa.list_(pa.float64())),
        pa.field("coherent_ff_y", pa.list_(pa.float64())),
        pa.field("coherent_anomalous_real_x", pa.list_(pa.float64())),
        pa.field("coherent_anomalous_real_y", pa.list_(pa.float64())),
        pa.field("coherent_anomalous_imag_x", pa.list_(pa.float64())),
        pa.field("coherent_anomalous_imag_y", pa.list_(pa.float64())),
        pa.field("incoherent_ff_x", pa.list_(pa.float64())),
        pa.field("incoherent_ff_y", pa.list_(pa.float64())),
    ],
    metadata={b"filetype": b"data_photon", b"version": b"4.0"},
)

SUBSHELLS_SCHEMA = pa.schema(
    [
        pa.field("designator", pa.utf8()),
        pa.field("binding_energy", pa.float64()),
        pa.field("num_electrons", pa.float64()),
        pa.field("xs", pa.list_(pa.float64())),
        pa.field("ln_xs", pa.list_(pa.float64())),
        pa.field("threshold_idx", pa.int32()),
        pa.field("transitions_data", pa.list_(pa.float64())),
        pa.field("transitions_shape", pa.list_(pa.int32())),
    ]
)

COMPTON_SCHEMA = pa.schema(
    [
        pa.field("num_electrons", pa.list_(pa.float64())),
        pa.field("binding_energy", pa.list_(pa.float64())),
        pa.field("pz", pa.list_(pa.float64())),
        pa.field("J_data", pa.list_(pa.float64())),
        pa.field("J_shape", pa.list_(pa.int32())),
        pa.field("J_cdf_data", pa.list_(pa.float64())),
        pa.field("J_cdf_shape", pa.list_(pa.int32())),
    ]
)

BREMSSTRAHLUNG_SCHEMA = pa.schema(
    [
        pa.field("I", pa.float64()),
        pa.field("electron_energy", pa.list_(pa.float64())),
        pa.field("photon_energy", pa.list_(pa.float64())),
        pa.field("num_electrons", pa.list_(pa.float64())),
        pa.field("ionization_energy", pa.list_(pa.float64())),
        pa.field("dcs_data", pa.list_(pa.float64())),
        pa.field("dcs_shape", pa.list_(pa.int32())),
    ]
)

# ---------------------------------------------------------------------------
# Depletion chain schemas
# ---------------------------------------------------------------------------

CHAIN_NUCLIDES_SCHEMA = pa.schema(
    [
        pa.field("name", pa.utf8(), nullable=False),
        pa.field("half_life", pa.float64(), nullable=True),
        pa.field("decay_energy", pa.float64(), nullable=False),
        pa.field("n_decay_modes", pa.int32(), nullable=False),
        pa.field("n_reactions", pa.int32(), nullable=False),
        pa.field("n_sources", pa.int32(), nullable=False),
        pa.field("has_fission_yields", pa.bool_(), nullable=False),
        pa.field("fission_yield_parent", pa.utf8(), nullable=True),
    ],
    metadata={b"filetype": b"depletion_chain", b"version": b"1.0"},
)

CHAIN_DECAYS_SCHEMA = pa.schema(
    [
        pa.field("nuclide", pa.utf8(), nullable=False),
        pa.field("type", pa.utf8(), nullable=False),
        pa.field("target", pa.utf8(), nullable=True),
        pa.field("branching_ratio", pa.float64(), nullable=False),
    ]
)

CHAIN_REACTIONS_SCHEMA = pa.schema(
    [
        pa.field("nuclide", pa.utf8(), nullable=False),
        pa.field("type", pa.utf8(), nullable=False),
        pa.field("target", pa.utf8(), nullable=True),
        pa.field("Q", pa.float64(), nullable=False),
        pa.field("branching_ratio", pa.float64(), nullable=False),
    ]
)

CHAIN_SOURCES_SCHEMA = pa.schema(
    [
        pa.field("nuclide", pa.utf8(), nullable=False),
        pa.field("particle", pa.utf8(), nullable=False),
        pa.field("type", pa.utf8(), nullable=False),
        pa.field("energies", pa.list_(pa.float64()), nullable=False),
        pa.field("intensities", pa.list_(pa.float64()), nullable=False),
    ]
)

CHAIN_FISSION_YIELDS_SCHEMA = pa.schema(
    [
        pa.field("nuclide", pa.utf8(), nullable=False),
        pa.field("energy", pa.float64(), nullable=False),
        pa.field("products", pa.list_(pa.utf8()), nullable=False),
        pa.field("yields", pa.list_(pa.float64()), nullable=False),
    ]
)
