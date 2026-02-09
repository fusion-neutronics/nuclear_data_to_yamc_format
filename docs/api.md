# API reference

## Conversion functions

High-level entry points that accept source files (ACE or ENDF) and write
simulation-ready `.arrow/` output.

```{autofunction} nuclear_data_to_yamc_format.convert_neutron
```

```{autofunction} nuclear_data_to_yamc_format.convert_photon
```

## Low-level writers

Accept OpenMC Python objects directly.  Useful when you have already loaded
data through your own pipeline.

```{autofunction} nuclear_data_to_yamc_format.neutron_writer.export_neutron_to_arrow
```

```{autofunction} nuclear_data_to_yamc_format.photon_writer.export_photon_to_arrow
```

## Readers

Read `.arrow/` directories back into Python dicts for inspection or
verification.

```{autofunction} nuclear_data_to_yamc_format.neutron_reader.read_neutron_from_arrow
```

```{autofunction} nuclear_data_to_yamc_format.photon_reader.read_photon_from_arrow
```

## Verification

Compare an OpenMC object against its Arrow export.  Uses tolerance-based
comparison (`np.allclose` with `rtol=1e-12`) and verifies synthesized MTs
and FastXSGrid consistency.

```{autofunction} nuclear_data_to_yamc_format.verify.verify_neutron
```

```{autofunction} nuclear_data_to_yamc_format.verify.verify_photon
```

## Synthesis module

MT synthesis and FastXSGrid construction algorithms.

```{autofunction} nuclear_data_to_yamc_format.synthesis.is_scattering_mt
```

```{autofunction} nuclear_data_to_yamc_format.synthesis.is_fission_mt
```

```{autofunction} nuclear_data_to_yamc_format.synthesis.synthesize_hierarchical_mts
```

```{autofunction} nuclear_data_to_yamc_format.synthesis.build_fast_xs
```

## Constants

```{autodata} nuclear_data_to_yamc_format.synthesis.N_LOG_BINS
```

```{autodata} nuclear_data_to_yamc_format.synthesis.SYNTHETIC_MTS
```

```{autodata} nuclear_data_to_yamc_format.synthesis.ALL_SCATTERING_MTS
```

```{autodata} nuclear_data_to_yamc_format.synthesis.FISSION_MTS
```

```{autodata} nuclear_data_to_yamc_format.synthesis.ABSORPTION_MTS
```
