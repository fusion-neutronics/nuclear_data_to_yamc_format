# API reference

## Conversion functions

These are the main entry points.  They accept source files (ACE or ENDF) and
write Arrow output.

```{autofunction} nuclear_data_to_yamc_format.convert_neutron
```

```{autofunction} nuclear_data_to_yamc_format.convert_photon
```

## Low-level writers

These accept OpenMC Python objects directly.  Useful when you have already
loaded data through your own pipeline.

```{autofunction} nuclear_data_to_yamc_format.neutron_writer.export_neutron_to_arrow
```

```{autofunction} nuclear_data_to_yamc_format.photon_writer.export_photon_to_arrow
```

## Readers

Read Arrow directories back into Python dicts for inspection or verification.

```{autofunction} nuclear_data_to_yamc_format.neutron_reader.read_neutron_from_arrow
```

```{autofunction} nuclear_data_to_yamc_format.photon_reader.read_photon_from_arrow
```

## Verification

Compare an OpenMC object against its Arrow export to confirm bit-identical
output.

```{autofunction} nuclear_data_to_yamc_format.verify.verify_neutron
```

```{autofunction} nuclear_data_to_yamc_format.verify.verify_photon
```
