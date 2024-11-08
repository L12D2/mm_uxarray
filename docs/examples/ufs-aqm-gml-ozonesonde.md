# GML ozonesonde and UFS-AQM

[GML ozonesondes](https://gml.noaa.gov/ozwv/ozsondes/) can be fetched and loaded with MONETIO.
For this example, we use a pre-prepared dataset of the 100-m data [^100m] and compare to
runs of developmental versions of the UFS-AQM regional model.

This example works with {doc}`the CLI </cli>`. (``melodies-monet run``).

````{admonition} Control file.
:class: dropdown

```{literalinclude} UFS-AQM_ozone_sonder.yaml
:caption:
:linenos:
```
````

[^100m]: Soundings regridded into a mostly 100-m vertical grid, with a single time assigned
    for the whole profile.

## Plots
