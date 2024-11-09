# GML ozonesonde and UFS-AQM

[GML ozonesondes](https://gml.noaa.gov/ozwv/ozsondes/) can be fetched and loaded with MONETIO.
For this example, we use a pre-prepared dataset of the 100-m data[^100m] and compare to
runs of developmental versions of the UFS-AQM regional model.

This example works with {doc}`the CLI </cli>`. (``melodies-monet run``).

````{admonition} Control file.
:class: dropdown

```{literalinclude} control_ufsaqm_ozonesonde.yaml
:caption:
:linenos:
```
````

[^100m]: Soundings regridded into a mostly 100-m vertical grid, with a single time assigned
    for the whole profile.

## Plots

```{figure} output/ufsaqm_ozonesonde/plot_grp1.vertical_single_date.o3.2023-06-24_00.2023-06-25_00.all.CONUS.png
:width: 85 %
:alt: vertical profiles

```

```{figure} output/ufsaqm_ozonesonde/plot_grp2.vertical_boxplot_os.o3.2023-06-24_00.2023-06-25_00.all.CONUS.png
:width: 85 %
:alt: vertical profiles with boxplots

````

```{figure} output/ufsaqm_ozonesonde/plot_grp3.density_scatter_plot_os.o3.2023-06-24_00.2023-06-25_00.all.CONUS.gml-ozonesondes_ufsaqm_cmaq52.ground-level.png
:width: 85 %
:alt: comparison scatter plot

```

```{figure} output/ufsaqm_ozonesonde/plot_grp3.density_scatter_plot_os.o3.2023-06-24_00.2023-06-25_00.all.CONUS.gml-ozonesondes_ufsaqm_cmaq54.ground-level.png
:width: 85 %
:alt: comparison scatter plot for the other model

```
