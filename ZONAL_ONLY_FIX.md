# Zonal-only execution patch

`write_rasters = FALSE` is supported when either `zonal_path` or `points_path` is provided. Raster writing is disabled without changing the modeled domain for zonal runs. Monthly climate output collection remains active for zonal CSVs, and point-only mode retains its reduced-cell optimization.

Example:

```r
cfg <- wb_config()
cfg$write_rasters <- FALSE
cfg$zonal_path <- "zones.gpkg"
cfg$zonal_id_field <- "zone_id"
wb_run(cfg)
```

No scientific equations were modified. This is a source patch; runtime testing is still required.
