# DEM latitude reader correction

`wb_latitude_raster()` now reads its input DEM directly rather than using `terra::rast(dem)` as its source. `rast(dem)` may return a geometry-only raster; calling `readValues` on it then yields no valid elevations. The output latitude raster is still created from a geometry-only template, which is appropriate for writing. No water-balance equations changed.

Runtime validation in R is required.
