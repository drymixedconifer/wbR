# Latitude transform fix

`wb_latitude_raster()` now transforms cell centers only where the DEM has finite elevation. NA cells in the rectangular raster extent are kept as NA. Invalid transformations of valid terrain cells still raise a descriptive CRS error. No scientific water-balance equations changed.

Validation needed in R: run initial preparation on the failing geography, then verify modeled latitude range and compare with `terra::project()` for selected valid cells.
