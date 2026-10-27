# Persistent zonal mapping cache (development patch)

`wb_run()` now caches the mapping from model cells to user-defined zones for local raster or vector files. No new configuration setting is required. The cache is stored under `<environment_cache>/zonal_mappings/` and is reused when the zone file signature, identifier field, model raster grid, and ordered modeled-cell indices agree.

The mapping is computed once per new zone/grid combination. It is already reused across years within a run; this change makes it reusable across **separate runs**. It does not optimize repeated monthly quantile grouping, nor does it change the model equations. In-memory `SpatRaster`/`SpatVector` objects are prepared each run.

Look for `Cached zonal cell mapping:` on the first run and `Reusing zonal cell mapping:` on subsequent runs. Compare zonal CSVs against the prior version and benchmark separately for first-run preparation and later-run reuse. No R runtime tests were available in this build environment.
