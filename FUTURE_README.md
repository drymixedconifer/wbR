# Private wbR future branch candidate

This is a full installable copy of wbR with one additional exported function, `wbR_future()`.
The v6.3 engine in `inst/engine/water_balance_v6_functions.R` is byte-for-byte the same as in the public candidate.

Do not publish this archive or push this branch to a public remote until ready.

```r
# from the private package directory, install it
install.packages(".", repos = NULL, type = "source")
library(wbR)
cfg <- wb_config(advanced = TRUE)
cfg$dem_path <- "/path/to/dem.tif"
cfg$awc_path <- "/path/to/awc.tif"
cfg$aoi_path <- "/path/to/aoi.shp"
cfg$prism_root <- "/path/to/prism"
cfg$output_root <- "/path/to/output"
cfg$base_years <- 1991:2020
cfg$years_annual <- 2021:2074
cfg$analysis_initialization <- "baseline_end"
cfg$months_out <- 1:12
cfg$outputs_to_write <- c("cwd", "aet")
cfg$monthly_output <- FALSE
cfg$write_rasters <- TRUE
cfg$stop_after_setup <- TRUE  # first test only
wbR_future(cfg, climate_dir = "/path/to/LOCA2",
           models = "INM-CM5-0", scenario = "ssp245", modes = "anchored")
```

Remove `cfg$stop_after_setup <- TRUE` to actually simulate. First use a small AOI and limited years consistent with your NetCDF coverage.

**Not runtime-tested.** The private climate extension was originally written for v6.2 and has been loaded onto the v6.3 engine in a static integration. Validate numerical results and restart behavior before operational West-wide runs.
