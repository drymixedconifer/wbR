# v6.3 integration candidate (not yet validated)

* Replaced both package engine snapshots with the uploaded standalone v6.3 engine (same bytes).
* Added missing soil-depth, soil-provenance, constant-AWC and PRISM-residency configuration defaults.
* Regenerated R/config_defaults.R from inst/config/defaults.R.
* Retained wb_config(), wb_run(), wb_engine_info() and the runtime source-loading architecture.
* Removed an interface check blocking the engine-supported `analysis_initialization` parameter; the engine itself validates chronology.
* Updated dev/sync_engine.R paths, previously pointing at nonexistent v5.5 installed snapshot.
* Future-specific code is NOT in this public candidate.

## Required on a Mac with R installed

```r
install.packages(c("terra","prism","Rcpp","RcppParallel","testthat"))
install.packages(".", repos = NULL, type = "source") # run from extracted public package directory
library(wbR)
wb_engine_info()
cfg <- wb_config(advanced = TRUE)
str(cfg[c("awc_depth_path", "awc_resample", "prism_check_resident")])
```

Point cfg paths at a **small** AOI and set `cfg$stop_after_setup <- TRUE` for the first validation.
Then use a small multi-year historical test. No runtime or data validation was possible in the build environment.
