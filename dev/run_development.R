# From the wbR project root after installing the package once:
library(wbR)

# Run the editable development engine (dev/water_balance_v6_functions.R) instead of the installed snapshot.
# After the edits are final, run dev/sync_engine.R to copy it into inst/engine/.
engine <- "dev/water_balance_v6_functions.R"

cfg <- wb_config(engine_path = engine)
cfg$base_dir <- "/path/to/project"
cfg$dem_path <- file.path(cfg$base_dir, "dem.tif")
cfg$aoi_path <- file.path(cfg$base_dir, "watershed.gpkg")
cfg$awc_path <- file.path(cfg$base_dir, "awc.tif")
cfg$prism_root <- file.path(cfg$base_dir, "PRISM")
cfg$output_root <- file.path(cfg$base_dir, "CWD_outputs_dev")
cfg$monthly_output <- FALSE

# THIS runs the editable source, not the installed snapshot.
wb_run(cfg, engine_path = engine)
