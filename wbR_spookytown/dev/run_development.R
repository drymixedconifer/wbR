# From the wbR project root after installing the package once:
library(wbR)

# Configure the run without modifying the installed engine.
cfg <- wb_config(engine_path = "dev/water_balance_v5.5_functions.R")
cfg$base_dir <- "/path/to/project"
cfg$dem_path <- file.path(cfg$base_dir, "dem.tif")
cfg$aoi_path <- file.path(cfg$base_dir, "watershed.gpkg")
cfg$awc_path <- file.path(cfg$base_dir, "awc.tif")
cfg$prism_root <- file.path(cfg$base_dir, "PRISM")
cfg$output_root <- file.path(cfg$base_dir, "CWD_outputs_dev")
cfg$monthly_output <- FALSE
cfg$zonal_path <- file.path(cfg$base_dir, "zones.tif")
cfg$zonal_monthly <- c("cwd", "aet", "ppt", "pack", "melt", "soil")

# THIS runs the editable source, not the installed snapshot.
wb_run(cfg, engine_path = "dev/water_balance_v5.5_functions.R")
