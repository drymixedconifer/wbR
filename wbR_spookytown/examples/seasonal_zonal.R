library(wbR)
cfg <- wb_config()
# Set cfg$dem_path, cfg$aoi_path, cfg$awc_path, cfg$prism_root,
# cfg$base_dir and cfg$output_root to your existing valid paths.
cfg$outputs_to_write <- c("aet", "pet", "snow", "soil", "runoff")
cfg$summarize <- TRUE
cfg$monthly_output <- FALSE
cfg$monthly_years <- NULL
cfg$zonal_path <- "/path/to/zones.gpkg"
cfg$zonal_id_field <- "HUC12"
cfg$zonal_monthly <- c("cwd", "aet", "pet", "ppt", "soil", "pack", "runoff", "melt")
cfg$zonal_annual <- TRUE
# Requires the full domain to fit in available memory; no tiling.
wb_run(cfg)
