library(wbR)

cfg <- wb_config()
cfg$base_dir <- "/path/to/project"
cfg$dem_path <- file.path(cfg$base_dir, "dem.tif")
cfg$aoi_path <- file.path(cfg$base_dir, "watershed.gpkg")
cfg$awc_path <- file.path(cfg$base_dir, "awc.tif")
cfg$prism_root <- file.path(cfg$base_dir, "PRISM")
cfg$output_root <- file.path(cfg$base_dir, "CWD_outputs")

# Monthly rasters OFF; monthly zonal histories ON.
cfg$monthly_output <- FALSE
cfg$zonal_path <- file.path(cfg$base_dir, "huc12.gpkg")
cfg$zonal_id_field <- "HUC12"
cfg$zonal_monthly <- c("cwd", "aet", "pet", "ppt", "water",
                        "soil", "pack", "snowfall", "melt", "tmean")
cfg$zonal_annual <- TRUE
cfg$outputs_to_write <- c("cwd", "aet")
cfg$write_rasters <- TRUE  # required by v5.5 for full-grid runs

wb_run(cfg)
