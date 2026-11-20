########################################################################
## ======== Functions for FUTURE water balance rasters (v6.2) ======== ##
########################################################################

# Companion to water_balance_future_v6.2.R. Keep it beside water_balance_v6.2_functions.R: the driver sources that file first and this one second, into the same environment. This file only defines functions; sourcing it runs nothing.

# What this does: The historical engine (terrain, soil, lapse rates, radiation, snow, Penman-Monteith, soil bookkeeping) is used exactly as water_balance_v6.2.R uses it. Only the monthly climate changes: instead of PRISM months, each simulated year reads 12 months of temperature, precipitation and vapour-pressure deficit from a rebuilt LOCA2 file (<MODEL>_<scenario>_monthly_1991-2074_west.nc, written by loca2_monthly_west_v2.py). Every model is run two ways:

#   "anchored"  Change from the PRISM starting point. Each LOCA2 month is turned into a departure from that model's own 1991-2020 mean for the calendar month (temperature: difference in C; precipitation and VPD: ratio), the departure is interpolated from the 6 km LOCA2 grid to the 800 m PRISM grid, and applied to the PRISM 1991-2020 normal. The 1991-2020 mean climate of the run is therefore the PRISM normal, and the maps are on the same footing as the 1950-2025 historical maps. Spin-up uses the PRISM normals.
#   "direct"    Absolute values as the climate model (after LOCA2 downscaling) gives them. LOCA2 months are fed to the engine on their own 6 km grid; elevation adjustment within each 6 km cell uses the PRISM-fitted lapse rates. Spin-up uses the model's own 1991-2020 monthly means. Change is measured against the same run's 1991-2020 mean.

# Unchanged from the historical run in both modes: cloud transmittance and dew-point normals (PRISM 1991-2020) and the wind climatology are static, so future change enters through temperature, precipitation and VPD only.

# The one function of the historical file that is replaced is wb_climate_year(): with a year and a site that carries
# LOCA2 information it returns LOCA2 months, otherwise it behaves exactly as before (PRISM normals).

fut_modes <- c("anchored", "direct")
fut_vars <- c(ppt = "pr_total", tmean = "tas_mean", vpd = "vpd")       # engine name = variable name in the LOCA2 file
fut_units <- c(pr_total = "mm", tas_mean = "degC", vpd = "kPa")

# ---------------- LOCA2 files ----------------

# The rebuilt files store longitude as degrees east plus 720 (594-619) and carry no coordinate system.
fut_geo <- function(x) {
  e <- as.vector(ext(x))
  if (e[3] < -90 || e[4] > 90) stop("LOCA2 file does not look like a longitude/latitude grid: ", sources(x)[1])
  e[1:2] <- e[1:2] - 360 * floor((mean(e[1:2]) + 180) / 360)
  ext(x) <- ext(e)
  if (!nzchar(crs(x))) crs(x) <- "EPSG:4326"
  x
}

fut_find_nc <- function(fut, model) {
  pattern <- paste0("^", gsub(".", "\\.", model, fixed = TRUE), "_", fut$scenario, "_monthly_[0-9]{4}-[0-9]{4}_west\\.nc$")
  f <- list.files(fut$climate_dir, pattern = pattern, full.names = TRUE)
  if (length(f) != 1L) stop("Expected exactly one LOCA2 file for ", model, " (", fut$scenario, ") in ", fut$climate_dir,
                            "; found ", length(f), ". Files are named <MODEL>_<scenario>_monthly_<first>-<last>_west.nc.")
  f
}

# Read a variable's NetCDF units attribute. terra::units() is not a SpatRaster method.
fut_nc_units <- function(path, variable) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  att <- ncdf4::ncatt_get(nc, variable, "units")
  if (!isTRUE(att$hasatt)) stop("Missing NetCDF units attribute for ", variable, ": ", path)
  as.character(att$value)
}

# Open the three variables, check units and the monthly calendar, put the grid on -180..180 and crop to the site.
fut_open <- function(path, model, site, cfg) {
  old <- wb_enter("Opening LOCA2 file", "not applicable", path); on.exit(.wb_log$context <- old, add = TRUE)
  vars <- fut_vars[c("ppt", "tmean", if (cfg$pet_method == "penman_monteith") "vpd")]
  h <- list(path = path, model = model)
  for (k in names(vars)) {
    r <- wb_step("Opening LOCA2 variable", path, rast(path, subds = vars[[k]]))
    u <- fut_nc_units(path, vars[[k]])
    if (length(u) != 1L || !identical(u, fut_units[[vars[[k]]]]))
      stop("LOCA2 variable ", vars[[k]], " should be in ", fut_units[[vars[[k]]]], " but the file says '", paste(u, collapse = "', '"), "': ", path)
    tm <- time(r)
    if (is.null(tm) || anyNA(tm)) stop("LOCA2 variable ", vars[[k]], " has no usable time axis: ", path)
    ym <- data.frame(year = as.integer(format(tm, "%Y")), month = as.integer(format(tm, "%m")))
    if (nrow(ym) %% 12L != 0L || any(ym$month != rep(1:12, nrow(ym) / 12L)) ||
        any(diff(ym$year[ym$month == 1L]) != 1L) || any(ym$year != rep(ym$year[ym$month == 1L], each = 12L)))
      stop("LOCA2 variable ", vars[[k]], " is not a run of complete, consecutive calendar years of monthly values: ", path)
    if (is.null(h$ym)) h$ym <- ym else if (!identical(h$ym, ym)) stop("LOCA2 variables cover different months: ", path)
    h[[k]] <- wb_crop_climate(fut_geo(r), site)
  }
  for (k in names(vars)[-1]) if (!compareGeom(h$ppt, h[[k]], stopOnError = FALSE)) stop("LOCA2 variables are on different grids: ", path)
  h$years <- unique(h$ym$year)
  h
}

# For every cell of a LOCA2 window, the cell number of a nearby cell that has data (itself where it has data). LOCA2 has
# no values over the ocean, and no humidity outside the United States, so cells on the coast and on the Canadian and
# Mexican borders would otherwise be interpolated from missing neighbours. Built once per file by spreading cell numbers
# outward one ring at a time.
fut_fill_index <- function(layer, label) {
  v <- values(layer, mat = FALSE)
  cells <- seq_len(ncell(layer)); cells[!is.finite(v)] <- NA
  n_missing <- sum(is.na(cells))
  if (n_missing == length(cells)) stop("The LOCA2 file has no ", label, " values over the study area.")
  idx <- rast(layer); values(idx) <- cells
  rings <- 0L
  while (anyNA(values(idx, mat = FALSE))) {
    idx <- focal(idx, w = 3, fun = "min", na.policy = "only", na.rm = TRUE)
    rings <- rings + 1L
    if (rings > sum(dim(layer)[1:2])) stop("Could not fill LOCA2 cells without ", label, " values.")
  }
  out <- as.integer(values(idx, mat = FALSE))
  attr(out, "missing") <- n_missing
  out
}

# Apply a fill index to a (multi-layer) window; the result is held in memory.
fut_fill <- function(x, idx) {
  v <- values(x)[idx, , drop = FALSE]
  if (anyNA(v)) stop("LOCA2 values are missing in cells that had data in the first month: ", sources(x)[1])
  out <- rast(x); values(out) <- v; names(out) <- names(x)
  out
}

# One model's file, ready to use: cropped variables, fill indices, and the model's own monthly means over the anchor
# years (1991-2020, the PRISM normal period).
fut_handle <- function(model, site, cfg, fut) {
  path <- fut_find_nc(fut, model)
  h <- fut_open(path, model, site, cfg)
  need <- sort(unique(c(cfg$base_years, cfg$years_annual, fut$anchor_years)))
  if (!all(need %in% h$years))
    stop("LOCA2 file for ", model, " covers ", min(h$years), "-", max(h$years), " but the run needs ",
         paste(range(need), collapse = "-"), ": ", path)
  old <- wb_enter("Preparing LOCA2 monthly means", paste(range(fut$anchor_years), collapse = "-"), path); on.exit(.wb_log$context <- old, add = TRUE)
  both <- h$ppt[[1]] + h$tmean[[1]]                       # data only where temperature and precipitation both exist
  h$idx <- list(ppt = fut_fill_index(both, "temperature and precipitation"))
  h$idx$tmean <- h$idx$ppt
  if (!is.null(h$vpd)) h$idx$vpd <- fut_fill_index(h$vpd[[1]], "vapour-pressure deficit")
  sel <- which(h$ym$year %in% fut$anchor_years)
  h$clim <- lapply(setNames(names(h$idx), names(h$idx)), function(k) {
    m <- tapp(h[[k]][[sel]], h$ym$month[sel], fun = "mean")
    if (nlyr(m) != 12L) stop("Monthly means for ", k, " did not produce 12 layers.")
    m <- fut_fill(m, h$idx[[k]]); names(m) <- paste0(k, "_", sprintf("%02d", 1:12)); m
  })
  msg <- sprintf("LOCA2 %s: %s; window %d x %d cells; filled from a neighbour: %d cells without temperature/precipitation%s.",
                 model, basename(path), nrow(h$ppt), ncol(h$ppt), attr(h$idx$ppt, "missing"),
                 if (is.null(h$vpd)) "" else sprintf(", %d without humidity", attr(h$idx$vpd, "missing")))
  cat(msg, "\n"); wb_note(msg)
  h
}

# Twelve months of one year on the LOCA2 window, gaps filled, with range checks.
fut_year <- function(h, year) {
  sel <- which(h$ym$year == year)
  if (length(sel) != 12L) stop("LOCA2 file has no complete year ", year, ": ", h$path)
  y <- lapply(setNames(names(h$idx), names(h$idx)), function(k) {
    x <- fut_fill(h[[k]][[sel]], h$idx[[k]]); names(x) <- paste0(k, "_", sprintf("%02d", 1:12)); x
  })
  rng <- function(x) range(as.numeric(as.matrix(global(x, c("min", "max")))))
  r <- rng(y$tmean); if (r[1] < -60 || r[2] > 60) stop("LOCA2 temperature for ", year, " is outside -60 to 60 C (", paste(round(r, 1), collapse = " to "), "): ", h$path)
  r <- rng(y$ppt); if (r[1] < 0 || r[2] > 5000) stop("LOCA2 precipitation for ", year, " is outside 0 to 5000 mm per month (", paste(round(r, 1), collapse = " to "), "): ", h$path)
  if (!is.null(y$vpd)) { r <- rng(y$vpd); if (r[1] < -0.001 || r[2] > 15) stop("LOCA2 vapour-pressure deficit for ", year, " is outside 0 to 15 kPa (", paste(round(r, 2), collapse = " to "), "): ", h$path) }
  y
}

# ---------------- The climate for one simulated year ----------------

# Anchored: PRISM normal + the model's departure from its own 1991-2020 monthly mean.
#   temperature    normal + (month - model mean)
#   precipitation  normal x (month + floor) / (model mean + floor)
#   VPD            normal x (month + floor) / (model mean + floor)
# The small floors (FUT_PR_FLOOR_MM, FUT_VPD_FLOOR_KPA) keep the ratio finite where the model's mean for a month is close
# to zero (desert summers); over the anchor years the ratio still averages exactly 1 in every cell and month.
fut_climate_year <- function(site, cfg, year) {
  f <- site$fut
  old <- wb_enter(paste("Reading LOCA2 climate,", f$mode), year, f$h$path); on.exit(.wb_log$context <- old, add = TRUE)
  y <- fut_year(f$h, year)
  if (f$mode == "direct") return(list(ppt = y$ppt, tmean = y$tmean, vpd = y$vpd))
  n <- site$normals; cl <- f$h$clim
  d <- c(y$tmean - cl$tmean, (y$ppt + f$pr_floor) / (cl$ppt + f$pr_floor))
  if (!is.null(y$vpd)) d <- c(d, (y$vpd + f$vpd_floor) / (cl$vpd + f$vpd_floor))
  d <- project(d, n$ppt, method = "bilinear")
  # The PRISM window is a rectangle around the study area and can reach past the edge of the LOCA2 file (corners of the
  # window in other states, Canada or Mexico). Those PRISM cells keep their normal; no modelled cell uses them, which
  # fut_check_extent() has confirmed for this file.
  tmean <- n$tmean + subst(d[[1:12]], NA, 0); names(tmean) <- paste0("tmean_", sprintf("%02d", 1:12))
  ppt <- n$ppt * subst(d[[13:24]], NA, 1); names(ppt) <- paste0("ppt_", sprintf("%02d", 1:12))
  vpd <- NULL
  if (!is.null(y$vpd)) { vpd <- n$vpd * subst(d[[25:36]], NA, 1); names(vpd) <- paste0("vpd_", sprintf("%02d", 1:12)) }
  list(ppt = ppt, tmean = tmean, vpd = vpd)
}

# Replacement for the historical reader. PRISM normals (year NULL) and any site without LOCA2 information go to the
# original function; a simulated year of a future run reads LOCA2.
wb_climate_year_prism <- wb_climate_year
wb_climate_year <- function(index, site, cfg, year = NULL) {
  if (is.null(year) || is.null(site$fut)) return(wb_climate_year_prism(index, site, cfg, year))
  fut_climate_year(site, cfg, year)
}

# ---------------- Direct runs: the site as seen from the LOCA2 grid ----------------

# Elevation difference between each DEM cell and the 6 km LOCA2 surface it is interpolated from (the historical run does
# the same against the 800 m PRISM grid), and the model's own 1991-2020 monthly means for spin-up.
fut_direct_site <- function(site, cfg, h) {
  old <- wb_enter("Locating the DEM on the LOCA2 grid", "not applicable", h$path); on.exit(.wb_log$context <- old, add = TRUE)
  low <- project(site$dem, h$clim$ppt[[1]], method = "bilinear")
  dz <- site$dem - project(low, site$dem, method = "bilinear")
  site$dz <- wb_repair_lapse_edges(site, cfg, low, dz)
  site
}

# Every modelled cell must be reached by the (gap-filled) LOCA2 window, in either mode.
fut_check_extent <- function(h, site, mem) {
  old <- wb_enter("Checking the LOCA2 file covers the study area", "not applicable", h$path); on.exit(.wb_log$context <- old, add = TRUE)
  reach <- values(project(h$clim$ppt[[1]], site$dem, method = "bilinear"), mat = FALSE)[mem$cells]
  lost <- sum(!is.finite(reach))
  if (lost > 0) stop(format(lost, big.mark = ","), " of ", format(mem$n, big.mark = ","),
                     " modelled cells lie outside the LOCA2 file; the study area extends beyond it: ", h$path)
  invisible(TRUE)
}

# How many modelled cells a climate grid reaches; cells it does not reach are not simulated.
fut_report_reach <- function(grid, mem, label) {
  n <- sum(grid$covered == 1L & is.finite(grid$dz))
  msg <- sprintf("%s grid reaches %s of %s modelled cells.", label, format(n, big.mark = ","), format(mem$n, big.mark = ","))
  cat(msg, "\n"); wb_note(msg)
  if (n < mem$n) warning(mem$n - n, " modelled cells are not reached by the ", label, " grid and will have no result.")
  invisible(n)
}

# ---------------- Configuration ----------------

fut_validate <- function(cfg, fut) {
  fut$climate_dir <- path.expand(fut$climate_dir)
  if (!dir.exists(fut$climate_dir)) stop("Missing CLIMATE_DIR: ", fut$climate_dir)
  if (!is.character(fut$models) || !length(fut$models) || anyDuplicated(fut$models)) stop("MODELS must be a vector of distinct model names.")
  if (!is.character(fut$modes) || !length(fut$modes) || any(!fut$modes %in% fut_modes) || anyDuplicated(fut$modes))
    stop('MODES must be "anchored", "direct" or both.')
  if (!is.character(fut$scenario) || length(fut$scenario) != 1L || !nzchar(fut$scenario)) stop("SCENARIO must be one name, e.g. 'ssp245'.")
  for (m in fut$models) fut_find_nc(fut, m)
  # Departures and direct spin-up use the model's means over the PRISM normal period.
  fut$anchor_years <- (cfg$prism_normal_end_year - 29L):cfg$prism_normal_end_year
  all_years <- c(cfg$base_years, cfg$years_annual)
  if (!is.list(fut$future_periods)) fut$future_periods <- list(fut$future_periods)
  fut$future_periods <- lapply(fut$future_periods, function(p) {
    if (!is.numeric(p) || length(p) != 2L || any(p != as.integer(p)) || p[1] > p[2] || !all(p[1]:p[2] %in% all_years))
      stop("FUTURE_PERIODS must be a list of c(first_year, last_year) pairs inside the simulated years ", min(all_years), "-", max(all_years), ".")
    as.integer(p)
  })
  for (field in c("pr_floor_mm", "vpd_floor_kpa"))
    if (!is.numeric(fut[[field]]) || length(fut[[field]]) != 1L || !is.finite(fut[[field]]) || fut[[field]] < 0) stop("Invalid ", toupper(field))
  if (is.null(fut$annual_datatype)) fut$annual_datatype <- "INT2S"
  if (!fut$annual_datatype %in% c("FLT4S", "INT2S")) stop("ANNUAL_DATATYPE must be 'INT2S' or 'FLT4S'.")
  if (cfg$analysis_initialization != "baseline_end") stop("Future runs carry soil and snow from the reference years into the future years: ANALYSIS_INITIALIZATION must be 'baseline_end'.")
  if (!identical(as.integer(cfg$base_years), fut$anchor_years))
    wb_note(sprintf("REFERENCE_YEARS (%d-%d) differ from the PRISM normal period (%d-%d) used to anchor the climate: change maps are relative to the reference years, the anchored climate equals PRISM over the normal period.",
                    min(cfg$base_years), max(cfg$base_years), min(fut$anchor_years), max(fut$anchor_years)))
  fut
}

# ---------------- Output ----------------

fut_label <- function(fut, model, mode) paste(model, fut$scenario, mode, sep = "_")
fut_period_tag <- function(p) paste0(p[1], "_", p[2])

fut_final_paths <- function(cfg, fut, model, mode) {
  lab <- fut_label(fut, model, mode)
  yrs <- range(c(cfg$base_years, cfg$years_annual))
  lapply(setNames(cfg$outputs_to_write, cfg$outputs_to_write), function(o) c(
    annual = file.path(cfg$output_root, sprintf("%s_%s_mm_%s_%d_%d.tif", toupper(o), cfg$month_tag, lab, yrs[1], yrs[2])),
    periods = file.path(cfg$output_root, sprintf("%s_%s_mm_%s_period_means.tif", toupper(o), cfg$month_tag, lab))))
}

# Per run: the annual stack over every simulated year (reference and future), and one file of period means and changes:
# the reference-period mean, then for each FUTURE_PERIODS entry the period mean and its change from the reference mean.
fut_summaries <- function(cfg, fut, model, mode, baseline, analysis) {
  old <- wb_enter("Writing final summaries"); on.exit(.wb_log$context <- old, add = TRUE)
  paths <- c(baseline$paths, analysis$paths); years <- c(cfg$base_years, cfg$years_annual)
  final <- fut_final_paths(cfg, fut, model, mode)
  for (o in cfg$outputs_to_write) {
    a <- rast(lapply(paths, function(f) wb_step("Reading annual result for summary", f, rast(f)[[o]])))
    names(a) <- paste0(toupper(o), "_", years)
    ref <- terra::app(a[[which(years %in% cfg$base_years)]], fun = "mean", na.rm = FALSE)
    out <- ref; nm <- paste0(o, "_mean_", fut_period_tag(range(cfg$base_years)))
    for (p in fut$future_periods) {
      m <- terra::app(a[[which(years >= p[1] & years <= p[2])]], fun = "mean", na.rm = FALSE)
      out <- c(out, m, m - ref)
      nm <- c(nm, paste0(o, c("_mean_", "_change_"), fut_period_tag(p)))
    }
    names(out) <- nm
    f_periods <- final[[o]][["periods"]]; f_annual <- final[[o]][["annual"]]
    wb_atomic_raster(out, f_periods, cfg)
    wb_atomic_raster(if (fut$annual_datatype == "INT2S") round(a) else a, f_annual, cfg, fut$annual_datatype)
    wb_note(c(paste("Annual result:", f_annual), paste("Period means and changes:", f_periods)))
    cat("Final stack:", f_annual, "\n")
  }
  invisible(final)
}

# Across models, per mode: mean of every period layer, plus the lowest and highest model for each change layer.
fut_ensemble <- function(cfg, fut) {
  old <- wb_enter("Writing ensemble summaries"); on.exit(.wb_log$context <- old, add = TRUE)
  for (mode in fut$modes) for (o in cfg$outputs_to_write) {
    files <- vapply(fut$models, function(m) {
      cr <- cfg; cr$output_root <- file.path(cfg$output_root, m, mode); cr$outputs_to_write <- o
      fut_final_paths(cr, fut, m, mode)[[o]][["periods"]] }, character(1))
    files <- files[file.exists(files)]
    if (length(files) < 2L) next
    x <- lapply(files, rast)
    layers <- names(x[[1]])
    if (!all(vapply(x, function(r) identical(names(r), layers), logical(1)))) { warning("Period layers differ between models; no ensemble summary for ", o, " (", mode, ")."); next }
    out <- NULL; nm <- character()
    for (l in layers) {
      s <- rast(lapply(x, function(r) r[[l]]))
      avg <- terra::app(s, fun = "mean", na.rm = FALSE)
      out <- if (is.null(out)) avg else c(out, avg)
      nm <- c(nm, paste0(l, "_model_mean"))
      if (grepl("_change_", l)) {
        lo <- terra::app(s, fun = "min", na.rm = FALSE)
        hi <- terra::app(s, fun = "max", na.rm = FALSE)
        out <- c(out, lo, hi)
        nm <- c(nm, paste0(l, c("_model_min", "_model_max")))
      }
    }
    names(out) <- nm
    path <- file.path(cfg$output_root, "ensemble", mode, sprintf("%s_%s_mm_ensemble_%s_%s_period_means.tif", toupper(o), cfg$month_tag, fut$scenario, mode))
    wb_atomic_raster(out, path, cfg)
    writeLines(c(paste("Models in this summary:", paste(names(files), collapse = ", ")), paste("Layers:", paste(nm, collapse = ", "))),
               file.path(dirname(path), sprintf("%s_ensemble_models.txt", toupper(o))))
    cat("Ensemble summary:", path, "\n"); wb_note(paste("Ensemble summary:", path))
  }
  invisible(NULL)
}

fut_readme <- function(cfg, fut) {
  ref <- fut_period_tag(range(cfg$base_years)); yrs <- range(c(cfg$base_years, cfg$years_annual))
  writeLines(c(
    "Future water balance, v6.2 engine (water_balance_future_v6.2.R).",
    sprintf("Climate: LOCA2 %s, monthly, from %s", fut$scenario, fut$climate_dir),
    sprintf("Simulated years %d-%d; reference period %s; months summed: %s.", yrs[1], yrs[2], gsub("_", "-", ref), cfg$month_tag),
    "",
    "Folders: <MODEL>/anchored and <MODEL>/direct, plus ensemble/<mode>.",
    "  anchored  PRISM 1991-2020 normal + each model's departure from its own 1991-2020 monthly mean (800 m climate grid).",
    "            Use for change from the PRISM starting point; comparable with the historical v6.2 maps.",
    "  direct    LOCA2 values fed in as they are (6 km climate grid, PRISM-fitted lapse rates within each cell).",
    "            Use for the climate models' own absolute values and their change from their own 1991-2020.",
    "",
    "Files in each run folder (CWD and AET, mm):",
    sprintf("  <VAR>_%s_mm_<MODEL>_%s_<mode>_%d_%d.tif   one layer per year (%s)", cfg$month_tag, fut$scenario, yrs[1], yrs[2],
            if (fut$annual_datatype == "INT2S") "whole millimetres, Int16" else "Float32"),
    sprintf("  <VAR>_%s_mm_<MODEL>_%s_<mode>_period_means.tif   layers: mean_%s, then for each future period mean_<period> and change_<period> (period mean minus the %s mean of the same run)",
            cfg$month_tag, fut$scenario, ref, gsub("_", "-", ref)),
    paste0("  Future periods: ", paste(vapply(fut$future_periods, function(p) paste(p, collapse = "-"), character(1)), collapse = ", ")),
    "  run_report.txt, SNOW_spinup_baseline_*.tif (direct runs), monthly/ and points/ when requested.",
    "ensemble/<mode>: across-model mean of every period layer; lowest and highest model for each change layer.",
    "",
    "Years 1991-2025 are each climate model's own weather, not the weather that occurred: compare period means with the",
    "historical run, never individual years.",
    "Static in every run: PRISM cloud-transmittance and dew-point normals, the wind climatology, soils and terrain."),
    file.path(cfg$output_root, "README.txt"))
}

# ---------------- One model, one mode ----------------

fut_run <- function(cfg, fut, index, site, mem, grids, points, h, mode, spin) {
  model <- h$model
  cr <- cfg; cr$output_root <- file.path(cfg$output_root, model, mode)
  cr$run_id <- wb_hash(list(base = cfg$run_id, model = model, mode = mode, climate = wb_file_stamp(h$path)))
  final <- fut_final_paths(cr, fut, model, mode)
  marker <- file.path(cr$output_root, "run_complete.rds")
  done <- tryCatch(cfg$resume && file.exists(marker) && identical(readRDS(marker)$run, cr$run_id) &&
                     (!cfg$write_rasters || all(file.exists(unlist(final)))), error = function(e) FALSE)
  if (isTRUE(done)) { cat("\n", model, " (", mode, "): already complete; skipped.\n", sep = ""); return(invisible(final)) }
  dir.create(cr$output_root, recursive = TRUE, showWarnings = FALSE)
  top_report <- .wb_log$report
  .wb_log$report <- file.path(cr$output_root, "run_report.txt"); on.exit(.wb_log$report <- top_report, add = TRUE)
  wb_note(c("", paste("===== Run started:", Sys.time(), "====="), paste("Model:", model), paste("Scenario:", fut$scenario),
            paste("Mode:", mode), paste("Climate file:", h$path), paste("Run fingerprint:", cr$run_id),
            paste("Shared setup (terrain, soil, lapse rates, PRISM normals) is logged in:", top_report)))
  cat("\n===== ", model, " (", mode, ") =====\n", sep = "")
  s <- site
  s$out_dir <- file.path(cfg$run_dir, model, mode)
  dir.create(s$out_dir, recursive = TRUE, showWarnings = FALSE)
  s$fut <- list(mode = mode, h = h, pr_floor = fut$pr_floor_mm, vpd_floor = fut$vpd_floor_kpa)
  if (mode == "direct") {
    # The model's own 1991-2020 monthly means stand in for the normals: spin-up starts from them.
    s$normals <- list(ppt = h$clim$ppt, tmean = h$clim$tmean, vpd = h$clim$vpd)
    fut_report_reach(wb_grid_for(s$normals$ppt, s$normals$tmean, s, mem, grids), mem, "LOCA2")
    spin <- wb_spinup(s, cr, mem, grids)
    wb_cleanup_temp(s)
  }
  baseline <- wb_run_sequence(s, cr, index, cr$base_years, spin, "baseline", mem, grids, points)
  analysis <- wb_run_sequence(s, cr, index, cr$years_annual, baseline$final, "analysis", mem, grids, points)
  if (cr$write_rasters) { fut_summaries(cr, fut, model, mode, baseline, analysis); wb_monthly_stacks(cr) }
  if (!is.null(points)) wb_points_write(cr, s, points, file.path(s$out_dir, c("baseline_years", "analysis_years")))
  for (f in c("spinup_status.txt", "spinup_diagnostics.csv")) {
    path <- file.path(s$out_dir, f)
    if (file.exists(path)) wb_note(c(paste("Diagnostics:", f), readLines(path)))
  }
  saveRDS(list(run = cr$run_id, finished = Sys.time()), marker)
  wb_note(paste("Completed:", Sys.time()))
  if (!cfg$diagnostics) unlink(file.path(s$out_dir, c("baseline_years", "analysis_years", "state")), recursive = TRUE)
  invisible(final)
}

# ---------------- Main ----------------

fut_main_impl <- function(cfg, fut) {
  wb_load_packages()
  cfg <- wb_validate(cfg)
  fut <- fut_validate(cfg, fut)
  wb_setup_runtime(cfg)
  jobs <- wb_jobs(cfg)
  if (length(jobs) != 1L) stop("Future runs model one combined area.")
  if (cfg$stop_after_setup) {
    cfg$run_dir <- file.path(cfg$work_root, "setup_check")
    site <- wb_prepare_site(cfg, jobs[[1]])
    for (m in fut$models) fut_open(fut_find_nc(fut, m), m, site, cfg)
    cat("Setup check finished: terrain, soil and the LOCA2 files open; no simulation was run.\n")
    return(invisible(cfg))
  }
  # PRISM is needed for its 1991-2020 normals only (lapse rates, anchoring, static cloud and dew-point fields).
  cfg_normals <- cfg; cfg_normals$base_years <- integer(); cfg_normals$years_annual <- integer()
  index <- wb_ensure_prism(cfg_normals)
  # Settings that change the science are part of the run fingerprint; the list of models and modes is not, so adding a
  # model later does not invalidate finished runs.
  cfg$future <- fut[c("scenario", "anchor_years", "pr_floor_mm", "vpd_floor_kpa", "code_md5")]
  cfg <- wb_start_run(cfg, index)
  fut_readme(cfg, fut)
  site <- wb_prepare_site(cfg, jobs[[1]])
  wb_cleanup_temp(site)
  site <- wb_add_lapse(site, cfg, index)
  wb_cleanup_temp(site)
  points <- if (!is.null(cfg$points_path)) wb_read_points(cfg, site) else NULL
  mem <- wb_memory_site(site, cfg, keep_cells = if (!cfg$write_rasters) unique(as.vector(points$nb_cells)) else NULL)
  mem$static <- wb_memory_static(mem, site, cfg)
  if (!is.null(points)) {
    points <- wb_points_rows(points, mem, cfg)
    cat("Points with modelled cells:", nrow(points$pts), if (!cfg$write_rasters) "(points-only run: no rasters will be written)" else "", "\n")
  }
  write.csv(mem$rad_summary, file.path(site$out_dir, "radiation_summary.csv"), row.names = FALSE)
  cat("PET method:", cfg$pet_method, "\n")
  grids <- list(anchored = new.env(), direct = new.env())
  spin_anchored <- NULL; site_direct <- NULL
  if ("anchored" %in% fut$modes) {
    # Anchored runs of every model share one spin-up: their 1991-2020 mean climate is the PRISM normal.
    fut_report_reach(wb_grid_for(site$normals$ppt, site$normals$tmean, site, mem, grids$anchored), mem, "PRISM")
    spin_anchored <- wb_spinup(site, cfg, mem, grids$anchored)
    wb_cleanup_temp(site)
    if (cfg$pet_diagnostics) { wb_pet_check(site, cfg, mem, grids$anchored, spin_anchored); wb_cleanup_temp(site) }
  }
  for (model in fut$models) {
    h <- fut_handle(model, site, cfg, fut)
    fut_check_extent(h, site, mem)
    for (mode in fut$modes) {
      if (mode == "direct" && is.null(site_direct)) site_direct <- fut_direct_site(site, cfg, h)
      fut_run(cfg, fut, index, if (mode == "direct") site_direct else site, mem, grids[[mode]], points, h, mode,
              if (mode == "anchored") spin_anchored else NULL)
      gc(verbose = FALSE)
    }
  }
  if (cfg$write_rasters) fut_ensemble(cfg, fut)
  for (f in c("site_metadata.txt", "lapse_coefficients.csv", "radiation_summary.csv", "spinup_status.txt", "spinup_diagnostics.csv")) {
    path <- file.path(site$out_dir, f)
    if (file.exists(path)) wb_note(c(paste("Diagnostics:", f), readLines(path)))
  }
  if (!cfg$diagnostics) {
    unlink(file.path(site$out_dir, "state"), recursive = TRUE)
    unlink(file.path(site$out_dir, "spinup_final_change_mm.tif"))
  }
  wb_note(paste("Completed:", Sys.time()))
  cat("\nCompleted. Outputs:", cfg$output_root, "\n")
  invisible(cfg)
}

# Same logging and warning handling as wb_main().
fut_main <- function(cfg, fut) {
  dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
  .wb_log$report <- file.path(cfg$output_root, "run_report.txt")
  .wb_log$context <- list(step = "startup", year = "not applicable", files = "not applicable")
  .wb_log$temp_dir <- NULL
  on.exit({
    if (!is.null(.wb_log$temp_dir)) unlink(.wb_log$temp_dir, recursive = TRUE)
    .wb_log$temp_dir <- NULL
    .wb_log$report <- NULL
  }, add = TRUE)
  wb_note(c("", paste("===== Future run started:", Sys.time(), "====="),
    "Shared setup is logged here; each model and mode has its own run_report.txt in its folder.",
    "Warnings are logged immediately. Raster read/write failures abort processing."))
  withCallingHandlers(fut_main_impl(cfg, fut),
    warning = function(w) {
      wb_condition(w, "WARNING")
      if (grepl("GDAL error|cannot read|cannot write|read failed|write failed|IReadBlock|TIFFRead|TIFFWrite|No space left|disk full|Operation timed out|Problem closing connection", conditionMessage(w), ignore.case = TRUE))
        stop(structure(list(message = paste("Raster I/O failure:", conditionMessage(w)), call = NULL),
                       class = c("wb_io_error", "error", "condition")))
      invokeRestart("muffleWarning")
    },
    error = function(e) wb_condition(e, "ERROR"))
}
