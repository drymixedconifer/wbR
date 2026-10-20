# Generated from inst/config/defaults.R. Keep synchronized with that source.
# Installed-package configuration does not depend on system.file() lookups.
.wbr_config_defaults <- function(.wb_functions_path) {
  # Original v5.5 driver settings, retained verbatim. Edit values via wb_config() or a project script.
  # ---- Paths and study settings ----
  BASE_DIR <- "/YOUR PATH"
  DEM_PATH <- file.path(BASE_DIR, "western_us_dem_500m.tif")
  AOI_PATH <- file.path(BASE_DIR, "westernstates.shp")
  AOI_BBOX <- NULL                # v5.4: set to c(xmin, xmax, ymin, ymax) in decimal degrees to model a box instead of the
                                  # shapefile, e.g. c(-122.6, -121.4, 43.5, 45.5) for the Oregon Cascades; NULL uses AOI_PATH.
                                  # Useful for calibration blocks: run a box with different snow settings and
                                  # validate_snow_snotel.R picks up the stations inside it. The box is part of the run
                                  # fingerprint, so give each block its own OUTPUT_ROOT.
  AWC_PATH <- "/YOUR PATH/POLARIS_AWC_western_states.tif"
  PRISM_ROOT <- "/YOUR PATH"
  PRISM_RES <- "800m"
  PRISM_MODE <- "full"   # "full": existing archive; "subset": compact buffered AOI cache; "local": no downloads
  PRISM_BUFFER_KM <- 25 # Additional buffer (lapse window + 4 PRISM cells are also respected)
  PRISM_NORMAL_END_YEAR <- 2020L  # Existing *_2020MM_avg_30y cache = 1991-2020 normals.
  
  # Final rasters go directly here; filenames follow MONTHS_OUT and the year ranges.
  # Compiled code and restart files live in an automatically managed local cache.
  OUTPUT_ROOT <- file.path(BASE_DIR, "CWD_outputs")
  PREPARED_ENVIRONMENT <- TRUE     # Reuse static DEM, slope, aspect, heat load, latitude, and soil grids.
  ENVIRONMENT_CACHE <- NULL       # NULL: hidden cache beside the input DEM; set a shared directory explicitly.
  
  BASE_YEARS <- 1948:1949        # Warm-up after spin-up; the baseline-mean rasters will be a 2-year average.
  YEARS_ANNUAL <- 1950:2025
  MONTHS_OUT <- 1:12
  OUTPUTS_TO_WRITE <- c("cwd", "aet")
  SUMMARIZE <- FALSE # No period means unless requested
  WRITE_BASELINE <- FALSE # Baseline used for initialization only
  RUN_AOIS_SEPARATELY <- FALSE   # Required for this flat layout (one combined geography).
  AOI_FEATURES_TO_RUN <- "all"
  AOI_LABEL_FIELD <- "name"
  
  # ---- Optional in-simulation zonal time series ----
  # Categorical/factor raster, polygon, or point layer. NULL disables collection.
  # This summarizes model arrays directly and DOES NOT require monthly rasters.
  ZONAL_PATH <- NULL
  ZONAL_ID_FIELD <- NULL       # Polygon/point attribute; NULL uses feature order
  ZONAL_MONTHLY <- c("cwd", "aet", "pet", "ppt", "water", "soil", "pack",
                     "snowfall", "melt", "tmean", "rs", "rnet", "vpd", "runoff", "balance")
  ZONAL_ANNUAL <- TRUE        # Annual CSV (monthly flux sums; state/radiation means)
  ZONAL_PLOTS <- FALSE        # Reserved, no automatic plotting in this version
  
  # ---- Snow (new in 5.3, revised in 5.4 and 5.5) ----
  # Rain/snow split: precipitation is all snow when the monthly mean temperature is at or below SNOW_T_ALL,
  # all rain at or above RAIN_T_ALL, and a straight line between. The transition describes months, not days:
  # a mountain month averaging +5 C still has snow days. Measured directly at 800 SNOTEL stations (140,000
  # November-March station-months, 1980-2025), observed snow retention falls smoothly from 1.0 at -6 C to
  # zero near +3.5 C once the winter melt below is accounted for; the values here reproduce that curve.
  # History: 4.0-5.2 used 0 and 3 C; 5.3 used -5 and 4 (five Malheur stations; 20% too little snow West-wide);
  # 5.4 used -2 and 5 with degree-day melt (peak SWE within 2%, melt-out a quarter month early, cold months
  # too snowy and mild months over-melted, cancelling). 5.5 fits phase and melt together (see below).
  SNOW_T_ALL <- -6                   # C. 5.4 value: -2. 5.2 behaviour: 0.
  RAIN_T_ALL <- 10                   # C. 5.4 value: 5. 5.2 behaviour: 3.
  
  # Snowmelt. "energy": melt (mm/day) = MELT_MM_PER_DEG_DAY x T + MELT_MM_PER_MJ x net radiation at the snow
  # surface, where net radiation = (1 - SNOW_ALBEDO) x incoming shortwave on the cell's slope under the month's
  # cloudiness, less net longwave loss; negative totals mean no melt (restricted degree-day form, Kustas et al.
  # 1994; Hock 2003). Winter melt is small even in a month averaging +1 C, late-spring melt is fast, and aspect
  # matters through the radiation term. Coefficients fitted jointly with the thresholds to 230,000 SNOTEL
  # station-months and checked by simulating 1950-2025 at every station: April 1 SWE -7%, peak -7% (SNOTEL
  # pillows read 5-10% above gauge precipitation, so this is the floor), melt-out within 0.01 month.
  # "degree_day": the 5.3/5.4 formulation, MELT_MM_PER_DEG_DAY x max(T, 0), with MELT_RADIATION_WEIGHT of the
  # factor scaled by the cell's clear-sky radiation relative to level ground; use -2, 5 and 3 with it.
  MELT_METHOD <- "energy"            # "energy" (5.5 default, fitted) or "degree_day" (5.3/5.4; thresholds -2 and 5, factor 3)
  MELT_MM_PER_DEG_DAY <- 2.3         # mm per C per day. Energy method: the temperature (turbulent-exchange) term (fitted 2.34).
                                     # Degree-day method: the whole melt factor on level ground (3 validated with -2/5).
  MELT_MM_PER_MJ <- 2.3              # Energy method only: mm of melt per MJ m-2 of net radiation (fitted 2.26; 3.0 = complete conversion).
  SNOW_ALBEDO <- 0.77                # Energy method only: snow-surface albedo (fitted 0.77; fresh 0.8-0.9, old/melting 0.5-0.6).
  MELT_RADIATION_WEIGHT <- 0.6       # Degree-day method only: share of the melt factor scaled with slope radiation (0 = plain 5.2 degree-day).
  
  # ---- Points (new in 5.5) ----
  # Monthly and annual values at a set of points, written as CSV (<OUTPUT_ROOT>/points/points_monthly.csv and
  # points_annual.csv, with a README). Points come from a CSV with longitude/latitude columns (WGS84 decimal
  # degrees; other columns are carried through) or any vector file terra can read (shapefile, GeoPackage, KML).
  # Every monthly variable listed under MONTHLY_OUTPUT below is included for every point and year, whatever
  # MONTHLY_OUTPUT says, plus the annual CWD and AET sums over MONTHS_OUT.
  POINTS_PATH <- NULL                # NULL: no point output; or e.g. file.path(BASE_DIR, "snotel_sites.csv")
  POINTS_ID_FIELD <- NULL            # Column/field holding the point id; NULL: first non-coordinate column (or a running number)
  POINTS_NEIGHBOURHOOD <- 1L         # 1: the point's own cell; 3 or 5: mean of the 3x3 or 5x5 cells around it
  WRITE_RASTERS <- TRUE              # FALSE (needs POINTS_PATH): simulate only the cells the points need and write no
                                     # rasters at all (no annual CWD/AET, monthly, spin-up or PET-check maps). Minutes
                                     # instead of hours for the West; the point values are identical to a full run.
  
  # ---- Monthly output (new in 5.2; four variables added in 5.5) ----
  MONTHLY_OUTPUT <- c("cwd", "aet", "pet", "water", "soil", "pack", "ppt", "snowfall", "melt", "tmean", "rs", "rnet")
                                     # FALSE for none; TRUE = c("cwd","aet"); or any subset of the twelve:
                                     #   cwd, aet, pet  monthly deficit, actual and potential ET (mm)
                                     #   water          rain + snowmelt reaching the soil (mm)
                                     #   soil, pack     end-of-month soil water and snowpack SWE (mm)
                                     #   ppt, snowfall  downscaled monthly precipitation and its snow share (mm)
                                     #   melt           snowmelt (mm)
                                     #   tmean          downscaled mean temperature (C)                       } always Float32,
                                     #   rs             incoming shortwave on the slope, cloud-adjusted (MJ m-2 d-1) } whatever
                                     #   rnet           net radiation at a snow surface (MJ m-2 d-1)         } MONTHLY_DATATYPE
  MONTHLY_YEARS <- list(cwd = 1950:2025, aet = 1950:2025, pet = 1950:2025, water = 1950:2025, soil = 1950:2025, pack = 1950:2025, ppt = 1950:2025, snowfall = 1950:2025, melt = 1950:2025, tmean = 1950:2025, rnet = 1980:2025, rs = 2000)
                                     # Per variable: a vector of years, "analysis" (default) or "all" (baseline years too).
                                     # rs is built from static radiation and cloud normals, so it is identical every year:
                                     # one year is enough. rnet varies only through temperature; the SNOTEL era suffices.
                                     # Disk at 500 m West-wide: roughly 1 GB per year for this list (about 90 GB in all);
                                     # trim the list or the years if that is more than you want to keep.
  MONTHLY_CHUNK <- 4L                # Monthly variables held in memory at once. The year is recomputed per chunk (fast,
                                     # identical results), so memory stays at about 450 + 96 x MONTHLY_CHUNK bytes per cell
                                     # whatever the length of MONTHLY_OUTPUT (West at 500 m: ~15 GB with 4, ~31 GB with 12).
  MONTHLY_DATATYPE <- "INT2S"        # "INT2S": whole millimetres, 5-10x smaller than Float32 (snowpack is mostly zero); "FLT4S": full precision.
                                     # zero); "FLT4S": full precision. INT2S is right for anything larger than a watershed.
  MONTHLY_STACK <- FALSE             # TRUE also writes one 12 x N-layer file per variable (large; e.g. 25 years at 30 m
                                     # over the Malheur is about 10 GB per variable). Per-year files are always written.
  # Files: <OUTPUT_ROOT>/monthly/<VAR>_monthly_<unit>_<year>.tif (unit mm, C or MJm2d), 12 layers named VAR_year_month. Monthly values
  # ignore MONTHS_OUT, which only governs the annual sums. Each extra variable adds 96 bytes per cell of RAM.
  
  # ---- Terrain and soil ----
  #DEM_AGGREGATE_FACT <- 3L      #10 m input -> roughly 30 m.
  DEM_AGGREGATE_FACT <- 1L       #Use native resolution (500 m West DEM)
  DEM_Z_UNITS <- "auto"         # Elevation-range heuristic; "meters" / "feet" override.
  METRIC_CRS <- "EPSG:5070"
  DEM_SOURCE_CRS <- NULL        # Only set to the VERIFIED source CRS to repair bad metadata.
  # Valid source CRSs are transformed to METRIC_CRS, never silently relabeled.
  
  AWC_UNITS <- "mm"             # Storage-unit assumption.
  AWC_DEPTH_MM <- NA_real_      # Optional KNOWN profile depth for plausibility checks.
  # A local POLARIS builder sums 0-200 cm storage in mm. If confirmed as the producer of your exact input TIFF, you can set AWC_DEPTH_MM <- 2000.
  # This depth does NOT rescale a raster already in mm.
  # Other explicit input options: "cm", "m", or "fraction". Only "fraction" requires a supplied depth and multiplies volumetric AWC by that depth.
  # Metadata conflicts stop the run. Suspicious value ranges produce warnings; values are never automatically reinterpreted from a guessed unit.
  
  # Heat load is divided by the flat-surface value at each cell's own latitude, so level ground has a modifier of exactly 1 and only slope and aspect change PET. Daylength uses each DEM cell's latitude. Snow equations and the monthly non-leap calendar are unchanged.
  
  # ---- Potential evapotranspiration and radiation (new in 5.0) ----
  PET_METHOD <- "penman_monteith"   # "priestley_taylor" or "hamon" (Hamon x heat load, exactly as 4.0).
  PT_ALPHA <- 1.26                   # Priestley-Taylor coefficient (1.26 is the standard well-watered value).
  SURFACE_ALBEDO <- 0.23             # Shortwave albedo for net radiation and terrain-reflected light (FAO reference surface).
  RAD_DIFFUSE_FRACTION <- 0.2        # Share of clear-sky radiation arriving as diffuse light (scaled by sky view on slopes).
  RAD_STEP_MINUTES <- 10             # Time step for integrating the direct beam over the mid-month day.
  # Shortwave on each cell = clear-sky radiation on its slope x monthly cloud transmittance.
  # Longwave uses mean monthly temperature, the transmittance, and dew point (FAO-56 eq. 39).
  CLOUD_TRANSMITTANCE <- "prism"     # PRISM 'soltrans' 1991-2020 monthly normals (downloaded once); or a constant 0-1 (e.g. 0.7).
  HUMIDITY_SOURCE <- "prism"         # PRISM 'tdmean' dew-point normals for the longwave term; or a constant relative humidity 0-1 (e.g. 0.5).
  PET_DIAGNOSTICS <- TRUE            # Write PET_check_*.csv/.tif (Hamon vs Priestley-Taylor over the normals) and RAD_terrain_modifier_by_month.tif.
  # Both normals are static monthly fields applied to every simulated year; year-to-year cloud and
  # humidity variation is not represented. Horizon shading by surrounding terrain is not included.
  
  # ---- Elevation adjustment of climate ----
  # Temperature and log-precipitation slopes against elevation are fitted to the
  # PRISM normals within a moving window this wide, so each region gets its own
  # slopes. Inf reproduces the earlier single West-wide fit. Slopes are limited
  # to the ranges below; LAPSE_regional_coefficients.tif maps the fitted values.
  LAPSE_WINDOW_KM <- 100
  LAPSE_T_LIMITS <- c(-0.0098, 0)       # C per m (dry adiabatic to isothermal).
  LAPSE_P_LIMITS <- c(-0.0005, 0.0015)  # Change in log precipitation per m.
  
  # ---- Initialization ----
  SPINUP_MAX_ITER <- 15L
  SOIL_TOL_MM <- 1
  PACK_TOL_MM <- 1
  SPINUP_MAX_PCT <- 0.1         # At most 0.1% of modeled cells may exceed each tolerance.
  # Convergence requires mean and fraction criteria only; snow maxima are reported separately.
  # Maxima, fractions, invalid-cell counts, and final residual maps are saved.
  # 0 requires every compared cell to meet the tolerance; 100 gives mean-only logic.
  # One two-layer SNOW_spinup_baseline_<years>.tif preserves snowpack and signed change (mm SWE).
  SPINUP_ACTION <- "stop"       # "warn": continue with explicit failed status; "stop": stop.
  # If not converged, inspect spinup_diagnostics.csv and spinup_final_change_mm.tif.
  # More iterations may help slowly stabilizing soil/snow but cannot guarantee equilibrium where this model accumulates snow indefinitely.
  
  
  # ---- Execution ----
  RESUME <- TRUE
  ALLOW_DOWNLOADS <- TRUE       # Only missing PRISM months; existing cache files are retained.
  DIAGNOSTICS <- FALSE          # TRUE retains local checkpoints/maps after success; FALSE cleans them up.
  STOP_AFTER_SETUP <- FALSE     # TRUE checks DEM/AOI/AWC only, with no climate download/year loop.
  LAPSE_SEED <- 20260908L       # Reproducible fitting if a very large AOI requires subsampling.
  WB_CORES <- local({
    n <- parallel::detectCores(logical = FALSE)
    if (!is.finite(n)) 1L else as.integer(max(1, min(6, n - 1)))
  })
  TERRA_MEMFRAC <- 0.6
  MEMORY_LIMIT_GB <- NA           # NA: stop if the in-memory engine would exceed 80% of physical RAM (macOS); or a number of GB.
  # The engine holds about 450 bytes per modeled cell (e.g. 75 million 10-m cells = 32 GB). If the estimate
  # printed at startup is too large, raise DEM_AGGREGATE_FACT (3 -> about 30 m) or shrink the AOI.
  PROJ_DATA_DIR <- "/opt/homebrew/share/proj"
  GDAL_DATA_DIR <- "/opt/homebrew/share/gdal"
  PROJ_NETWORK <- FALSE        # Use installed datum grids; avoid failed on-demand downloads.
  GDAL_OPTS <- c("COMPRESS=LZW", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  # Each finished year writes its CWD/AET raster and one restart file holding the
  # carried soil and snow (full precision), so an interrupted run resumes there.
  # Final CWD/AET stacks and baseline means use Float32, as in v2.2.
  # Final baseline means require all baseline years.
  # Warnings and diagnostic summaries are saved immediately in run_report.txt.
  # Failed runs retain local restart files. After success, only ending states are
  # retained by default; rerunning a completed run then recomputes the annual series.
  # Changing years changes filenames. Repeating the same months/years replaces
  # those result rasters; the report appends a dated run record. No AOI names are used.
  
  # v6 atmospheric inputs (all other original package controls remain unchanged)
  WIND_PATH <- NULL
  WIND_CONSTANT <- 2.0
  AERO_SCALE <- 1.0
  
  cfg <- list(
    base_dir = BASE_DIR, dem_path = DEM_PATH, aoi_path = AOI_PATH, aoi_bbox = AOI_BBOX, awc_path = AWC_PATH,
    prism_root = PRISM_ROOT, prism_res = PRISM_RES, prism_mode = PRISM_MODE,
    prism_buffer_km = PRISM_BUFFER_KM, prism_normal_end_year = PRISM_NORMAL_END_YEAR,
    output_root = OUTPUT_ROOT, prepared_environment = PREPARED_ENVIRONMENT,
    environment_cache = ENVIRONMENT_CACHE, base_years = BASE_YEARS, years_annual = YEARS_ANNUAL,
    months_out = MONTHS_OUT, outputs_to_write = OUTPUTS_TO_WRITE, summarize = SUMMARIZE, write_baseline = WRITE_BASELINE,
    run_aois_separately = RUN_AOIS_SEPARATELY, aoi_features_to_run = AOI_FEATURES_TO_RUN,
    aoi_label_field = AOI_LABEL_FIELD, dem_aggregate_fact = DEM_AGGREGATE_FACT,
    dem_z_units = DEM_Z_UNITS, metric_crs = METRIC_CRS, dem_source_crs = DEM_SOURCE_CRS,
    awc_units = AWC_UNITS, awc_depth_mm = AWC_DEPTH_MM,
    lapse_window_km = LAPSE_WINDOW_KM, lapse_t_limits = LAPSE_T_LIMITS, lapse_p_limits = LAPSE_P_LIMITS,
    spinup_max_iter = SPINUP_MAX_ITER, soil_tol_mm = SOIL_TOL_MM, pack_tol_mm = PACK_TOL_MM,
    spinup_max_pct = SPINUP_MAX_PCT, spinup_action = SPINUP_ACTION,
    pet_method = PET_METHOD, wind_path = WIND_PATH, wind_constant = WIND_CONSTANT, aero_scale = AERO_SCALE, pt_alpha = PT_ALPHA, surface_albedo = SURFACE_ALBEDO,
    rad_diffuse_fraction = RAD_DIFFUSE_FRACTION, rad_step_minutes = RAD_STEP_MINUTES,
    cloud_transmittance = CLOUD_TRANSMITTANCE, humidity_source = HUMIDITY_SOURCE,
    pet_diagnostics = PET_DIAGNOSTICS, memory_limit_gb = MEMORY_LIMIT_GB,
    monthly_output = MONTHLY_OUTPUT, monthly_years = MONTHLY_YEARS, monthly_stack = MONTHLY_STACK,
    monthly_datatype = MONTHLY_DATATYPE, monthly_chunk = MONTHLY_CHUNK,
    snow_t_all = SNOW_T_ALL, rain_t_all = RAIN_T_ALL, melt_mm_per_deg_day = MELT_MM_PER_DEG_DAY,
    melt_radiation_weight = MELT_RADIATION_WEIGHT,
    melt_method = MELT_METHOD, melt_mm_per_mj = MELT_MM_PER_MJ, snow_albedo = SNOW_ALBEDO,
    points_path = POINTS_PATH, points_id_field = POINTS_ID_FIELD, points_neighbourhood = POINTS_NEIGHBOURHOOD,
    write_rasters = WRITE_RASTERS,
    zonal_path = ZONAL_PATH, zonal_id_field = ZONAL_ID_FIELD,
    zonal_monthly = ZONAL_MONTHLY, zonal_annual = ZONAL_ANNUAL, zonal_plots = ZONAL_PLOTS,
    resume = RESUME, allow_downloads = ALLOW_DOWNLOADS, diagnostics = DIAGNOSTICS,
    stop_after_setup = STOP_AFTER_SETUP, lapse_seed = LAPSE_SEED, wb_cores = WB_CORES,
    terra_memfrac = TERRA_MEMFRAC, proj_data_dir = PROJ_DATA_DIR, gdal_data_dir = GDAL_DATA_DIR, proj_network = PROJ_NETWORK,
    gdal_opts = GDAL_OPTS,
    # Scientific settings are fingerprinted separately; performance-only edits here do not invalidate checkpoints. The full functions/engine file is hashed.
    code_md5 = unname(tools::md5sum(.wb_functions_path))
  )
  
  cfg
}
