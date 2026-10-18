###########################################
#### Validate modeled snowpack against ####
## SNOTEL stations and NRCS snow courses ##
###########################################

# This is going to compare the water-balance model's end-of-month snow water equivalent (SWE) from the monthly PACK rasters (v5.2+) with measurements inside the model area to:
#   SNOTEL: automated daily SWE (and precipitation), mostly 1980 onward.
#   Snow courses: manual SWE measured near the first of each month, typically February-May or June, many from the 1930s

# A first-of-month course value is compared with the model's end-of-previous-month value.
# Station data come from the NRCS web service: SNOTEL through the snotelr package (install.packages("snotelr")), snow courses through the NRCS report generator (read directly). Both are cached under the output folder so later runs do not download again. Local CSV folders can be used instead (see
# SNOTEL_LOCAL_DIR / SNOWCOURSE_LOCAL_DIR; column lists below).

# Metrics per station and water year (Oct-Sep):
#   * April 1 SWE (station 1 April vs model end of March), the standard metric;
#     May 1, same. Both sources.
#   * SWE climatology: mean SWE by month, model vs station. Both sources.
#   * Peak SWE and melt-out month (SNOTEL only; courses are too sparse in time).
#   * Precipitation: gauge Nov-Mar and water-year totals vs the model's downscaled precipitation (SNOTEL only; needs "ppt" in MONTHLY_OUTPUT).
# Peak SWE mainly reflects the rain/snow split (given correct precipitation); melt-out mainly reflects the melt rate. A decade-by-decade table shows whether the snow settings hold across the record, e.g. 1950s vs 2010s.
# Outputs go to <OUTPUT_ROOT>/validation_snotel/.

suppressPackageStartupMessages(library(terra))
terra::projNetwork(FALSE)

# ---- Settings ---- #
BASE_DIR <- "/YOUR PATH"
OUTPUT_ROOT <- file.path(BASE_DIR, "CWD_outputs")     # water-balance output folder (monthly/PACK_monthly_mm_<year>.tif)
YEARS <- 1950:2025                                     # calendar years to use (those with PACK files are used)
USE_SNOTEL <- TRUE
USE_SNOW_COURSES <- TRUE
SNOTEL_LOCAL_DIR <- NULL        # NULL: download; or a folder of CSVs with site_id, site_name, latitude, longitude, elev, date, snow_water_equivalent[, precipitation] (mm)
SNOWCOURSE_LOCAL_DIR <- NULL    # NULL: download; or a folder of CSVs with site_id, site_name, latitude, longitude, elev, date (first of month), snow_water_equivalent (mm)
STATION_IDS <- NULL             # NULL: every station inside the model area; or a vector of site_id values (either source)
MIN_YEARS <- 5                  # stations with fewer usable water years are dropped
MELTOUT_MM <- 10                # SWE below this counts as melted out
NEIGHBOURHOOD <- 3              # model cells per side averaged around the station (3 = 3x3)
OUT_DIR <- file.path(OUTPUT_ROOT, "validation_snotel")

# ---- Model snowpack ---- #
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(OUT_DIR, "snotel_cache"); dir.create(cache, showWarnings = FALSE)
pack_file <- function(y) file.path(OUTPUT_ROOT, "monthly", paste0("PACK_monthly_mm_", y, ".tif"))
ppt_file <- function(y) file.path(OUTPUT_ROOT, "monthly", paste0("PPT_monthly_mm_", y, ".tif"))
years <- YEARS[file.exists(pack_file(YEARS))]
if (!length(years)) stop("No PACK_monthly_mm_<year>.tif files under ", file.path(OUTPUT_ROOT, "monthly"),
                         ". Run the model with \"pack\" in MONTHLY_OUTPUT first.")
cat("Model years with monthly snowpack:", paste(range(years), collapse = "-"), "(", length(years), "years )\n")
ppt_years <- years[file.exists(ppt_file(years))]
tmpl <- rast(pack_file(years[1]))[[1]]
model_crs <- crs(tmpl)
bbox <- ext(project(as.polygons(ext(tmpl), crs = model_crs), "EPSG:4326"))
inside_bbox <- function(lon, lat) lon >= xmin(bbox) & lon <= xmax(bbox) & lat >= ymin(bbox) & lat <= ymax(bbox)
std_cols <- c("site_id", "site_name", "latitude", "longitude", "elev")

# ---- NRCS station lists ---- #
# The NRCS "yearcount" page lists every station of a network with coordinates and elevation (feet). snotelr parses the same page for SNOTEL; this parser also keeps alphanumeric snow-course ids such as 22G21.
nrcs_station_list <- function(network) {
  f <- file.path(cache, paste0("stations_", network, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  for (p in c("httr", "rvest")) if (!requireNamespace(p, quietly = TRUE))
    stop("Install the ", p, " package (it comes with snotelr) or supply local CSV folders.")
  page <- httr::GET("https://wcc.sc.egov.usda.gov/nwcc/yearcount",
                    query = list(network = network, counttype = "listwithdiscontinued"))
  tabs <- rvest::html_table(rvest::read_html(page))
  if (!length(tabs)) stop("No station table returned for network ", network)
  df <- as.data.frame(tabs[[which.max(sapply(tabs, nrow))]])
  names(df) <- tolower(trimws(names(df)))
  pick <- function(...) { for (n in c(...)) if (n %in% names(df)) return(df[[n]]); rep(NA, nrow(df)) }
  raw_name <- as.character(pick("site_name", "site name", "name"))
  code <- regmatches(raw_name, regexpr("\\(([^)]*)\\)[^(]*$", raw_name))
  code <- gsub("[()]", "", code)
  out <- data.frame(
    site_id = ifelse(nzchar(code), code, NA),
    site_name = tolower(trimws(sub("\\([^)]*\\)[^(]*$", "", raw_name))),
    state = as.character(pick("state", "st")),
    latitude = suppressWarnings(as.numeric(pick("lat", "latitude"))),
    longitude = suppressWarnings(as.numeric(pick("lon", "longitude"))),
    elev = round(suppressWarnings(as.numeric(gsub(",", "", pick("elev", "elevation")))) * 0.3048),
    start = as.character(pick("start")), end = as.character(pick("enddate", "end")),
    stringsAsFactors = FALSE)
  out <- out[!is.na(out$site_id) & is.finite(out$latitude) & is.finite(out$longitude), ]
  saveRDS(out, f); out
}
year_of <- function(x) suppressWarnings(as.integer(regmatches(x, regexpr("[12][0-9]{3}", x))))

# ---- SNOTEL (daily) ---- #
read_local <- function(dir, need) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No CSV files in ", dir)
  do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE); names(d) <- tolower(names(d))
    miss <- setdiff(need, names(d))
    if (length(miss)) stop(basename(f), " lacks columns: ", paste(miss, collapse = ", "))
    d$date <- as.Date(d$date)
    if (!"precipitation" %in% names(d)) d$precipitation <- NA_real_
    d$site_id <- as.character(d$site_id)
    d[, c(std_cols, "date", "snow_water_equivalent", "precipitation")]
  }))
}
fetch_snotel <- function() {
  if (!requireNamespace("snotelr", quietly = TRUE))
    stop("Install the snotelr package (install.packages('snotelr')) or set SNOTEL_LOCAL_DIR / USE_SNOTEL <- FALSE.")
  meta <- nrcs_station_list("sntl")
  meta <- meta[inside_bbox(meta$longitude, meta$latitude), ]
  if (!is.null(STATION_IDS)) meta <- meta[meta$site_id %in% as.character(STATION_IDS), ]
  cat("SNOTEL stations in the model area:", nrow(meta), "\n")
  out <- list()
  for (i in seq_len(nrow(meta))) {
    if (i %% 50 == 0) cat("  ...", i, "of", nrow(meta), "SNOTEL stations\n")
    id <- meta$site_id[i]
    f <- file.path(cache, paste0("site_", id, "_v2.rds"))
    d <- if (file.exists(f)) readRDS(f) else tryCatch({
      x <- snotelr::snotel_download(site_id = as.numeric(id), internal = TRUE)
      names(x) <- tolower(names(x))
      x <- x[, intersect(c("date", "snow_water_equivalent", "precipitation", "precipitation_cumulative"), names(x))]
      saveRDS(x, f); x
    }, error = function(e) { message("  SNOTEL download failed for site ", id, ": ", conditionMessage(e)); NULL })
    if (is.null(d) || !nrow(d)) next
    d$date <- as.Date(d$date)
    if (!"precipitation" %in% names(d)) d$precipitation <- NA_real_
    if ("precipitation_cumulative" %in% names(d)) {
      d <- d[order(d$date), ]
      inc <- c(NA, diff(d$precipitation_cumulative)); inc[inc < 0 | inc > 300] <- NA
      d$precipitation <- ifelse(is.na(d$precipitation), inc, d$precipitation)
    }
    d <- d[!is.na(d$snow_water_equivalent), c("date", "snow_water_equivalent", "precipitation")]
    if (!nrow(d)) next
    out[[length(out) + 1]] <- data.frame(meta[i, std_cols], d, row.names = NULL)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ---- Snow courses (monthly, first of month) ---- #
# NRCS report generator: monthly start-of-period SWE for one course, period of record.
course_url <- function(id, state) sprintf(
  "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/monthly/start_of_period/%s:%s:SNOW%%7Cid=%%22%%22%%7Cname/POR_BEGIN,POR_END/WTEQ::value",
  utils::URLencode(id, reserved = TRUE), state)
parse_course_csv <- function(txt) {
  lines <- txt[!grepl("^#", txt) & nzchar(txt)]
  if (length(lines) < 2) return(NULL)
  d <- read.csv(text = lines, stringsAsFactors = FALSE, check.names = FALSE)
  names(d)[1] <- "date"
  swe_col <- grep("snow water equivalent", names(d), ignore.case = TRUE)
  if (!length(swe_col)) return(NULL)
  dt <- as.Date(paste0(d$date, "-01"), "%Y-%m-%d")                       # "1950-04"
  if (all(is.na(dt))) dt <- as.Date(paste("01", d$date), "%d %b %Y")       # "Apr 1950"
  if (all(is.na(dt))) dt <- as.Date(d$date)                                # "1950-04-01"
  data.frame(date = dt, snow_water_equivalent = suppressWarnings(as.numeric(d[[swe_col[1]]])) * 25.4,
             precipitation = NA_real_)[!is.na(dt), ]
}
fetch_courses <- function() {
  meta <- nrcs_station_list("snow")
  meta <- meta[inside_bbox(meta$longitude, meta$latitude), ]
  if (!is.null(STATION_IDS)) meta <- meta[meta$site_id %in% as.character(STATION_IDS), ]
  cat("Snow courses in the model area:", nrow(meta), "\n")
  out <- list()
  for (i in seq_len(nrow(meta))) {
    if (i %% 50 == 0) cat("  ...", i, "of", nrow(meta), "snow courses\n")
    id <- meta$site_id[i]
    f <- file.path(cache, paste0("course_", gsub("[^A-Za-z0-9]", "_", id), ".rds"))
    d <- if (file.exists(f)) readRDS(f) else tryCatch({
      u <- course_url(id, meta$state[i])
      txt <- tryCatch(readLines(url(u), warn = FALSE), error = function(e) { Sys.sleep(3); readLines(url(u), warn = FALSE) })  # one retry for transient NRCS refusals
      x <- parse_course_csv(txt)
      if (is.null(x)) stop("no SWE column in the NRCS response")
      saveRDS(x, f); x
    }, error = function(e) { message("  snow course download failed for ", id, ": ", conditionMessage(e)); NULL })
    if (is.null(d) || !nrow(d)) next
    d <- d[!is.na(d$snow_water_equivalent), ]
    if (!nrow(d)) next
    out[[length(out) + 1]] <- data.frame(meta[i, std_cols], d, row.names = NULL)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ---- Assemble observations as station-month rows ---- #
# Each row: one station, one calendar year-month, the SWE at the END of that month (SNOTEL: mean of the last three days; course: the first-of-next-month reading),  and the station's monthly precipitation (SNOTEL only).
obs_cache <- file.path(cache, "observations_monthly.rds")   # delete this file to rebuild from the station downloads
obs_rows <- list()
if (file.exists(obs_cache)) {
  obs_rows <- list(readRDS(obs_cache)); cat("Using cached station-month observations:", obs_cache, "\n")
} else {
if (USE_SNOTEL) {
  sn <- if (is.null(SNOTEL_LOCAL_DIR)) fetch_snotel() else read_local(SNOTEL_LOCAL_DIR, c(std_cols, "date", "snow_water_equivalent"))
  if (!is.null(sn) && nrow(sn)) {
    # Group ~10 million daily rows by station-month with rowsum() (fast, little memory).
    sn <- sn[!duplicated(sn[, c("site_id", "date")]), ]            # one row per station and day
    ym <- format(sn$date, "%Y-%m"); dom <- as.integer(format(sn$date, "%d"))
    key <- paste(sn$site_id, ym)
    late <- dom >= 26 & !is.na(sn$snow_water_equivalent)
    swe_sum <- rowsum(sn$snow_water_equivalent[late], key[late]); swe_n <- rowsum(rep(1, sum(late)), key[late])
    swe <- data.frame(key = rownames(swe_sum), swe_obs = swe_sum[, 1] / swe_n[, 1], stringsAsFactors = FALSE)
    hasp <- !is.na(sn$precipitation)
    if (any(hasp)) {
      p_sum <- rowsum(sn$precipitation[hasp], key[hasp]); p_n <- rowsum(rep(1, sum(hasp)), key[hasp])
      pm <- data.frame(key = rownames(p_sum), ppt_obs = ifelse(p_n[, 1] >= 25, p_sum[, 1], NA_real_), stringsAsFactors = FALSE)
      swe <- merge(swe, pm, by = "key", all.x = TRUE)
    } else swe$ppt_obs <- NA_real_
    swe$site_id <- sub(" [0-9]{4}-[0-9]{2}$", "", swe$key); swe$ym <- substr(swe$key, nchar(swe$key) - 6, nchar(swe$key))
    m <- merge(swe[, c("site_id", "ym", "swe_obs", "ppt_obs")], unique(sn[, std_cols]), by = "site_id")
    m$source <- "snotel"
    obs_rows[[length(obs_rows) + 1]] <- m
    rm(sn, ym, dom, key, late, swe); invisible(gc())
  }
}
if (USE_SNOW_COURSES) {
  sc <- if (is.null(SNOWCOURSE_LOCAL_DIR)) fetch_courses() else read_local(SNOWCOURSE_LOCAL_DIR, c(std_cols, "date", "snow_water_equivalent"))
  if (!is.null(sc) && nrow(sc)) {
    # a reading dated the 1st of month M describes the end of month M-1
    prev <- as.Date(format(sc$date, "%Y-%m-01")) - 1
    sc$ym <- format(prev, "%Y-%m")
    key <- paste(sc$site_id, sc$ym)
    s_sum <- rowsum(sc$snow_water_equivalent, key); s_n <- rowsum(rep(1, nrow(sc)), key)
    m <- data.frame(site_id = sub(" [0-9]{4}-[0-9]{2}$", "", rownames(s_sum)), ym = substr(rownames(s_sum), nchar(rownames(s_sum)) - 6, nchar(rownames(s_sum))),
                    swe_obs = s_sum[, 1] / s_n[, 1], ppt_obs = NA_real_, stringsAsFactors = FALSE)
    m <- merge(m, unique(sc[, std_cols]), by = "site_id")
    m$source <- "course"
    obs_rows[[length(obs_rows) + 1]] <- m
  }
}
if (!length(obs_rows)) stop("No station data available.")
saveRDS(do.call(rbind, obs_rows), obs_cache)
}
obs <- do.call(rbind, obs_rows)
obs$site_id <- paste0(ifelse(obs$source == "snotel", "SNTL_", "SNOW_"), obs$site_id)
obs$year <- as.integer(substr(obs$ym, 1, 4)); obs$month <- as.integer(substr(obs$ym, 6, 7))
obs <- obs[obs$year %in% years & inside_bbox(obs$longitude, obs$latitude), ]
stations <- unique(obs[, c(std_cols, "source")])
cat("Stations with data in the model years:", sum(stations$source == "snotel"), "SNOTEL,", sum(stations$source == "course"), "snow courses\n")

# ---- Model values at the stations ---- #
pts <- project(vect(stations, geom = c("longitude", "latitude"), crs = "EPSG:4326"), model_crs)
cells <- cellFromXY(tmpl, crds(pts))
in_model <- !is.na(cells) & !is.na(values(tmpl, mat = FALSE)[cells])
if (any(!in_model)) cat("Outside the modelled area (dropped):", sum(!in_model), "stations\n")
stations <- stations[in_model, ]; pts <- pts[in_model, ]
stopifnot(nrow(pts) == nrow(stations))
if (!nrow(stations)) stop("No station falls on a modelled cell.")
# Mean of the NEIGHBOURHOOD x NEIGHBOURHOOD cells around each station, read directly by cell number (no focal pass over the full raster). Each year's values are cached, so a crash or a rerun against the same model output only extracts the years not yet done.
half <- (NEIGHBOURHOOD - 1) %/% 2
rc <- rowColFromCell(tmpl, cells[in_model])
nb_cells <- do.call(cbind, lapply(-half:half, function(dr) sapply(-half:half, function(dc) {
  r <- rc[, 1] + dr; c <- rc[, 2] + dc
  ok <- r >= 1 & r <= nrow(tmpl) & c >= 1 & c <= ncol(tmpl)
  out <- rep(NA_integer_, nrow(rc)); out[ok] <- cellFromRowCol(tmpl, r[ok], c[ok]); out })))
nb_cells <- matrix(nb_cells, nrow = nrow(rc))
nb_ok <- !is.na(nb_cells)
at_pts <- function(r) {
  vals <- r[as.vector(ifelse(nb_ok, nb_cells, 1L))]  # rows = station x neighbour cell, columns = layers
  out <- as.data.frame(lapply(seq_len(nlyr(r)), function(k) { v <- vals[[k]]; v[!as.vector(nb_ok)] <- NA
    rowMeans(matrix(v, nrow = nrow(nb_cells)), na.rm = TRUE) }))
  names(out) <- names(r); out[is.nan(as.matrix(out))] <- NA
  stopifnot(nrow(out) == nrow(stations)); out }
model_cache <- file.path(cache, "model_at_stations"); dir.create(model_cache, showWarnings = FALSE)
model_stamp <- function(y) { f <- c(pack_file(y), if (y %in% ppt_years) ppt_file(y)); paste(basename(f), file.size(f), file.mtime(f), collapse = ";") }
station_key <- paste(nrow(stations), paste(range(stations$latitude), collapse = ","), paste(range(stations$longitude), collapse = ","))
model_rows <- list()
for (y in years) {
  f <- file.path(model_cache, paste0("model_", y, ".rds"))
  cached <- if (file.exists(f)) readRDS(f) else NULL
  if (!is.null(cached) && identical(cached$stamp, model_stamp(y)) && identical(cached$stations, station_key)) {
    model_rows[[length(model_rows) + 1]] <- cached$rows; next }
  if (which(years == y) %% 10 == 1) cat("  extracting model values at stations:", y, "\n")
  v <- at_pts(rast(pack_file(y)))
  vp <- if (y %in% ppt_years) at_pts(rast(ppt_file(y))) else NULL
  rows <- do.call(rbind, lapply(1:12, function(m) data.frame(
    site_id = stations$site_id, year = y, month = m, swe_model = v[[m]], ppt_model = if (is.null(vp)) NA_real_ else vp[[m]])))
  saveRDS(list(stamp = model_stamp(y), stations = station_key, rows = rows), f)
  model_rows[[length(model_rows) + 1]] <- rows
}
mod <- do.call(rbind, model_rows)
d <- merge(mod, obs[, c("site_id", "year", "month", "swe_obs", "ppt_obs", "source")], by = c("site_id", "year", "month"))
d <- merge(d, stations[, std_cols], by = "site_id")
d$wyear <- ifelse(d$month >= 10, d$year + 1L, d$year)
d$decade <- paste0(floor(d$wyear / 10) * 10, "s")
write.csv(d, file.path(OUT_DIR, "swe_station_month.csv"), row.names = FALSE)

# ---- Water-year metrics ---- #
wy_order <- function(m) ifelse(m >= 10, m - 9, m + 3)
metrics <- do.call(rbind, lapply(split(d, list(d$site_id, d$wyear), drop = TRUE), function(x) {
  x <- x[order(wy_order(x$month)), ]
  src <- x$source[1]
  get <- function(m, col) { v <- x[[col]][x$month == m]; if (length(v)) v[1] else NA_real_ }
  r <- data.frame(site_id = x$site_id[1], source = src, wyear = x$wyear[1], n_months = nrow(x),
                  apr1_obs = get(3, "swe_obs"), apr1_model = get(3, "swe_model"),
                  may1_obs = get(4, "swe_obs"), may1_model = get(4, "swe_model"),
                  peak_obs = NA_real_, peak_model = NA_real_, meltout_obs = NA_real_, meltout_model = NA_real_,
                  ppt_novmar_obs = NA_real_, ppt_novmar_model = NA_real_, ppt_wy_obs = NA_real_, ppt_wy_model = NA_real_)
  if (src == "snotel" && nrow(x) >= 10) {
    r$peak_obs <- max(x$swe_obs); r$peak_model <- max(x$swe_model)
    meltout <- function(v) { i <- which.max(v); j <- which(seq_along(v) > i & v < MELTOUT_MM); if (length(j)) wy_order(x$month[min(j)]) else NA }
    r$meltout_obs <- meltout(x$swe_obs); r$meltout_model <- meltout(x$swe_model)
    win <- x$month %in% c(11, 12, 1, 2, 3)
    if (sum(win) == 5 && all(is.finite(x$ppt_obs[win])) && all(is.finite(x$ppt_model[win]))) {
      r$ppt_novmar_obs <- sum(x$ppt_obs[win]); r$ppt_novmar_model <- sum(x$ppt_model[win]) }
    if (nrow(x) == 12 && all(is.finite(x$ppt_obs)) && all(is.finite(x$ppt_model))) {
      r$ppt_wy_obs <- sum(x$ppt_obs); r$ppt_wy_model <- sum(x$ppt_model) }
  }
  r
}))
metrics$meltout_diff_months <- metrics$meltout_model - metrics$meltout_obs
metrics$decade <- paste0(floor(metrics$wyear / 10) * 10, "s")
write.csv(metrics, file.path(OUT_DIR, "swe_station_wateryear.csv"), row.names = FALSE)

# ---- Station summary ---- #
safe_cor <- function(a, b) if (sum(complete.cases(a, b)) > 3) cor(a, b, use = "complete.obs") else NA
pct <- function(m, o) 100 * (mean(m, na.rm = TRUE) / mean(o, na.rm = TRUE) - 1)
summ <- do.call(rbind, lapply(split(metrics, metrics$site_id), function(x) {
  x <- x[is.finite(x$apr1_obs) & is.finite(x$apr1_model), ]
  if (nrow(x) < MIN_YEARS) return(NULL)
  dd <- d[d$site_id == x$site_id[1], ]
  data.frame(site_id = x$site_id[1], source = x$source[1], water_years = nrow(x),
    first_year = min(x$wyear), last_year = max(x$wyear),
    apr1_obs_mm = mean(x$apr1_obs), apr1_model_mm = mean(x$apr1_model), apr1_bias_pct = pct(x$apr1_model, x$apr1_obs),
    apr1_cor_years = safe_cor(x$apr1_obs, x$apr1_model),
    may1_bias_pct = pct(x$may1_model, x$may1_obs),
    peak_bias_pct = pct(x$peak_model, x$peak_obs),
    meltout_diff_months = mean(x$meltout_diff_months, na.rm = TRUE),
    ppt_novmar_bias_pct = pct(x$ppt_novmar_model, x$ppt_novmar_obs),
    monthly_cor = safe_cor(dd$swe_obs, dd$swe_model),
    monthly_rmse_mm = sqrt(mean((dd$swe_model - dd$swe_obs)^2, na.rm = TRUE)))
}))
summ <- merge(stations[, std_cols], summ, by = "site_id")
summ <- summ[order(summ$source, summ$elev), ]
write.csv(summ, file.path(OUT_DIR, "swe_station_summary.csv"), row.names = FALSE)

# ---- Decade summary (does the calibration hold through time?) ---- #
by_decade <- do.call(rbind, lapply(split(metrics, list(metrics$source, metrics$decade), drop = TRUE), function(x) {
  x <- x[is.finite(x$apr1_obs) & is.finite(x$apr1_model), ]
  if (nrow(x) < 10) return(NULL)
  data.frame(source = x$source[1], decade = x$decade[1], station_years = nrow(x), stations = length(unique(x$site_id)),
    apr1_obs_mm = mean(x$apr1_obs), apr1_model_mm = mean(x$apr1_model), apr1_bias_pct = pct(x$apr1_model, x$apr1_obs),
    apr1_cor = safe_cor(x$apr1_obs, x$apr1_model),
    peak_bias_pct = pct(x$peak_model, x$peak_obs), meltout_diff_months = mean(x$meltout_diff_months, na.rm = TRUE))
}))
if (is.null(by_decade)) by_decade <- data.frame(source = character(), decade = character(), station_years = integer(), stations = integer(),
  apr1_obs_mm = numeric(), apr1_model_mm = numeric(), apr1_bias_pct = numeric(), apr1_cor = numeric(), peak_bias_pct = numeric(), meltout_diff_months = numeric())
by_decade <- by_decade[order(by_decade$source, by_decade$decade), ]
write.csv(by_decade, file.path(OUT_DIR, "swe_by_decade.csv"), row.names = FALSE)

overall <- do.call(rbind, lapply(split(summ, summ$source), function(s) data.frame(
  source = s$source[1], stations = nrow(s), station_years = sum(s$water_years),
  apr1_bias_pct_median = median(s$apr1_bias_pct[is.finite(s$apr1_bias_pct)]), apr1_bias_pct_mean_abs = mean(abs(s$apr1_bias_pct[is.finite(s$apr1_bias_pct)])),
  apr1_cor_years_median = median(s$apr1_cor_years, na.rm = TRUE),
  peak_bias_pct_median = median(s$peak_bias_pct, na.rm = TRUE), meltout_diff_months_mean = mean(s$meltout_diff_months, na.rm = TRUE),
  ppt_novmar_bias_pct_median = median(s$ppt_novmar_bias_pct, na.rm = TRUE), monthly_cor_median = median(s$monthly_cor, na.rm = TRUE))))
write.csv(overall, file.path(OUT_DIR, "swe_overall.csv"), row.names = FALSE)

# ---- Region, elevation and month-of-year breakdowns (where does a bias live?) ---- #
# Regions are coarse longitude/latitude boxes, good enough to separate maritime from continental snow climates.
summ$region <- with(summ, ifelse(longitude < -118 & latitude >= 42, "Cascades/PNW",
                          ifelse(longitude < -116 & latitude < 42, "Sierra/CA-NV",
                          ifelse(longitude >= -116 & latitude >= 44, "N Rockies",
                          ifelse(longitude >= -116 & latitude >= 39, "C Rockies/Wasatch",
                          ifelse(longitude >= -116, "S Rockies/CO-NM-AZ", "Interior NW"))))))
summ$elev_band <- cut(summ$elev, c(-Inf, 1000, 1500, 2000, 2500, 3000, Inf),
                      labels = c("<1000 m", "1000-1500", "1500-2000", "2000-2500", "2500-3000", ">3000 m"))
fmed <- function(x) { x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
group_table <- function(by) do.call(rbind, lapply(split(summ, list(summ$source, summ[[by]]), drop = TRUE), function(s) data.frame(
  source = s$source[1], group = as.character(s[[by]][1]), stations = nrow(s),
  apr1_bias_pct_median = fmed(s$apr1_bias_pct), apr1_cor_years_median = fmed(s$apr1_cor_years),
  peak_bias_pct_median = fmed(s$peak_bias_pct), meltout_diff_months_median = fmed(s$meltout_diff_months),
  ppt_novmar_bias_pct_median = fmed(s$ppt_novmar_bias_pct))))
by_region <- group_table("region"); by_region <- by_region[order(by_region$source, by_region$group), ]
by_elev <- group_table("elev_band"); by_elev <- by_elev[order(by_elev$source, by_elev$group), ]
write.csv(by_region, file.path(OUT_DIR, "swe_by_region.csv"), row.names = FALSE)
write.csv(by_elev, file.path(OUT_DIR, "swe_by_elevation.csv"), row.names = FALSE)
# Mean SWE by month of the water year, all station-months: shows whether SWE is lost in accumulation or in melt.
by_month <- do.call(rbind, lapply(split(d, d$source), function(x) {
  a <- aggregate(cbind(swe_obs, swe_model) ~ month, x, mean); a$source <- x$source[1]
  a$ratio_model_obs <- a$swe_model / a$swe_obs; a[order(wy_order(a$month)), c("source", "month", "swe_obs", "swe_model", "ratio_model_obs")] }))
write.csv(by_month, file.path(OUT_DIR, "swe_by_month.csv"), row.names = FALSE)
# Snow retention: the month's change in SWE divided by the month's precipitation, summed over all
# SNOTEL station-months with gauge data. In months averaging below 0 C the model can only lose
# snow through the rain/snow split, so a model retention below the stations' in Dec-Feb points at
# SNOW_T_ALL / RAIN_T_ALL; a shortfall only in Apr-May points at the melt rate. Pillows gain
# 5-10% more than gauges catch (drifting snow, undercatch), so station values above 1 are normal.
retention <- NULL
dp <- d[d$source == "snotel" & is.finite(d$ppt_obs) & is.finite(d$ppt_model), ]
if (nrow(dp)) {
  dp <- dp[order(dp$site_id, dp$wyear, wy_order(dp$month)), ]
  lag1 <- function(v) c(NA, head(v, -1))
  dp$d_obs <- dp$swe_obs - ave(dp$swe_obs, dp$site_id, dp$wyear, FUN = lag1)
  dp$d_mod <- dp$swe_model - ave(dp$swe_model, dp$site_id, dp$wyear, FUN = lag1)
  dp <- dp[is.finite(dp$d_obs) & is.finite(dp$d_mod) & dp$month %in% c(11, 12, 1, 2, 3, 4, 5), ]
  dp <- merge(dp, summ[, c("site_id", "region")], by = "site_id")
  if (nrow(dp)) {
    retention <- aggregate(cbind(d_obs, d_mod, ppt_obs, ppt_model) ~ region + month, dp, sum)
    retention$station_retained <- retention$d_obs / retention$ppt_obs
    retention$model_retained <- retention$d_mod / retention$ppt_model
    retention <- retention[order(retention$region, wy_order(retention$month)), c("region", "month", "station_retained", "model_retained")]
    write.csv(retention, file.path(OUT_DIR, "swe_retention_by_region_month.csv"), row.names = FALSE)
  }
}

cat("\nOverall by source (bias = model relative to station; melt-out difference in months, + = model later):\n")
print(overall, digits = 3, row.names = FALSE)
cat("\nBy decade:\n"); print(by_decade, digits = 3, row.names = FALSE)
cat("\nBy region (station medians):\n"); print(by_region, digits = 3, row.names = FALSE)
cat("\nBy elevation band (station medians):\n"); print(by_elev, digits = 3, row.names = FALSE)
cat("\nMean SWE by month, all station-months (ratio < 1 = model low):\n"); print(by_month, digits = 3, row.names = FALSE)
if (!is.null(retention)) {
  cat("\nSnow retention (change in SWE / precipitation in the month), SNOTEL, station vs model.",
      "\nDec-Feb shortfall = rain/snow split too warm; Apr-May shortfall only = melt rate:\n")
  for (r in unique(retention$region)) { x <- retention[retention$region == r, ]
    cat(sprintf("%-20s", r), sprintf("%s %5.2f/%5.2f ", month.abb[x$month], x$station_retained, x$model_retained), "\n") }
}
cat("\nStations (first", min(40, nrow(summ)), "of", nrow(summ), "; full table in swe_station_summary.csv):\n")
print(head(summ[, c("site_name", "source", "elev", "water_years", "first_year", "last_year", "apr1_obs_mm", "apr1_model_mm",
                    "apr1_bias_pct", "apr1_cor_years", "peak_bias_pct", "meltout_diff_months", "ppt_novmar_bias_pct")], 40),
      digits = 3, row.names = FALSE)
cat("\nHow to read: compare peak_bias_pct with ppt_novmar_bias_pct. If winter precipitation is already high where",
    "\nsnowpack is high, the excess is precipitation, not phase. If precipitation matches (within ~10%, gauges undercatch)",
    "\nbut snowpack is high, lower SNOW_T_ALL. A consistent meltout_diff_months points at the melt rate.",
    "\nThe decade table shows whether a bias is constant through time or drifts, which would point at the climate inputs.\n")

# ---- Figures ---- #
png(file.path(OUT_DIR, "swe_april1_scatter.png"), 1500, 550, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (src in c("snotel", "course")) {
  m <- metrics[metrics$source == src & is.finite(metrics$apr1_obs) & is.finite(metrics$apr1_model), ]
  if (!nrow(m)) { plot.new(); title(paste("no", src, "data")); next }
  lim <- c(0, max(c(m$apr1_obs, m$apr1_model)))
  plot(m$apr1_obs, m$apr1_model, xlim = lim, ylim = lim, pch = 16, col = adjustcolor("grey20", 0.35), cex = 0.6,
       xlab = "Station 1 April SWE (mm)", ylab = "Model end-of-March SWE (mm)",
       main = sprintf("April 1 SWE, %s (%d station-years)", if (src == "snotel") "SNOTEL" else "snow courses", nrow(m)))
  abline(0, 1, lty = 2)
}
m <- metrics[metrics$source == "snotel", ]
if (any(is.finite(m$meltout_diff_months))) {
  hist(m$meltout_diff_months, breaks = seq(-12.5, 12.5, 1), xlim = c(-6, 6), col = "grey80",
       xlab = "Melt-out difference (months, model - station)", main = "Melt-out timing (SNOTEL)"); abline(v = 0, lty = 2)
} else { plot.new(); title("no melt-out data") }
dev.off()
# Bias by region and elevation, and the mean seasonal cycle
png(file.path(OUT_DIR, "swe_bias_by_region_elevation.png"), 1500, 550, res = 110)
par(mfrow = c(1, 3), mar = c(7, 4, 3, 1))
for (tab in list(by_region, by_elev)) {
  b <- tab[tab$source == "snotel", ]; if (!nrow(b)) b <- tab
  if (!nrow(b)) { plot.new(); next }
  barplot(rbind(b$apr1_bias_pct_median, b$peak_bias_pct_median), beside = TRUE, names.arg = b$group, las = 2,
          col = c("grey40", "firebrick"), ylab = "Median station bias (%)", main = paste("April 1 (grey) and peak (red) SWE bias,", b$source[1]))
  abline(h = 0, lty = 2)
}
par(mar = c(4, 4, 3, 1))
bm <- by_month[by_month$source == "snotel", ]; if (!nrow(bm)) bm <- by_month
if (nrow(bm)) {
  plot(wy_order(bm$month), bm$swe_obs, type = "b", pch = 16, ylim = c(0, max(c(bm$swe_obs, bm$swe_model)) * 1.1), xaxt = "n",
       xlab = "", ylab = "Mean SWE at month end (mm)", main = paste("Seasonal cycle, all", bm$source[1], "station-months"))
  axis(1, at = 1:12, labels = c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S"))
  lines(wy_order(bm$month), bm$swe_model, type = "b", pch = 17, col = "firebrick")
  legend("topleft", c("stations", "model"), pch = c(16, 17), col = c("black", "firebrick"), bty = "n")
}
dev.off()
if (nrow(by_decade)) {
  png(file.path(OUT_DIR, "swe_april1_by_decade.png"), 1000, 550, res = 110)
  par(mar = c(4, 4, 3, 1))
  plot(NA, xlim = range(as.integer(sub("s", "", by_decade$decade))) + c(-5, 15), ylim = c(-50, 50), xlab = "Decade",
       ylab = "April 1 SWE bias, model vs station (%)", main = "Snowpack bias by decade"); abline(h = 0, lty = 2)
  for (src in unique(by_decade$source)) { b <- by_decade[by_decade$source == src, ]
    lines(as.integer(sub("s", "", b$decade)) + 5, b$apr1_bias_pct, type = "b", pch = if (src == "snotel") 16 else 17,
          col = if (src == "snotel") "black" else "firebrick") }
  legend("topright", c("SNOTEL", "snow courses"), pch = c(16, 17), col = c("black", "firebrick"), bty = "n")
  dev.off()
}
clim <- aggregate(cbind(swe_obs, swe_model) ~ site_id + month, d, mean)
clim <- merge(clim, stations[, c("site_id", "site_name", "elev", "source")])
clim$wm <- wy_order(clim$month)
show <- summ[order(summ$source, -summ$water_years), ]; show <- show[!duplicated(show$site_id), ]
show <- rbind(head(show[show$source == "snotel", ], 12), head(show[show$source == "course", ], 12))
if (nrow(show)) {
  nc <- min(4, nrow(show)); nr <- ceiling(nrow(show) / nc)
  png(file.path(OUT_DIR, "swe_climatology_by_station.png"), 400 * nc, 320 * nr + 60, res = 110)
  par(mfrow = c(nr, nc), mar = c(3.5, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (s in show$site_id) {
    x <- clim[clim$site_id == s, ]; x <- x[order(x$wm), ]
    if (!any(is.finite(c(x$swe_obs, x$swe_model)))) { plot.new(); title(x$site_name[1]); next }
    plot(x$wm, x$swe_obs, type = "b", pch = 16, ylim = c(0, max(c(x$swe_obs, x$swe_model), na.rm = TRUE) * 1.1), xlim = c(1, 12),
         xaxt = "n", xlab = "", ylab = "SWE at month end (mm)", main = paste0(x$site_name[1], " (", x$elev[1], " m, ", x$source[1], ")"))
    axis(1, at = 1:12, labels = c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S"))
    lines(x$wm, x$swe_model, type = "b", pch = 17, col = "firebrick")
  }
  mtext("Mean end-of-month SWE: station (black) vs model (red); longest records shown", outer = TRUE, cex = 1)
  dev.off()
}
cat("\nWritten to", OUT_DIR, "\n")
