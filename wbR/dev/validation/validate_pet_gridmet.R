############################################
######## Validate water-balance PET ########
######## against gridMET reference ET ######

# Compares the model's monthly PET (the PET monthly rasters of v5.5, or AET + CWD from older runs) with FAO-56 / ASCE Penman-Monteith grass reference ET (ETo)
# computed from the gridMET daily inputs (tmmx, tmmn, srad, rmin, rmax, vs). gridMET's own "pet" product uses the same formulation and inputs.

# Also computes a Priestley-Taylor ET from gridMET's observed radiation and  humidity (PT_obs). Model PET vs PT_obs tests the model's radiation, transmittance and humidity normals; PT_obs vs ETo is the expected gap between a radiation-only formulation and the full Penman-Monteith in dry air.

# Model rasters are averaged up to the 4-km gridMET cells, so terrain cancels and the comparison is about level and seasonality. The comparison covers the overlap of the model area and the gridMET files: a West-wide run checked against gridMET files that cover the Malheur is compared over the Malheur.

# Everything is written to <OUTPUT_ROOT>/validation_gridmet/.

# Works on any run with monthly "pet" (v5.5) in MONTHLY_OUTPUT. Requires terra. Runtime: 1-2 minutes per year.

suppressPackageStartupMessages(library(terra))
terra::projNetwork(FALSE)

# ---- Settings ---- #
BASE_DIR <- "/YOUR PATH"
OUTPUT_ROOT <- file.path(BASE_DIR, "CWD_outputs")                      # the water-balance output folder
GRIDMET_DIR <- file.path(BASE_DIR, "FireClimateWeather_MODIS500m_gridMET/cache/gridmet_nc")
DEM_PATH <- file.path(BASE_DIR, "western_us_dem_500m.tif")              # same DEM as the model run (mal_dem.tif for a Malheur run)
DEM_Z_UNITS <- "auto"                                                   # "feet", "meters", or "auto" (feet if the maximum exceeds 5,000)
YEARS <- 2001:2025                                                      # years with monthly model files and gridMET inputs
PT_ALPHA <- 1.26                                                        # as in the model run
ELEV_BANDS <- c(0, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 4500)  # m, for the elevation tables (empty bands are dropped)
OUT_DIR <- file.path(OUTPUT_ROOT, "validation_gridmet")

# ---- Helpers ---- #
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
month_of <- function(d) as.integer(format(d, "%m"))
es_kpa <- function(Tc) 0.6108 * exp(17.27 * Tc / (Tc + 237.3))
model_monthly <- function(year) {
  # v5.5 writes PET directly; older runs (or runs without "pet") give PET = AET + CWD, an exact identity in the model.
  f_pet <- file.path(OUTPUT_ROOT, "monthly", paste0("PET_monthly_mm_", year, ".tif"))
  if (file.exists(f_pet)) return(rast(f_pet))
  f <- file.path(OUTPUT_ROOT, "monthly", paste0(c("CWD", "AET"), "_monthly_mm_", year, ".tif"))
  if (!all(file.exists(f))) return(NULL)
  rast(f[1]) + rast(f[2])
}
gm <- function(var, year) {
  f <- file.path(GRIDMET_DIR, paste0(var, "_", year, ".nc"))
  if (!file.exists(f)) stop("Missing gridMET file: ", f)
  rast(f)
}

# ---- Grids ---- #
first <- NULL
for (y in YEARS) { first <- model_monthly(y); if (!is.null(first)) break }
if (is.null(first)) stop("No monthly PET (or CWD + AET) files found under ", file.path(OUTPUT_ROOT, "monthly"),
                         ". Run the model with \"pet\" (or \"cwd\" and \"aet\") in MONTHLY_OUTPUT.")
model_crs <- crs(first)
# gridMET window: the overlap of the model extent (in lon/lat, plus one cell margin) and the gridMET files.
e <- ext(project(as.polygons(ext(first), crs = model_crs), "EPSG:4326"))
g1 <- gm("tmmx", YEARS[1])[[1]]
win <- intersect(ext(e) + 0.05, ext(g1))
if (is.null(win)) stop("The gridMET files do not overlap the model area.")
tmpl <- crop(g1, win, snap = "out")
tmpl <- rast(tmpl); values(tmpl) <- 1
cat("gridMET window:", nrow(tmpl), "x", ncol(tmpl), "cells;",
    if (win == ext(e) + 0.05) "the gridMET files cover the whole model area" else "the gridMET files cover part of the model area (comparison limited to that part)", "\n")

# Elevation per gridMET cell (mean of the DEM), latitude per cell.
dem <- crop(rast(DEM_PATH), project(as.polygons(ext(tmpl) + 0.1, crs = "EPSG:4326"), crs(rast(DEM_PATH))), snap = "out")
fact <- max(1L, round(100 / mean(res(dem)))); if (fact > 1) dem <- aggregate(dem, fact = fact, fun = "mean", na.rm = TRUE)
if (DEM_Z_UNITS == "auto") DEM_Z_UNITS <- if (isTRUE(global(dem, "max", na.rm = TRUE)[1, 1] > 5000)) "feet" else "meters"
if (DEM_Z_UNITS == "feet") dem <- dem * 0.3048
cat("DEM elevation units taken as", DEM_Z_UNITS, "\n")
elev <- project(dem, tmpl, method = "average")
lat <- init(tmpl, "y")
# Only gridMET cells substantially covered by the model (>= 50% of cell area).
cover <- project(crop(!is.na(first[[1]]), project(as.polygons(ext(tmpl) + 0.1, crs = "EPSG:4326"), model_crs), snap = "out"), tmpl, method = "average")
keep <- cover >= 0.5
cat("gridMET cells with >= 50% model coverage:", sum(values(keep), na.rm = TRUE), "\n")

z <- values(elev, mat = FALSE); phi <- values(lat, mat = FALSE) * pi / 180
P_kpa <- 101.3 * ((293 - 0.0065 * z) / 293)^5.26
gam <- 0.000665 * P_kpa

# ---- Daily ETo and PT from gridMET, summed to months ---- #
eto_year <- function(year) {
  tx <- values(crop(gm("tmmx", year), tmpl)) - 273.15
  tn <- values(crop(gm("tmmn", year), tmpl)) - 273.15
  rs <- values(crop(gm("srad", year), tmpl)) * 0.0864      # W m-2 -> MJ m-2 d-1
  rhx <- values(crop(gm("rmax", year), tmpl)); rhn <- values(crop(gm("rmin", year), tmpl))
  u10 <- values(crop(gm("vs", year), tmpl)); u2 <- u10 * 4.87 / log(67.8 * 10 - 5.42)
  dates <- time(gm("tmmx", year)); doy <- as.integer(format(dates, "%j"))
  nd <- ncol(tx); nc <- nrow(tx)
  tm <- (tx + tn) / 2
  es <- (es_kpa(tx) + es_kpa(tn)) / 2
  ea <- (es_kpa(tn) * pmin(rhx, 100) + es_kpa(tx) * pmin(rhn, 100)) / 200
  ea <- pmin(ea, es)
  delta <- 4098 * es_kpa(tm) / (tm + 237.3)^2
  # Extraterrestrial and clear-sky radiation (FAO-56 eq. 21, 37)
  dr <- 1 + 0.033 * cos(2 * pi * doy / 365); dec <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws <- acos(pmax(-1, pmin(1, -outer(tan(phi), tan(dec)))))
  ra <- 24 * 60 / pi * 0.0820 * matrix(dr, nc, nd, byrow = TRUE) *
        (ws * outer(sin(phi), sin(dec)) + outer(cos(phi), cos(dec)) * sin(ws))
  rso <- (0.75 + 2e-5 * z) * ra
  ratio <- pmax(0.3, pmin(1, rs / rso))
  rnl <- 4.903e-9 * ((tx + 273.16)^4 + (tn + 273.16)^4) / 2 * (0.34 - 0.14 * sqrt(ea)) * (1.35 * ratio - 0.35)
  rn <- 0.77 * rs - rnl
  G <- matrix(gam, nc, nd)
  eto <- (0.408 * delta * rn + G * 900 / (tm + 273) * u2 * (es - ea)) / (delta + G * (1 + 0.34 * u2))
  pt <- PT_ALPHA * delta / (delta + G) * rn / 2.45
  eto[eto < 0] <- 0; pt[pt < 0] <- 0
  m <- month_of(dates)
  list(eto = sapply(1:12, function(k) rowSums(eto[, m == k, drop = FALSE])),
       pt = sapply(1:12, function(k) rowSums(pt[, m == k, drop = FALSE])))
}

# ---- Loop over years ---- #
rows <- list()
for (year in YEARS) {
  mod <- model_monthly(year)
  if (is.null(mod)) { cat("No model monthly files for", year, "- skipped\n"); next }
  cat("Year", year, "... ")
  g <- eto_year(year)
  mod <- crop(mod, project(as.polygons(ext(tmpl) + 0.1, crs = "EPSG:4326"), model_crs), snap = "out")  # only the window, not the whole run
  pm <- values(project(mod, tmpl, method = "average"))   # model PET on the gridMET grid, 12 columns
  ok <- which(values(keep, mat = FALSE) == 1 & is.finite(z))
  rows[[length(rows) + 1]] <- data.frame(
    cell = rep(ok, 12), year = year, month = rep(1:12, each = length(ok)),
    elev_m = rep(z[ok], 12),
    pet_model = as.vector(pm[ok, ]), eto_pm = as.vector(g$eto[ok, ]), pt_obs = as.vector(g$pt[ok, ]))
  cat("done\n")
}
d <- do.call(rbind, rows)
d <- d[is.finite(d$pet_model) & is.finite(d$eto_pm), ]
write.csv(d, file.path(OUT_DIR, "pet_vs_gridmet_cell_month.csv"), row.names = FALSE)

# ---- Tables ---- #
agg <- function(by, dd = d) {
  a <- aggregate(cbind(pet_model, eto_pm, pt_obs) ~ ., data = cbind(dd[c("pet_model", "eto_pm", "pt_obs")], by), FUN = mean)
  a$model_over_eto <- a$pet_model / a$eto_pm
  a$model_over_pt_obs <- a$pet_model / a$pt_obs
  a$pt_obs_over_eto <- a$pt_obs / a$eto_pm
  a
}
monthly <- agg(data.frame(month = d$month))
d$elev_band <- droplevels(cut(d$elev_m, ELEV_BANDS))
warm <- d[d$month %in% 4:9, ]
ann <- aggregate(cbind(pet_model, eto_pm, pt_obs) ~ cell + year + elev_band, d, sum)
by_elev <- agg(data.frame(elev_band = ann$elev_band), ann)
by_elev_warm <- agg(data.frame(elev_band = warm$elev_band), aggregate(cbind(pet_model, eto_pm, pt_obs) ~ cell + year + elev_band, warm, sum))
write.csv(monthly, file.path(OUT_DIR, "pet_vs_gridmet_by_month.csv"), row.names = FALSE)
write.csv(by_elev, file.path(OUT_DIR, "pet_vs_gridmet_by_elevation_annual.csv"), row.names = FALSE)
write.csv(by_elev_warm, file.path(OUT_DIR, "pet_vs_gridmet_by_elevation_AprSep.csv"), row.names = FALSE)

cat("\nYears used:", paste(range(d$year), collapse = "-"), "; gridMET cells:", length(unique(d$cell)), "\n")
cat("\nMonthly means over all cells and years (mm):\n"); print(monthly, digits = 3, row.names = FALSE)
cat("\nAnnual totals by elevation band (mm):\n"); print(by_elev, digits = 3, row.names = FALSE)
cat("\nApril-September totals by elevation band (mm):\n"); print(by_elev_warm, digits = 3, row.names = FALSE)
cat("\nHow to read: model_over_eto is the model's PET as a fraction of Penman-Monteith grass reference ET.",
    "\nA radiation-only PET is expected to sit below ETo in dry air (roughly 0.75-0.90 in the warm season).",
    "\nmodel_over_pt_obs isolates the model's radiation/transmittance/humidity normals: near 1.0 means they",
    "\nreproduce what gridMET observed. Drift of either ratio with month or elevation is the signal to act on.\n")

# ---- Figure ---- #
png(file.path(OUT_DIR, "pet_vs_gridmet.png"), 1500, 1000, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(monthly$month, monthly$eto_pm, type = "b", pch = 16, ylim = c(0, max(monthly$eto_pm) * 1.1),
     xlab = "Month", ylab = "mm per month", main = "Area-mean monthly ET demand")
lines(monthly$month, monthly$pet_model, type = "b", pch = 17, col = "firebrick")
lines(monthly$month, monthly$pt_obs, type = "b", pch = 15, col = "steelblue")
legend("topleft", c("gridMET Penman-Monteith ETo", "Model PET (Priestley-Taylor)", "PT from gridMET radiation"),
       col = c("black", "firebrick", "steelblue"), pch = c(16, 17, 15), bty = "n")
plot(monthly$month, monthly$model_over_eto, type = "b", pch = 17, col = "firebrick", ylim = c(0, 1.6),
     xlab = "Month", ylab = "ratio", main = "Model PET / gridMET ETo (red), / PT_obs (blue)")
lines(monthly$month, monthly$model_over_pt_obs, type = "b", pch = 15, col = "steelblue"); abline(h = 1, lty = 2)
band_mid <- function(lv) sapply(strsplit(gsub("[][()]", "", as.character(lv)), ","), function(x) mean(as.numeric(x)))
mids <- band_mid(by_elev$elev_band)
plot(mids, by_elev$eto_pm, type = "b", pch = 16, ylim = c(0, max(by_elev$eto_pm) * 1.1),
     xlab = "Elevation band midpoint (m)", ylab = "mm per year", main = "Annual ET demand by elevation")
lines(mids, by_elev$pet_model, type = "b", pch = 17, col = "firebrick")
lines(mids, by_elev$pt_obs, type = "b", pch = 15, col = "steelblue")
a <- aggregate(cbind(pet_model, eto_pm) ~ cell + year, d, sum)
plot(a$eto_pm, a$pet_model, pch = ".", col = adjustcolor("grey30", 0.4), xlab = "gridMET ETo, annual (mm)",
     ylab = "Model PET, annual (mm)", main = "Cell-years"); abline(0, 1, lty = 2)
abline(lm(pet_model ~ eto_pm, a), col = "firebrick")
dev.off()
cat("\nWritten to", OUT_DIR, "\n")
