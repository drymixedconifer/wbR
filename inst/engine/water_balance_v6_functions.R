##############################################
##############################################
############### Functions to #################
######## create water balance rasters ########
##############################################
##############################################

# Companion to water_balance_v6.0.R. Keep these two R files together.
# C++ is embedded below: no separate .cpp file is needed.
# Dependencies: terra, prism, Rcpp, RcppParallel; an R-compatible C++ compiler.
# This file defines functions; sourcing it does not launch the analysis.

# Some notes
#  * Radiation. The single McCune-Keon heat-load scalar of Hamon is replaced by a 12-month clear-sky shortwave radiation term computed for every cell from its DEM slope, aspect, latitude and elevation (direct beam integrated over the day on the tilted surface, isotropic diffuse scaled by sky view, terrain-reflected diffuse). No horizon shading yet; that may come later in a fine-resolution option. Heat load is still computed and reported, and is still used by the optional Hamon method.
#  * PET. Priestley-Taylor (default) driven by that radiation, PRISM monthly cloud-transmittance normals (soltrans) and PRISM dew-point normals (tdmean) for the longwave term; or Hamon (PET_METHOD = "hamon"), which reproduces v4.
#  * Static normals. soltrans and tdmean are 1991-2020 monthly normals, applied every year; each can be replaced by a constant (see the driver).
#  * PET self-check. Optionally writes Hamon vs Priestley-Taylor monthly PET over the normals so a new study area can be sanity-checked.
#  * Memory. The per-cell monthly radiation table adds 12 doubles per modeled cell (about 1.9 GB for the 11 western states at 400 m).
# v5.5: points (POINTS_PATH): monthly and annual values at user points written as CSV, from a full run or, with WRITE_RASTERS = FALSE, from a run on the point cells only. Four more monthly variables: tmean, rs, rnet, melt.
#    Snow defaults re-fitted to SNOTEL station-months (SNOW_T_ALL -6, RAIN_T_ALL 10, energy melt 2.3 mm/C/day + 2.3 mm/MJ, snow albedo 0.77); no engine change.
#    Monthly variables are written in chunks (MONTHLY_CHUNK) so memory does not grow with the number requested; the year is recomputed per chunk with identical results.
# v6.0: Penman-Monteith PET (PET_METHOD = "penman_monteith", the default). The radiation term is the one Priestley-Taylor
#    uses (clear-sky on the slope x PRISM transmittance normal, FAO-56 longwave with the dew-point normal) with alpha 1 and the
#    wind-dependent denominator; the aerodynamic term uses PRISM monthly vapour-pressure deficit ((vpdmax + vpdmin)/2, downloaded
#    like ppt and tmean, 1948-2025 plus the 1991-2020 normals for spin-up) downscaled to each cell by conserving relative
#    humidity through the temperature lapse, and a static 2 m wind climatology (WIND_PATH, 12 monthly layers on any grid) or a
#    constant (WIND_CONSTANT). AERO_SCALE multiplies the aerodynamic term (1 = FAO-56 as written). "priestley_taylor" reproduces
#    5.5 exactly and "hamon" 4.0; the three share one radiation chain. New monthly variable vpd (kPa). The transmittance and
#    dew-point normals are now prepared under every PET method, so the energy melt option is defined under Hamon too.

# Context follows nested processing calls, including lazy raster reads.
.wb_log <- new.env(parent = emptyenv())
.wb_log$context <- list(step = "startup", year = "not applicable", files = "not applicable")
.wb_log$report <- NULL
wb_enter <- function(step, year = NULL, files = NULL) {
  old <- .wb_log$context
  .wb_log$context$step <- step
  if (!is.null(year)) .wb_log$context$year <- paste(year, collapse = ", ")
  if (!is.null(files)) .wb_log$context$files <- paste(unique(files), collapse = " | ")
  old
}
wb_note <- function(text) {
  if (!is.null(.wb_log$report)) cat(paste(text, collapse = "\n"), "\n",
                                   file = .wb_log$report, append = TRUE)
}
wb_condition <- function(c, level) {
  ctx <- .wb_log$context
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " [", level,
    "] ", conditionMessage(c), "\n  step=", ctx$step, "; year=", ctx$year,
    "\n  call=", paste(deparse(conditionCall(c)), collapse = " "),
    "\n  file(s)=", ctx$files)
  wb_note(msg)
  cat(msg, "\n", file = stderr())
}
# Delete only this run's unreferenced terra files, never inputs/checkpoints.
wb_cleanup_temp <- function(keep = list()) {
  td <- .wb_log$temp_dir
  if (is.null(td) || !dir.exists(td)) return(invisible(NULL))
  referenced <- function(x) {
    if (inherits(x, "SpatRaster")) return(sources(x))
    if (is.list(x)) return(unlist(lapply(x, referenced), use.names = FALSE))
    character()
  }
  live <- unique(referenced(keep))
  live <- normalizePath(live[nzchar(live)], mustWork = FALSE)
  gc(verbose = FALSE)
  files <- list.files(td, full.names = TRUE)
  raster <- files[grepl("\\.tif$", files, ignore.case = TRUE)]
  dead <- raster[!normalizePath(raster, mustWork = FALSE) %in% live]
  unlink(c(dead, paste0(dead, ".aux.xml")))
  invisible(NULL)
}
wb_check_space <- function(path, minimum_gib = 10) {
  if (Sys.info()[["sysname"]] != "Darwin") return(invisible(NULL))
  z <- system2("df", c("-kP", shQuote(path)), stdout = TRUE, stderr = FALSE)
  fields <- strsplit(trimws(tail(z, 1)), "[[:space:]]+")[[1]]
  free <- suppressWarnings(as.numeric(fields[4])) / 1024^2
  if (!is.finite(free)) stop("Could not check available disk space for ", path)
  if (free < minimum_gib) stop(sprintf("Low disk space: %.1f GiB free at %s; at least %.1f GiB required before the next step.", free, path, minimum_gib))
  invisible(free)
}

wb_step <- function(step, files, expr) {
  old <- wb_enter(step, files = files)
  on.exit(.wb_log$context <- old)
  force(expr)
}

wb_load_packages <- function() {
  pkgs <- c("terra", "prism", "Rcpp", "RcppParallel")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(terra))
}

wb_sanitize_label <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) x <- "aoi"
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

wb_month_tag <- function(months_vec) {
  month_abbr_vec <- c(
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  )
  
  months_vec <- sort(unique(as.integer(months_vec)))
  
  if (length(months_vec) == 12L && identical(months_vec, 1:12)) {
    return("Jan_Dec")
  }
  
  if (length(months_vec) > 1L && all(diff(months_vec) == 1L)) {
    return(paste0(
      month_abbr_vec[min(months_vec)],
      "_",
      month_abbr_vec[max(months_vec)]
    ))
  }
  
  paste(month_abbr_vec[months_vec], collapse = "_")
}

wb_infer_dem_z_units <- function(dem_stats) {
  
  zmin  <- as.numeric(dem_stats[1, "min"])
  zmax  <- as.numeric(dem_stats[1, "max"])
  zmean <- as.numeric(dem_stats[1, "mean"])
  
  if (!all(is.finite(c(zmin, zmax, zmean)))) {
    stop("DEM summary statistics are not finite; cannot infer z-units.")
  }
  
  if (zmax <= 4500 && zmean <= 3000) {
    return("meters")
  }
  
  if (zmax > 4500 || zmean > 3000) {
    return("feet")
  }
  
  stop(
    "Could not infer DEM z-units. Set cfg$dem_z_units to 'meters' or 'feet'."
  )
}

wb_daylength <- function(lat_deg) {
  
  doy_mid <- c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  lat_rad <- lat_deg * pi / 180
  
  decl_rad <- 0.409 * sin(2 * pi * (doy_mid - 80) / 365)
  omega_s  <- acos(pmax(-1, pmin(1, -tan(lat_rad) * tan(decl_rad))))
  
  24 / pi * omega_s
}

wb_latitude_raster <- function(dem) {
  template <- rast(dem)
  from <- crs(template)
  latitude <- app(init(template, "cell"), function(cells) {
    xy <- xyFromCell(template, as.vector(cells))
    ll <- terra::project(xy, from, "EPSG:4326")
    lat <- ll[,2]
    if (any(!is.finite(lat)) || any(abs(lat) > 90))
      stop("Could not transform DEM cell centers to valid latitudes.")
    lat
  })
  names(latitude) <- "latitude_degrees"
  mask(latitude, dem)
}

wb_make_12 <- function(x) {
  do.call(c, rep(list(x), 12))
}

wb_cpp_source <- r"---(
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace RcppParallel;

static const double WB_PI = 3.14159265358979323846;
static const double GSC = 0.0820;          // solar constant, MJ m-2 min-1 (FAO-56)
static const double SIGMA = 4.903e-9;      // Stefan-Boltzmann, MJ K-4 m-2 d-1
static const double LAMBDA = 2.45;         // latent heat of vaporization, MJ kg-1
static const double MID_DOY[12] = {15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349};

// ---------------------------------------------------------------------------
// Clear-sky shortwave radiation on each cell's tilted surface, by month.
//
// For the mid-month day the direct beam is integrated through the daylight
// hours on both the horizontal plane (Ra_h) and the tilted plane (Ra_s) with
// the standard incidence-angle geometry (Duffie & Beckman). Daylight is where
// the sun is above the horizontal plane: no horizon shading by surrounding
// terrain. The clear-sky atmosphere transmits a fraction 0.75 + 2e-5 z
// (FAO-56), of which `diffuse_frac` arrives as isotropic diffuse light. On the
// slope: beam scaled by Ra_s/Ra_h, diffuse by the sky-view fraction
// (1 + cos b)/2, plus ground-reflected light albedo * (1 - cos b)/2.
// Level ground reproduces the FAO-56 clear-sky value exactly.
// ---------------------------------------------------------------------------
struct RadWorker : public Worker {
  const RVector<double> slope;     // radians
  const RVector<double> aspect;    // degrees clockwise from north
  const RVector<double> lat;       // degrees
  const RVector<double> elev;      // m
  const RVector<double> decl;      // 12 solar declinations, radians
  const double diffuse_frac;
  const double albedo;
  const double step_min;
  const bool want_flat;
  RMatrix<double> out;             // cells x 12, MJ m-2 d-1 (clear sky, slope)
  RMatrix<double> flat;            // cells x 12, MJ m-2 d-1 (clear sky, level); 1 x 12 dummy if !want_flat

  RadWorker(const NumericVector slope, const NumericVector aspect, const NumericVector lat,
            const NumericVector elev, const NumericVector decl, double diffuse_frac,
            double albedo, double step_min, bool want_flat, NumericMatrix out, NumericMatrix flat) :
    slope(slope), aspect(aspect), lat(lat), elev(elev), decl(decl),
    diffuse_frac(diffuse_frac), albedo(albedo), step_min(step_min), want_flat(want_flat), out(out), flat(flat) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int nstep = (int) std::ceil(1440.0 / step_min);
    const double dw = 2.0 * WB_PI / nstep;              // hour-angle step, radians
    for (std::size_t i = begin; i < end; i++) {
      const double phi = lat[i] * WB_PI / 180.0;
      double b = slope[i];
      if (!std::isfinite(b) || b < 0.0) b = 0.0;
      if (b > WB_PI / 2.0) b = WB_PI / 2.0;
      double a = aspect[i];
      if (!std::isfinite(a)) a = 0.0;
      // Surface azimuth: 0 due south, west positive, east negative.
      const double gamma = (a - 180.0) * WB_PI / 180.0;
      double z = elev[i];
      if (!std::isfinite(z)) z = 0.0;
      z = std::max(0.0, std::min(z, 5000.0));
      const double tau = 0.75 + 2e-5 * z;
      const double sinphi = std::sin(phi), cosphi = std::cos(phi);
      const double sinb = std::sin(b), cosb = std::cos(b);
      const double sing = std::sin(gamma), cosg = std::cos(gamma);
      const bool ok = std::isfinite(phi) && std::abs(lat[i]) <= 90.0;
      for (int m = 0; m < 12; m++) {
        if (!ok) { out(i, m) = NA_REAL; if (want_flat) flat(i, m) = NA_REAL; continue; }
        const double d = decl[m];
        const double dr = 1.0 + 0.033 * std::cos(2.0 * WB_PI * MID_DOY[m] / 365.0);
        const double sind = std::sin(d), cosd = std::cos(d);
        double ra_h = 0.0, ra_s = 0.0;
        for (int k = 0; k < nstep; k++) {
          const double w = -WB_PI + (k + 0.5) * dw;        // hour angle, midpoint rule
          const double cosz = sinphi * sind + cosphi * cosd * std::cos(w);
          if (cosz <= 0.0) continue;                        // sun below the horizontal plane
          ra_h += cosz;
          const double cost = sind * sinphi * cosb
                            - sind * cosphi * sinb * cosg
                            + cosd * cosphi * cosb * std::cos(w)
                            + cosd * sinphi * sinb * cosg * std::cos(w)
                            + cosd * sinb * sing * std::sin(w);
          if (cost > 0.0) ra_s += cost;
        }
        // Each step spans dw radians = dw * 720/pi minutes.
        const double minutes = dw * 720.0 / WB_PI;
        ra_h *= GSC * dr * minutes;
        ra_s *= GSC * dr * minutes;
        const double rso_h = tau * ra_h;
        double rso_s = tau * ((1.0 - diffuse_frac) * ra_s + diffuse_frac * ra_h * (1.0 + cosb) / 2.0)
                     + albedo * rso_h * (1.0 - cosb) / 2.0;
        if (!std::isfinite(rso_s) || rso_s < 0.0) rso_s = 0.0;
        out(i, m) = rso_s;
        if (want_flat) flat(i, m) = rso_h;
      }
    }
  }
};

// [[Rcpp::export]]
List wb_radiation_cpp(NumericVector slope, NumericVector aspect, NumericVector lat,
                      NumericVector elev, NumericVector declination,
                      double diffuse_frac, double albedo, double step_min, bool want_flat) {
  const R_xlen_t n = slope.size();
  if (aspect.size() != n || lat.size() != n || elev.size() != n)
    stop("Terrain vectors must describe the same cells.");
  if (declination.size() != 12) stop("Expected 12 monthly declinations.");
  if (!(diffuse_frac >= 0.0 && diffuse_frac <= 1.0) || !(albedo >= 0.0 && albedo <= 1.0) ||
      !(step_min > 0.0 && step_min <= 60.0))
    stop("Invalid radiation settings.");
  NumericMatrix out(n, 12), flat(want_flat ? n : 1, 12);
  RadWorker worker(slope, aspect, lat, elev, declination, diffuse_frac, albedo, step_min, want_flat, out, flat);
  parallelFor(0, n, worker);
  return List::create(_["slope"] = out, _["flat"] = flat);
}


// ---------------------------------------------------------------------------
// One simulated year for every cell, entirely in memory. For each cell:
// interpolate the 12 monthly values from the coarse climate grid, apply the
// elevation adjustments, compute PET, then run the monthly water balance.
// ---------------------------------------------------------------------------
struct WBYearWorker : public Worker {

  const RMatrix<double> Pc;        // coarse precipitation, cells x 12 (row-major cells)
  const RMatrix<double> Tc;        // coarse temperature, cells x 12
  const RMatrix<double> Vc;        // coarse vapour-pressure deficit, kPa, cells x 12 (1 x 12 of NA unless Penman-Monteith) (v6.0)
  const int nrc;
  const int ncc;
  const RVector<double> rowf;      // each cell's fractional row on the coarse grid
  const RVector<double> colf;      // each cell's fractional column on the coarse grid
  const RVector<int> covered;      // 1 where the coarse grid reaches the cell
  const RVector<double> dz;        // cell elevation minus coarse-grid elevation, m
  const RMatrix<double> lapse_t;   // temperature slope, C per m, cells x 12
  const RMatrix<double> beta_p;    // log-precipitation slope, per m, cells x 12
  const RVector<double> awcv;
  const RVector<double> hlv;       // heat load (Hamon only)
  const RVector<double> soilv;
  const RVector<double> packv;
  const RVector<double> latv;
  const RVector<double> elevv;     // m, for atmospheric pressure
  const RMatrix<double> radv;      // clear-sky radiation on the slope, cells x 12, MJ m-2 d-1
  // Static normals on their own coarse grid (cloud transmittance, dew point).
  const RMatrix<double> Sc;        // transmittance, static cells x 12 (may be 1 x 12 of NA)
  const RMatrix<double> Ec;        // dew point C, static cells x 12
  const int nrs;
  const int ncs;
  const RVector<double> rowf_s;
  const RVector<double> colf_s;
  const RVector<int> covered_s;
  const RVector<double> trans_fill; // 12 monthly fallback transmittances
  const double rh_fill;             // fallback relative humidity (fraction)
  // Wind climatology on its own grid (v6.0): 2 m wind speed, m s-1, grid cells x 12; 1 x 12 of NA when a constant is used.
  const RMatrix<double> Wc;
  const int nrw;
  const int ncw;
  const RVector<double> rowf_w;
  const RVector<double> colf_w;
  const RVector<int> covered_w;
  const double wind_fill;           // 2 m wind speed where the climatology does not reach (or everywhere, when constant)
  const double aero_scale;          // multiplier on the aerodynamic term (1 = FAO-56 as written)
  const RVector<int> keep;
  const RVector<double> declination;
  const RVector<double> month_days;
  const int method;                 // 1 = Hamon x heat load (v4), 2 = Priestley-Taylor, 3 = Penman-Monteith (v6.0)
  const double pt_alpha;
  const double albedo;
  const int diag;                   // 0: none; 1: annual PET column; 2: 12 monthly PET columns
  // Snow (v5.3)
  const double snow_t_all;          // monthly mean T (C) at or below which precipitation is all snow
  const double rain_t_all;          // monthly mean T (C) at or above which precipitation is all rain
  const double melt_factor;         // mm of melt per degree-day above 0 C on level ground
  const double melt_rad_weight;     // share of the melt factor scaled by the cell's radiation relative to level ground (0 = v5.2)
  // Snow (v5.4)
  const int melt_method;            // 1 = degree-day (v5.3), 2 = energy: c_T * T + c_R * (1 - snow_albedo) * Rs
  const double melt_rad_factor;     // mm of melt per MJ m-2 of shortwave absorbed by the snow (energy method)
  const double snow_albedo;         // snow surface albedo for the energy method
  const int monthly;                // bitmask of monthly outputs: 1 cwd, 2 aet, 4 pet, 8 water, 16 soil, 32 pack, 64 ppt, 128 snowfall,
                                    //   256 tmean, 512 rs, 1024 rnet, 2048 melt (v5.5), 4096 vpd (v6.0)
  const int mcol0;                  // first monthly output column

  RMatrix<double> out;

  WBYearWorker(
    const NumericMatrix Pc, const NumericMatrix Tc, const NumericMatrix Vc, int nrc, int ncc,
    const NumericVector rowf, const NumericVector colf, const IntegerVector covered,
    const NumericVector dz, const NumericMatrix lapse_t, const NumericMatrix beta_p,
    const NumericVector awcv, const NumericVector hlv, const NumericVector soilv,
    const NumericVector packv, const NumericVector latv, const NumericVector elevv,
    const NumericMatrix radv, const NumericMatrix Sc, const NumericMatrix Ec, int nrs, int ncs,
    const NumericVector rowf_s, const NumericVector colf_s, const IntegerVector covered_s,
    const NumericVector trans_fill, double rh_fill,
    const NumericMatrix Wc, int nrw, int ncw, const NumericVector rowf_w, const NumericVector colf_w,
    const IntegerVector covered_w, double wind_fill, double aero_scale, const IntegerVector keep,
    const NumericVector declination, const NumericVector month_days, int method,
    double pt_alpha, double albedo, int diag, double snow_t_all, double rain_t_all,
    double melt_factor, double melt_rad_weight, int melt_method, double melt_rad_factor, double snow_albedo,
    int monthly, int mcol0, NumericMatrix out
  ) :
    Pc(Pc), Tc(Tc), Vc(Vc), nrc(nrc), ncc(ncc), rowf(rowf), colf(colf), covered(covered),
    dz(dz), lapse_t(lapse_t), beta_p(beta_p), awcv(awcv), hlv(hlv), soilv(soilv),
    packv(packv), latv(latv), elevv(elevv), radv(radv), Sc(Sc), Ec(Ec), nrs(nrs), ncs(ncs),
    rowf_s(rowf_s), colf_s(colf_s), covered_s(covered_s), trans_fill(trans_fill),
    rh_fill(rh_fill), Wc(Wc), nrw(nrw), ncw(ncw), rowf_w(rowf_w), colf_w(colf_w), covered_w(covered_w),
    wind_fill(wind_fill), aero_scale(aero_scale), keep(keep), declination(declination), month_days(month_days),
    method(method), pt_alpha(pt_alpha), albedo(albedo), diag(diag), snow_t_all(snow_t_all), rain_t_all(rain_t_all),
    melt_factor(melt_factor), melt_rad_weight(melt_rad_weight), melt_method(melt_method),
    melt_rad_factor(melt_rad_factor), snow_albedo(snow_albedo), monthly(monthly), mcol0(mcol0), out(out) {}

  // Bilinear weights of the four coarse cells around (r, c); neighbours outside
  // the grid are dropped. Returns the number of usable neighbours.
  static int neighbours(double r, double c, int nr, int nc, int* idx, double* w) {
    const int r0 = (int) std::floor(r);
    const int c0 = (int) std::floor(c);
    const double fy = r - r0;
    const double fx = c - c0;
    int k = 0;
    for (int dr = 0; dr < 2; dr++) {
      for (int dc = 0; dc < 2; dc++) {
        const int rr = r0 + dr;
        const int cc = c0 + dc;
        const double ww = (dr ? fy : 1.0 - fy) * (dc ? fx : 1.0 - fx);
        if (rr >= 0 && rr < nr && cc >= 0 && cc < nc && ww > 0.0) {
          idx[k] = rr * nc + cc;
          w[k] = ww;
          k++;
        }
      }
    }
    return k;
  }

  static double sat_vp(double T) {               // kPa, Tetens (FAO-56 eq. 11)
    return 0.6108 * std::exp(17.27 * T / (T + 237.3));
  }

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; i++) {

      // -------------------------------------------------------
      // Inputs
      // -------------------------------------------------------
      double awc  = awcv[i];
      double hl   = hlv[i];
      double soil = soilv[i];
      double pack = packv[i];
      const double latitude = latv[i];
      const double tan_lat = std::tan(latitude * WB_PI / 180.0);

      double Pm[12];
      double Tm[12];
      double TRm[12];   // cloud transmittance
      double EAm[12];   // actual vapour pressure, kPa
      double VPDm[12];  // vapour-pressure deficit at the cell, kPa (Penman-Monteith)
      double U2m[12];   // 2 m wind speed, m s-1 (Penman-Monteith)

      bool ok = covered[i] == 1 &&
                std::isfinite(awc) && awc > 0.0 &&
                std::isfinite(latitude) && std::abs(latitude) <= 90.0;

      // -------------------------------------------------------
      // Monthly climate at this cell
      // -------------------------------------------------------
      if (ok) {
        int idx[4]; double w[4];
        const int k = neighbours(rowf[i], colf[i], nrc, ncc, idx, w);

        for (int m = 0; m < 12; m++) {
          double sp = 0.0, wp = 0.0, st = 0.0, wt = 0.0, sv = 0.0, wv = 0.0;
          for (int q = 0; q < k; q++) {
            const double pv = Pc(idx[q], m);
            const double tv = Tc(idx[q], m);
            if (std::isfinite(pv)) { sp += w[q] * pv; wp += w[q]; }
            if (std::isfinite(tv)) { st += w[q] * tv; wt += w[q]; }
            if (method == 3) {
              const double vv = Vc(idx[q], m);
              if (std::isfinite(vv)) { sv += w[q] * vv; wv += w[q]; }
            }
          }
          if (wp <= 0.0 || wt <= 0.0) { ok = false; break; }

          // Elevation adjustments
          const double t_coarse = st / wt;
          Tm[m] = t_coarse + dz[i] * lapse_t(i, m);
          Pm[m] = sp / wp * std::exp(dz[i] * beta_p(i, m));
          if (!std::isfinite(Pm[m]) || !std::isfinite(Tm[m])) { ok = false; break; }
          if (method == 3) {
            // Vapour-pressure deficit downscaled by conserving relative humidity through the
            // temperature lapse: VPD_cell = VPD_coarse x es(T_cell) / es(T_coarse) (v6.0).
            if (wv <= 0.0) { ok = false; break; }
            double vpd = sv / wv * sat_vp(Tm[m]) / sat_vp(t_coarse);
            if (!std::isfinite(vpd) || vpd < 0.0) vpd = 0.0;
            VPDm[m] = vpd;
          }
        }
      }

      // Static normals (transmittance, dew point): interpolated where the
      // static grid reaches the cell, otherwise the configured fallbacks. Needed by
      // Priestley-Taylor and Penman-Monteith PET and by the energy melt option, so
      // they are prepared for every PET method (v6.0; 5.5 skipped them under Hamon).
      if (ok) {
        int idx[4]; double w[4];
        int k = 0;
        if (covered_s[i] == 1) k = neighbours(rowf_s[i], colf_s[i], nrs, ncs, idx, w);
        for (int m = 0; m < 12; m++) {
          double ss = 0.0, ws = 0.0, se = 0.0, we = 0.0;
          for (int q = 0; q < k; q++) {
            const double sv = Sc(idx[q], m);
            const double ev = Ec(idx[q], m);
            if (std::isfinite(sv)) { ss += w[q] * sv; ws += w[q]; }
            if (std::isfinite(ev)) { se += w[q] * ev; we += w[q]; }
          }
          double tr = ws > 0.0 ? ss / ws : trans_fill[m];
          if (!std::isfinite(tr)) tr = trans_fill[m];
          TRm[m] = std::max(0.05, std::min(tr, 1.0));
          const double es = sat_vp(Tm[m]);
          double ea = we > 0.0 ? sat_vp(se / we) : rh_fill * es;
          if (!std::isfinite(ea)) ea = rh_fill * es;
          EAm[m] = std::max(0.001, std::min(ea, es));   // dew point cannot exceed air temperature
        }
      }
      // Wind climatology (Penman-Monteith): interpolated where its grid reaches the cell,
      // otherwise the configured constant (v6.0).
      if (ok && method == 3) {
        int idx[4]; double w[4];
        int k = 0;
        if (covered_w[i] == 1) k = neighbours(rowf_w[i], colf_w[i], nrw, ncw, idx, w);
        for (int m = 0; m < 12; m++) {
          double su = 0.0, wu = 0.0;
          for (int q = 0; q < k; q++) {
            const double uv = Wc(idx[q], m);
            if (std::isfinite(uv)) { su += w[q] * uv; wu += w[q]; }
          }
          double u = wu > 0.0 ? su / wu : wind_fill;
          if (!std::isfinite(u)) u = wind_fill;
          U2m[m] = std::max(0.1, std::min(u, 20.0));
        }
      }

      if (!ok) {
        for (int c = 0; c < (int) out.ncol(); c++) out(i, c) = NA_REAL;
        continue;
      }

      // -------------------------------------------------------
      // Initialize state
      // -------------------------------------------------------
      if (!std::isfinite(hl)) hl = 1.0;
      if (hl < 0.0) hl = 0.0;
      if (!std::isfinite(soil)) soil = awc;
      if (!std::isfinite(pack)) pack = 0.0;
      soil = std::max(0.0, std::min(soil, awc));
      pack = std::max(0.0, pack);

      double def_sum = 0.0;
      double aet_sum = 0.0;
      double pet_sum = 0.0;

      // Atmospheric pressure and psychrometric constant from elevation (FAO-56 eq. 7-8).
      double z = elevv[i];
      if (!std::isfinite(z)) z = 0.0;
      const double Patm = 101.3 * std::pow((293.0 - 0.0065 * z) / 293.0, 5.26);
      const double gam = 0.000665 * Patm;

      // -------------------------------------------------------
      // Monthly water balance
      // -------------------------------------------------------
      for (int m = 0; m < 12; m++) {

        double P = Pm[m];
        double T = Tm[m];
        if (P < 0.0) P = 0.0;
        if (T < -80.0) T = -80.0;

        const double nd = month_days[m];
        // Shortwave on this cell (clear-sky on its slope x cloud transmittance) and FAO-56
        // eq. 39 net longwave loss, MJ m-2 d-1; used by Priestley-Taylor PET, the energy
        // melt option and the monthly radiation outputs (v5.5).
        const double Rs  = TRm[m] * radv(i, m);
        const double Tk  = T + 273.16;
        const double ratio = std::max(0.3, std::min(TRm[m], 1.0));
        const double Rnl = SIGMA * Tk * Tk * Tk * Tk *
                           (0.34 - 0.14 * std::sqrt(EAm[m])) * (1.35 * ratio - 0.35);
        const double Rn_snow = (1.0 - snow_albedo) * Rs - Rnl;   // net radiation at a snow surface
        double PET = 0.0;

        if (method == 1) {
          // Hamon (as v4), scaled by the heat-load modifier.
          const double cos_omega = std::max(-1.0, std::min(1.0, -tan_lat * std::tan(declination[m])));
          const double dl = 24.0 / WB_PI * std::acos(cos_omega);
          const double Wt = 4.95 * std::exp(0.062 * T) / 100.0;
          const double D  = dl / 12.0;
          PET = 13.97 * nd * D * D * Wt;
          if (!std::isfinite(PET) || PET < 0.0) PET = 0.0;
          PET *= hl;
        } else if (method == 2) {
          // Priestley-Taylor. Shortwave: clear-sky radiation on this cell's
          // slope times the month's cloud transmittance. Longwave: FAO-56
          // eq. 39 with Rs/Rso = transmittance and mean monthly temperature.
          // Soil heat flux: FAO-56 eq. 43 for monthly periods.
          const double Rn  = (1.0 - albedo) * Rs - Rnl;
          const double G   = 0.07 * (Tm[(m + 1) % 12] - Tm[(m + 11) % 12]);
          const double es  = sat_vp(T);
          const double delta = 4098.0 * es / ((T + 237.3) * (T + 237.3));
          const double pet_d = pt_alpha * delta / (delta + gam) * (Rn - G) / LAMBDA;
          PET = pet_d * nd;
          if (!std::isfinite(PET) || PET < 0.0) PET = 0.0;
        } else {
          // Penman-Monteith reference demand from monthly means (v6.0; FAO-56 eq. 6 with the
          // reference-surface constants 900 and 0.34). Radiation term exactly as Priestley-Taylor
          // uses it but with alpha 1 and the wind-dependent denominator; aerodynamic term from the
          // downscaled PRISM vapour-pressure deficit and the 2 m wind climatology.
          const double Rn  = (1.0 - albedo) * Rs - Rnl;
          const double G   = 0.07 * (Tm[(m + 1) % 12] - Tm[(m + 11) % 12]);
          const double es  = sat_vp(T);
          const double delta = 4098.0 * es / ((T + 237.3) * (T + 237.3));
          const double u2  = U2m[m];
          const double rad  = delta * (Rn - G) / LAMBDA;
          const double aero = aero_scale * gam * 900.0 / (T + 273.0) * u2 * VPDm[m];
          const double pet_d = (rad + aero) / (delta + gam * (1.0 + 0.34 * u2));
          PET = pet_d * nd;
          if (!std::isfinite(PET) || PET < 0.0) PET = 0.0;
        }

        pet_sum += PET;
        if (diag == 2) out(i, 4 + m) = PET;

        // -----------------------------------------------------
        // Rain / snow
        // -----------------------------------------------------
        // Linear transition between all snow (T <= snow_t_all) and all rain (T >= rain_t_all).
        // A wide transition stands in for the mix of cold and warm days within a month.
        double rain_frac = (T - snow_t_all) / (rain_t_all - snow_t_all);
        rain_frac = std::max(0.0, std::min(rain_frac, 1.0));
        const double RAIN = P * rain_frac;
        const double SNOW = P - RAIN;

        // -----------------------------------------------------
        // Snow melt
        // -----------------------------------------------------
        double melt_potential;
        if (melt_method == 2) {
          // Energy (restricted degree-day) melt, v5.4: a temperature term standing for the
          // turbulent exchange plus a radiation term from the net radiation the snow
          // receives (clear-sky on this slope x cloud transmittance, less what snow reflects,
          // less net longwave loss). Both use monthly means, so a cold month can offset
          // strong sun (no melt) and a mild overcast winter month melts little, while a
          // clear May month melts fast. melt_rad_factor = 3.0 mm per MJ converts all net
          // radiation to melt (latent heat of fusion 0.334 MJ kg-1).
          // Net radiation at the snow surface: absorbed shortwave less net longwave loss
          // (FAO-56 eq. 39 with the month's cloudiness and vapour pressure), MJ m-2 d-1.
          melt_potential = std::max(0.0, melt_factor * T + melt_rad_factor * Rn_snow) * nd;
        } else {
        // Degree-day melt (v5.3). A share melt_rad_weight of the factor scales with the cell's
        // clear-sky radiation relative to level ground at the same latitude and elevation
        // (FAO-56 eq. 21 and 37), so north slopes melt slower and south slopes faster while
        // level ground keeps exactly the degree-day rate.
        double rad_scale = 1.0;
        if (melt_rad_weight > 0.0) {
          const double dr = 1.0 + 0.033 * std::cos(2.0 * WB_PI * MID_DOY[m] / 365.0);
          const double phi = latitude * WB_PI / 180.0;
          const double ws = std::acos(std::max(-1.0, std::min(1.0, -tan_lat * std::tan(declination[m]))));
          const double ra_flat = 24.0 * 60.0 / WB_PI * GSC * dr *
            (ws * std::sin(phi) * std::sin(declination[m]) + std::cos(phi) * std::cos(declination[m]) * std::sin(ws));
          const double rso_flat = (0.75 + 2e-5 * std::max(0.0, std::min(z, 5000.0))) * ra_flat;
          double ratio = rso_flat > 0.0 ? radv(i, m) / rso_flat : 1.0;
          if (!std::isfinite(ratio)) ratio = 1.0;
          ratio = std::max(0.2, std::min(ratio, 2.5));
          rad_scale = 1.0 - melt_rad_weight + melt_rad_weight * ratio;
        }
        melt_potential = std::max(T, 0.0) * melt_factor * rad_scale * nd;
        }
        double MELT = std::min(pack + SNOW, melt_potential);
        if (!std::isfinite(MELT) || MELT < 0.0) MELT = 0.0;
        pack = pack + SNOW - MELT;
        if (!std::isfinite(pack) || pack < 0.0) pack = 0.0;

        const double W = RAIN + MELT;

        // -----------------------------------------------------
        // Wet / dry balance
        // -----------------------------------------------------
        double AET = 0.0;
        double DEF = 0.0;
        double RUNOFF = 0.0;

        if ((W - PET) >= 0.0) {
          // Wet month
          AET = PET;
          soil += W - PET;
          RUNOFF = std::max(0.0, soil - awc);
          soil = std::max(0.0, std::min(soil, awc));
        } else {
          // Dry month
          double dSOIL = soil * (1.0 - std::exp((W - PET) / awc));
          if (!std::isfinite(dSOIL) || dSOIL < 0.0) dSOIL = 0.0;
          dSOIL = std::min(dSOIL, soil);

          AET = W + dSOIL;
          AET = std::max(0.0, std::min(AET, PET));

          double soil_withdraw = AET - W;
          if (!std::isfinite(soil_withdraw) || soil_withdraw < 0.0) soil_withdraw = 0.0;
          soil_withdraw = std::min(soil_withdraw, soil);
          soil -= soil_withdraw;
          soil = std::max(0.0, std::min(soil, awc));

          DEF = PET - AET;
          if (!std::isfinite(DEF) || DEF < 0.0) DEF = 0.0;
        }

        // -----------------------------------------------------
        // Accumulate requested months
        // -----------------------------------------------------
        if (keep[m]) {
          def_sum += DEF;
          aet_sum += AET;
        }

        if (monthly) {
          int col = mcol0;
          if (monthly & 1)  { out(i, col + m) = DEF;  col += 12; }
          if (monthly & 2)  { out(i, col + m) = AET;  col += 12; }
          if (monthly & 4)  { out(i, col + m) = PET;  col += 12; }
          if (monthly & 8)  { out(i, col + m) = W;    col += 12; }
          if (monthly & 16) { out(i, col + m) = soil; col += 12; }
          if (monthly & 32) { out(i, col + m) = pack; col += 12; }
          if (monthly & 64) { out(i, col + m) = P;    col += 12; }
          if (monthly & 128) { out(i, col + m) = SNOW; col += 12; }
          if (monthly & 256) { out(i, col + m) = T;    col += 12; }   // downscaled mean temperature, C (v5.5)
          if (monthly & 512) { out(i, col + m) = Rs;   col += 12; }   // incoming shortwave on the slope, MJ m-2 d-1
          if (monthly & 1024) { out(i, col + m) = Rn_snow; col += 12; } // net radiation at a snow surface, MJ m-2 d-1
          if (monthly & 2048) { out(i, col + m) = MELT; col += 12; }  // snowmelt, mm
          if (monthly & 4096) { out(i, col + m) = method == 3 ? VPDm[m] : NA_REAL; col += 12; }
          if (monthly & 8192) { out(i, col + m) = RUNOFF; col += 12; } // excess water/runoff (mm)
        }
      }

      // -------------------------------------------------------
      // Output
      // -------------------------------------------------------
      out(i, 0) = def_sum;
      out(i, 1) = aet_sum;
      out(i, 2) = soil;
      out(i, 3) = pack;
      if (diag == 1) out(i, 4) = pet_sum;
    }
  }
};


// [[Rcpp::export]]
NumericMatrix wb_year_cpp(
    NumericMatrix Pc,
    NumericMatrix Tc,
    NumericMatrix Vc,
    int nrc,
    int ncc,
    NumericVector rowf,
    NumericVector colf,
    IntegerVector covered,
    NumericVector dz,
    NumericMatrix lapse_t,
    NumericMatrix beta_p,
    NumericVector awc,
    NumericVector hl,
    NumericVector soil,
    NumericVector pack,
    NumericVector latitude,
    NumericVector elevation,
    NumericMatrix rad,
    NumericMatrix Sc,
    NumericMatrix Ec,
    int nrs,
    int ncs,
    NumericVector rowf_s,
    NumericVector colf_s,
    IntegerVector covered_s,
    NumericVector trans_fill,
    double rh_fill,
    NumericMatrix Wc,
    int nrw,
    int ncw,
    NumericVector rowf_w,
    NumericVector colf_w,
    IntegerVector covered_w,
    double wind_fill,
    double aero_scale,
    IntegerVector keep,
    NumericVector declination,
    NumericVector month_days,
    int method,
    double pt_alpha,
    double albedo,
    int diag,
    double snow_t_all,
    double rain_t_all,
    double melt_factor,
    double melt_rad_weight,
    int melt_method,
    double melt_rad_factor,
    double snow_albedo,
    int monthly
) {

  const R_xlen_t n = awc.size();

  if (Pc.ncol() != 12 || Tc.ncol() != 12 || nrc < 1 || ncc < 1 ||
      Pc.nrow() != nrc * ncc || Tc.nrow() != nrc * ncc)
    stop("Coarse climate must hold 12 months for every cell of the climate grid.");
  if (rowf.size() != n || colf.size() != n || covered.size() != n || dz.size() != n ||
      hl.size() != n || soil.size() != n || pack.size() != n || latitude.size() != n ||
      elevation.size() != n || rad.nrow() != n || rad.ncol() != 12 ||
      rowf_s.size() != n || colf_s.size() != n || covered_s.size() != n ||
      lapse_t.nrow() != n || beta_p.nrow() != n || lapse_t.ncol() != 12 || beta_p.ncol() != 12)
    stop("Per-cell inputs must all describe the same cells.");
  if (Sc.ncol() != 12 || Ec.ncol() != 12 || nrs < 1 || ncs < 1 ||
      Sc.nrow() != nrs * ncs || Ec.nrow() != nrs * ncs)
    stop("Static normals must hold 12 months for every cell of their grid.");
  if (method == 3 && (Vc.ncol() != 12 || Vc.nrow() != nrc * ncc))
    stop("Penman-Monteith needs monthly vapour-pressure deficit on the climate grid.");
  if (Wc.ncol() != 12 || nrw < 1 || ncw < 1 || Wc.nrow() != nrw * ncw ||
      rowf_w.size() != n || colf_w.size() != n || covered_w.size() != n)
    stop("Wind climatology must hold 12 months for every cell of its grid.");
  if (!(wind_fill > 0.0 && wind_fill <= 20.0) || !(aero_scale > 0.0 && aero_scale <= 2.0))
    stop("Invalid wind or aerodynamic settings.");
  if (keep.size() != 12 || declination.size() != 12 || month_days.size() != 12 || trans_fill.size() != 12)
    stop("Expected 12 monthly settings.");
  for (int m = 0; m < 12; ++m) {
    if (!std::isfinite(declination[m]) || std::abs(declination[m]) > 1.5707963267948966 ||
        !std::isfinite(month_days[m]) || month_days[m] <= 0 || (keep[m] != 0 && keep[m] != 1) ||
        !std::isfinite(trans_fill[m]) || trans_fill[m] <= 0.0 || trans_fill[m] > 1.0)
      stop("Invalid monthly settings.");
  }
  if (method != 1 && method != 2 && method != 3) stop("method must be 1 (Hamon), 2 (Priestley-Taylor) or 3 (Penman-Monteith).");
  if (!(pt_alpha > 0.0 && pt_alpha < 3.0) || !(albedo >= 0.0 && albedo < 1.0) ||
      !(rh_fill > 0.0 && rh_fill <= 1.0))
    stop("Invalid PET settings.");
  if (!std::isfinite(snow_t_all) || !std::isfinite(rain_t_all) || rain_t_all <= snow_t_all ||
      !(melt_factor > 0.0) || !(melt_rad_weight >= 0.0 && melt_rad_weight <= 1.0))
    stop("Invalid snow settings.");
  if ((melt_method != 1 && melt_method != 2) || !(melt_rad_factor >= 0.0) ||
      !(snow_albedo >= 0.0 && snow_albedo < 1.0))
    stop("Invalid snowmelt settings (melt_method 1 or 2, melt_rad_factor >= 0, snow_albedo in [0, 1)).");

  if (diag < 0 || diag > 2) stop("diag must be 0, 1 or 2.");
  if (monthly < 0 || monthly > 16383) stop("monthly must be a bitmask in 0..16383.");
  const int mcol0 = diag == 2 ? 16 : (diag == 1 ? 5 : 4);
  int nmonthly = 0;
  for (int bit = 0; bit < 14; bit++) if (monthly & (1 << bit)) nmonthly++;
  const int ncol_out = mcol0 + 12 * nmonthly;
  NumericMatrix out(n, ncol_out);

  WBYearWorker worker(
      Pc, Tc, Vc, nrc, ncc, rowf, colf, covered, dz, lapse_t, beta_p,
      awc, hl, soil, pack, latitude, elevation, rad, Sc, Ec, nrs, ncs,
      rowf_s, colf_s, covered_s, trans_fill, rh_fill,
      Wc, nrw, ncw, rowf_w, colf_w, covered_w, wind_fill, aero_scale, keep, declination, month_days,
      method, pt_alpha, albedo, diag, snow_t_all, rain_t_all, melt_factor, melt_rad_weight,
      melt_method, melt_rad_factor, snow_albedo, monthly, mcol0, out
  );

  parallelFor(0, n, worker);

  CharacterVector nm(ncol_out);
  nm[0] = "cwd"; nm[1] = "aet"; nm[2] = "soil_end"; nm[3] = "pack_end";
  if (diag == 1) nm[4] = "pet";
  if (diag == 2) {
    for (int m = 0; m < 12; m++) {
      char buf[8];
      std::snprintf(buf, sizeof(buf), "pet_%02d", m + 1);
      nm[4 + m] = buf;
    }
  }
  if (monthly) {
    const char* labels[14] = {"cwd", "aet", "pet", "water", "soil", "pack", "ppt", "snowfall", "tmean", "rs", "rnet", "melt", "vpd", "runoff"};
    int col = mcol0;
    for (int bit = 0; bit < 14; bit++) {
      if (!(monthly & (1 << bit))) continue;
      for (int m = 0; m < 12; m++) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%s_m%02d", labels[bit], m + 1);
        nm[col + m] = buf;
      }
      col += 12;
    }
  }
  colnames(out) = nm;
  return out;
}
)---"

wb_compile <- function(cfg) {
  # Cached shared libraries still need the RcppParallel/TBB dependency loaded.
  loadNamespace("RcppParallel")
  dir.create(cfg$rcpp_cache_dir, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(code = wb_cpp_source, env = environment(wb_compile),
                  cacheDir = cfg$rcpp_cache_dir, showOutput = FALSE)
  RcppParallel::setThreadOptions(numThreads = cfg$wb_cores)
}


wb_fit_lapse_from_normals <- function(tmean_stack, ppt_stack, z_lo_m) {
  
  lapse_t <- rep(NA_real_, 12)
  beta_p  <- rep(NA_real_, 12)
  
  for (m in 1:12) {
    
    cat("Fitting lapse month ", m, "...\n", sep = "")
    
    t_lo <- values(tmean_stack[[m]])
    ok_t <- is.finite(t_lo) & is.finite(z_lo_m)
    idx <- which(ok_t)
    
    if (length(idx) < 10) stop("Too few valid PRISM cells for tmean month ", m)
    if (length(idx) > 200000) idx <- sample(idx, 200000)
    
    fit_t <- lm(t_lo[idx] ~ z_lo_m[idx])
    lapse_t[m] <- coef(fit_t)[2]
    
    if (!is.finite(lapse_t[m]) || lapse_t[m] > -0.001 || lapse_t[m] < -0.012) {
      warning(sprintf("Suspicious tmean lapse month %d: %.5f C/m", m, lapse_t[m]))
    }
    
    p_lo <- values(ppt_stack[[m]])
    ok_p <- is.finite(p_lo) & is.finite(z_lo_m) & p_lo > 0
    idx <- which(ok_p)
    
    if (length(idx) < 10) stop("Too few valid PRISM cells for ppt month ", m)
    if (length(idx) > 200000) idx <- sample(idx, 200000)
    
    fit_p <- lm(log(p_lo[idx]) ~ z_lo_m[idx])
    beta_p[m] <- coef(fit_p)[2]
    
    if (!is.finite(beta_p[m]) || beta_p[m] < -0.0005 || beta_p[m] > 0.005) {
      warning(sprintf("Suspicious ppt log-elev beta month %d: %.6f per m", m, beta_p[m]))
    }
  }
  
  list(lapse_t = lapse_t, beta_p = beta_p)
}

# Regional elevation slopes. Within a square moving window (cfg$lapse_window_km wide) stepped across the PRISM normals in ~10-km blocks, regress each month's normal on elevation: temperature (C per m) and log precipitation (per m).
# A window with little relief cannot identify a slope, so every fit is pulled smoothly toward the single West-wide slope (`west`); that prior carries the weight of 10% of a full window whose elevations have a 100-m standard deviation. Slopes are finally held within cfg$lapse_t_limits/lapse_p_limits.
wb_fit_lapse_local <- function(tmean_stack, ppt_stack, dem_low, west, cfg) {
  cell_km <- if (is.lonlat(dem_low)) res(dem_low)[2] * 111.32 else res(dem_low)[2] / 1000
  fact <- max(1L, as.integer(round(10 / cell_km)))
  nwin <- max(3L, as.integer(round(cfg$lapse_window_km / (fact * cell_km))))
  if (nwin %% 2L == 0L) nwin <- nwin + 1L
  w <- matrix(1, nwin, nwin)
  window_sum <- function(x) focal(aggregate(x, fact = fact, fun = "sum", na.rm = TRUE),
                                  w = w, fun = "sum", na.rm = TRUE)
  prior <- 0.1 * (nwin * fact)^2 * 100^2
  # Centered elevations keep the sums of squares numerically stable.
  z0 <- dem_low - global(dem_low, "mean", na.rm = TRUE)[1, 1]
  fit <- function(y, b0, limits) {
    ok <- is.finite(y) & is.finite(z0)
    z <- mask(z0, ok, maskvalues = 0)
    y <- mask(y, ok, maskvalues = 0)
    n <- window_sum(ok); sz <- window_sum(z); sy <- window_sum(y)
    sxx <- window_sum(z * z) - sz * sz / n
    sxy <- window_sum(z * y) - sz * sy / n
    b <- (sxy + prior * b0) / (sxx + prior)
    clamp(ifel(is.finite(b), b, b0), limits[1], limits[2], values = TRUE)
  }
  lapse_t <- rast(lapply(1:12, function(m) {
    cat("Fitting regional slopes month ", m, "...\n", sep = "")
    fit(tmean_stack[[m]], west$lapse_t[m], cfg$lapse_t_limits)
  }))
  beta_p <- rast(lapply(1:12, function(m) fit(log(ppt_stack[[m]]), west$beta_p[m], cfg$lapse_p_limits)))
  names(lapse_t) <- paste0("lapse_C_per_m_", sprintf("%02d", 1:12))
  names(beta_p) <- paste0("ppt_log_per_m_", sprintf("%02d", 1:12))
  # Blocks that contain study-area cells, for reporting and mapping only.
  inside <- aggregate(is.finite(dem_low), fact = fact, fun = "sum", na.rm = TRUE) > 0
  list(lapse_t = lapse_t, beta_p = beta_p, inside = inside)
}

# ---------------- In-memory yearly engine ----------------
# The yearly calculation runs on plain vectors with one value per cell that has
# soil, not on rasters. For the 11 western states at 400 m that is about
# 19.5 million cells and roughly 8 GB of RAM. Rasters are used only for setup
# and for writing results.

# Fixed inputs for every cell with soil, including the 12-month clear-sky
# radiation table computed from each cell's slope, aspect, latitude and
# elevation (see RadWorker in the C++ source). Only the sloped-surface table is
# kept for all cells; the terrain modifier (slope / level ground) is summarized
# from a random subsample and mapped separately, block by block
# (wb_write_terrain_modifier), so no second full-size table is ever held.
wb_memory_site <- function(site, cfg, keep_cells = NULL) {
  awc <- values(site$SOILmax, mat = FALSE)
  cells <- which(is.finite(awc) & awc > 0)
  # Points-only runs (WRITE_RASTERS FALSE) simulate just the cells the points need (v5.5).
  if (!is.null(keep_cells)) cells <- cells[cells %in% keep_cells]
  if (!length(cells)) stop(if (is.null(keep_cells)) "No cells with positive soil water capacity." else
                           "None of the points falls on a modelled cell (inside the AOI, with soil).")
  wb_check_memory(length(cells), cfg)
  slopes <- function(x) if (inherits(x, "SpatRaster")) values(x)[cells, , drop = FALSE] else
    matrix(rep(as.numeric(x), each = length(cells)), ncol = 12)
  pick <- function(x) values(x, mat = FALSE)[cells]
  mem <- list(cells = cells, n = length(cells), template = rast(site$dem),
       xy = xyFromCell(site$dem, cells), awc = awc[cells],
       hl = pick(site$HL), latitude = pick(site$latitude), elevation = pick(site$dem),
       lapse_t = slopes(site$lapse_t), beta_p = slopes(site$beta_p))
  rm(awc); gc(verbose = FALSE)
  old <- wb_enter("Computing clear-sky radiation by month"); on.exit(.wb_log$context <- old, add = TRUE)
  cat("Computing monthly clear-sky radiation for ", format(mem$n, big.mark = ","), " cells...\n", sep = "")
  slope <- pick(site$slope_rad); aspect <- pick(site$aspect_deg)
  mem$rad <- wb_radiation_cpp(slope, aspect, mem$latitude, mem$elevation,
                              site$solar_declination, cfg$rad_diffuse_fraction, cfg$surface_albedo,
                              cfg$rad_step_minutes, FALSE)$slope
  # Terrain modifier statistics from a reproducible subsample.
  set.seed(cfg$lapse_seed)
  sub <- if (mem$n > 500000L) sort(sample.int(mem$n, 500000L)) else seq_len(mem$n)
  r <- wb_radiation_cpp(slope[sub], aspect[sub], mem$latitude[sub], mem$elevation[sub],
                        site$solar_declination, cfg$rad_diffuse_fraction, cfg$surface_albedo,
                        cfg$rad_step_minutes, TRUE)
  ratio <- r$slope / r$flat
  mem$rad_summary <- data.frame(month = 1:12,
    clear_sky_level_MJ_m2_d = colMeans(r$flat, na.rm = TRUE),
    clear_sky_slope_MJ_m2_d = colMeans(mem$rad, na.rm = TRUE),
    terrain_modifier_p05 = apply(ratio, 2, quantile, probs = .05, na.rm = TRUE),
    terrain_modifier_median = apply(ratio, 2, median, na.rm = TRUE),
    terrain_modifier_p95 = apply(ratio, 2, quantile, probs = .95, na.rm = TRUE),
    cells_sampled = length(sub))
  rm(r, ratio, slope, aspect); gc(verbose = FALSE)
  mem
}

# Rough size of the in-memory engine. Stops before any large allocation when
# the estimate exceeds the configured limit (or physical memory on macOS).
wb_check_memory <- function(n, cfg) {
  bytes_per_cell <- 8 * (11 + 3 * 12 + 2 + 5 + 4 + 3) + 8 +  # vectors, three 12-column tables, outputs, state, wind-grid positions (v6.0)
    8 * 12 * min(length(cfg$monthly_output), cfg$monthly_chunk)   # monthly output columns held at once
  need_gb <- n * bytes_per_cell / 1024^3
  limit <- cfg$memory_limit_gb
  if (!is.finite(limit) && Sys.info()[["sysname"]] == "Darwin") {
    phys <- suppressWarnings(as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE)))
    if (is.finite(phys)) limit <- 0.8 * phys / 1024^3
  }
  msg <- sprintf("In-memory engine: %s cells, about %.1f GB (plus raster working space).",
                 format(n, big.mark = ","), need_gb)
  cat(msg, "\n"); wb_note(msg)
  if (is.finite(limit) && need_gb > limit) {
    fact <- ceiling(sqrt(need_gb / limit))
    stop(sprintf("%s That exceeds the %.1f GB available for this run. Increase DEM_AGGREGATE_FACT (currently %d; at least %d times coarser would fit), use a smaller AOI, or raise MEMORY_LIMIT_GB if you know the machine can take it.",
                 msg, limit, as.integer(cfg$dem_aggregate_fact), as.integer(cfg$dem_aggregate_fact * fact)))
  }
  invisible(need_gb)
}

# 12-layer terrain modifier (clear-sky radiation on the slope / on level ground
# at the same place), written block by block as Float32 so the full table never
# exists in memory.
wb_write_terrain_modifier <- function(site, cfg, path) {
  old <- wb_enter("Writing terrain modifier map", files = path); on.exit(.wb_log$context <- old, add = TRUE)
  out <- rast(site$dem, nlyrs = 12, names = paste0("terrain_modifier_", sprintf("%02d", 1:12)))
  tmp <- tempfile(pattern = "partial_", tmpdir = dirname(path), fileext = ".tif")
  on.exit(unlink(c(tmp, paste0(tmp, ".aux.xml"))), add = TRUE)
  src <- c(site$dem, site$slope_rad, site$aspect_deg, site$latitude)
  readStart(src); on.exit(readStop(src), add = TRUE)
  b <- writeStart(out, tmp, overwrite = TRUE, wopt = list(datatype = "FLT4S", gdal = cfg$gdal_opts))
  for (i in seq_len(b$n)) {
    v <- readValues(src, row = b$row[i], nrows = b$nrows[i], mat = TRUE)
    res <- matrix(NA_real_, nrow(v), 12)
    ok <- which(is.finite(v[, 1]))
    if (length(ok)) {
      r <- wb_radiation_cpp(v[ok, 2], v[ok, 3], v[ok, 4], v[ok, 1], site$solar_declination,
                            cfg$rad_diffuse_fraction, cfg$surface_albedo, cfg$rad_step_minutes, TRUE)
      res[ok, ] <- r$slope / r$flat
    }
    writeValues(out, res, b$row[i], b$nrows[i])
  }
  out <- writeStop(out)
  if (!file.rename(tmp, path)) stop("Could not atomically replace ", path)
  invisible(rast(path))
}

# Static monthly normals (cloud transmittance, dew point) located on their own
# coarse grid, or constant fallbacks when a normal is not used.
wb_memory_static <- function(mem, site, cfg) {
  none <- list(Sc = matrix(NA_real_, 1, 12), Ec = matrix(NA_real_, 1, 12), nr = 1L, nc = 1L,
               rowf = numeric(mem$n), colf = numeric(mem$n), covered = integer(mem$n))
  trans_fill <- if (is.numeric(cfg$cloud_transmittance)) rep(cfg$cloud_transmittance, 12) else rep(0.7, 12)
  rh_fill <- if (is.numeric(cfg$humidity_source)) cfg$humidity_source else 0.5
  # Wind climatology on its own grid (v6.0), or the constant everywhere.
  wind <- list(Wc = matrix(NA_real_, 1, 12), nrw = 1L, ncw = 1L, rowf_w = numeric(mem$n), colf_w = numeric(mem$n),
               covered_w = integer(mem$n), wind_fill = cfg$wind_constant)
  if (!is.null(site$wind)) {
    gw <- wb_memory_grid(mem, site, site$wind[[1]], site$dz)
    wind <- list(Wc = values(site$wind), nrw = gw$nr, ncw = gw$nc, rowf_w = gw$rowf, colf_w = gw$colf,
                 covered_w = gw$covered, wind_fill = cfg$wind_constant)
    cat("Wind climatology reaches ", format(sum(gw$covered), big.mark = ","), " of ", format(mem$n, big.mark = ","),
        " cells; the rest use ", cfg$wind_constant, " m/s.\n", sep = "")
  }
  layers <- site$static
  if (is.null(layers)) return(c(none, list(trans_fill = trans_fill, rh_fill = rh_fill), wind))
  ref <- layers[[1]]
  grid <- wb_memory_grid(mem, site, rast(lapply(layers, function(x) x[[1]])), site$dz)
  if (!is.null(layers$soltrans)) {
    trans_fill <- as.numeric(global(layers$soltrans, "mean", na.rm = TRUE)[, 1])
    trans_fill[!is.finite(trans_fill)] <- 0.7
    trans_fill <- pmin(pmax(trans_fill, 0.05), 1)
  }
  c(list(Sc = if (is.null(layers$soltrans)) matrix(NA_real_, ncell(ref), 12) else values(layers$soltrans),
         Ec = if (is.null(layers$tdmean)) matrix(NA_real_, ncell(ref), 12) else values(layers$tdmean),
         nr = grid$nr, nc = grid$nc, rowf = grid$rowf, colf = grid$colf, covered = grid$covered,
         trans_fill = trans_fill, rh_fill = rh_fill), wind)
}

# The per-cell parts of the static inputs for a subset of engine rows (points, PET check).
wb_static_subset <- function(s, idx) {
  list(Sc = s$Sc, Ec = s$Ec, nr = s$nr, nc = s$nc, rowf = s$rowf[idx], colf = s$colf[idx], covered = s$covered[idx],
       trans_fill = s$trans_fill, rh_fill = s$rh_fill,
       Wc = s$Wc, nrw = s$nrw, ncw = s$ncw, rowf_w = s$rowf_w[idx], colf_w = s$colf_w[idx], covered_w = s$covered_w[idx],
       wind_fill = s$wind_fill)
}

# Identifies a climate grid by its dimensions, extent and coordinate system.
# Extents are compared to 5 decimals of the grid's units: PRISM precipitation
# and temperature files differ in the seventh decimal of a degree.
wb_grid_key <- function(x) paste(c(nrow(x), ncol(x), round(as.vector(ext(x)), 5), crs(x, proj = TRUE)), collapse = " ")

# Where every cell sits on one coarse climate grid, and its elevation
# difference (dz) from that grid. `layers` holds one or more data layers of the
# grid. A cell is covered where terra's bilinear interpolation of every one of
# those layers reaches it, so the same cells are simulated as when each month
# was interpolated with project().
wb_memory_grid <- function(mem, site, layers, dz) {
  old <- wb_enter("Locating cells on the climate grid", files = sources(layers)); on.exit(.wb_log$context <- old, add = TRUE)
  pos <- terra::project(mem$xy, crs(site$dem), crs(layers))
  reach <- values(project(layers, site$dem, method = "bilinear"))[mem$cells, , drop = FALSE]
  list(key = wb_grid_key(layers), nr = nrow(layers), nc = ncol(layers),
       rowf = (ymax(layers) - pos[, 2]) / yres(layers) - 0.5,
       colf = (pos[, 1] - xmin(layers)) / xres(layers) - 0.5,
       covered = as.integer(rowSums(is.finite(reach)) == ncol(reach)),
       dz = values(dz, mat = FALSE)[mem$cells])
}

# One year for every cell: interpolate the 12 monthly grids, apply the
# elevation adjustments and run the monthly water balance. Returns columns
# cwd, aet, soil_end and pack_end; cells that cannot be simulated are NA.
wb_year <- function(mem, grid, ppt, tmean, soil, pack, site, awc = mem$awc,
                    method = site$pet_method, diag = 0L, monthly = 0L, vpd = NULL) {
  old <- wb_enter("Monthly water-balance calculation"); on.exit(.wb_log$context <- old, add = TRUE)
  if (nlyr(ppt) != 12 || nlyr(tmean) != 12 || wb_grid_key(ppt) != grid$key ||
      !compareGeom(ppt, tmean, stopOnError = FALSE))
    stop("Monthly climate does not match the prepared climate grid.")
  method_code <- switch(method, hamon = 1L, priestley_taylor = 2L, penman_monteith = 3L,
                        stop("Unknown PET method: ", method))
  # Vapour-pressure deficit (v6.0): on the climate grid when Penman-Monteith runs; a 1 x 12 dummy otherwise.
  Vc <- matrix(NA_real_, 1, 12)
  if (method_code == 3L) {
    if (is.null(vpd)) stop("Penman-Monteith needs the monthly vapour-pressure deficit.")
    if (nlyr(vpd) != 12 || !compareGeom(ppt, vpd, stopOnError = FALSE))
      stop("Monthly vapour-pressure deficit does not match the climate grid.")
    Vc <- values(vpd)
  }
  s <- mem$static
  wb_year_cpp(values(ppt), values(tmean), Vc, grid$nr, grid$nc, grid$rowf, grid$colf,
              grid$covered, grid$dz, mem$lapse_t, mem$beta_p, awc, mem$hl, soil, pack,
              mem$latitude, mem$elevation, mem$rad, s$Sc, s$Ec, s$nr, s$nc, s$rowf, s$colf,
              s$covered, s$trans_fill, s$rh_fill,
              s$Wc, s$nrw, s$ncw, s$rowf_w, s$colf_w, s$covered_w, s$wind_fill, site$aero_scale,
              site$keep, site$solar_declination,
              site$month_days, method_code, site$pt_alpha,
              site$surface_albedo, as.integer(diag), site$snow_t_all, site$rain_t_all,
              site$melt_factor, site$melt_rad_weight, site$melt_method, site$melt_rad_factor,
              site$snow_albedo, as.integer(monthly))
}

# Per-cell values back onto the model grid (one layer per column). Layers are
# filled one at a time so the peak is one full-grid vector, not a matrix.
wb_to_raster <- function(mem, x, layer_names) {
  x <- as.matrix(x)
  layers <- lapply(seq_len(ncol(x)), function(j) {
    v <- rep(NA_real_, ncell(mem$template))
    v[mem$cells] <- x[, j]
    r <- rast(mem$template); values(r) <- v
    r
  })
  r <- rast(layers)
  names(r) <- layer_names
  r
}

wb_change_stats <- function(previous, current, tolerance) {
  delta <- abs(current - previous)
  ok <- !is.na(delta)
  list(mean = if (any(ok)) mean(delta[ok]) else NA_real_,
       max = if (any(ok)) max(delta[ok]) else NA_real_,
       pct = if (any(ok)) 100 * mean(delta[ok] > tolerance) else NA_real_,
       count = sum(ok), lost = sum(!is.na(previous) & is.na(current)), delta = delta)
}

# Carried soil and snow states are saved whole, then renamed into place, so an
# interrupted write never leaves a half-written restart file.
wb_save_state <- function(path, state) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(state, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Could not commit restart state: ", path)
  invisible(state)
}

wb_read_state <- function(path, key) {
  tryCatch({
    if (!file.exists(path)) return(NULL)
    state <- readRDS(path)
    if (identical(state$key, key)) state else NULL
  }, error = function(e) NULL)
}

# Crop a national climate raster to the model area plus a four-cell margin.
wb_crop_climate <- function(x, site) {
  e <- as.vector(ext(project(as.polygons(ext(site$dem), crs = crs(site$dem)), x)))
  crop(x, ext(e + c(-4, 4, -4, 4) * rep(res(x), each = 2)), snap = "out")
}


# ---------------- Configuration, runtime, and output integrity ----------------
wb_validate <- function(cfg) {
  for (field in c("base_dir", "dem_path", "aoi_path", "awc_path", "prism_root", "output_root"))
    cfg[[field]] <- path.expand(cfg[[field]])
  # v5.4: AOI_BBOX (xmin, xmax, ymin, ymax in decimal degrees) replaces the AOI shapefile when set.
  if (is.null(cfg$aoi_bbox) || identical(cfg$aoi_bbox, NA) || !length(cfg$aoi_bbox)) cfg$aoi_bbox <- NULL
  if (!is.null(cfg$aoi_bbox)) {
    b <- cfg$aoi_bbox
    if (!is.numeric(b) || length(b) != 4 || any(!is.finite(b)) || b[1] >= b[2] || b[3] >= b[4] ||
        b[1] < -180 || b[2] > 180 || b[3] < -90 || b[4] > 90)
      stop("AOI_BBOX must be c(xmin, xmax, ymin, ymax) in decimal degrees (longitude first), e.g. c(-122.5, -121.5, 44, 45.5).")
    cfg$aoi_bbox <- as.numeric(b)
  }
  for (field in c("dem_path", if (is.null(cfg$aoi_bbox)) "aoi_path", "awc_path"))
    if (!file.exists(cfg[[field]])) stop("Missing ", field, ": ", cfg[[field]])
  check_years <- function(x, label) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x != as.integer(x)) ||
        any(diff(x) != 1)) stop(label, " must be increasing, consecutive integer years.")
    as.integer(x)
  }
  cfg$base_years <- check_years(cfg$base_years, "BASE_YEARS")
  cfg$years_annual <- check_years(cfg$years_annual, "YEARS_ANNUAL")
  if (length(intersect(cfg$base_years, cfg$years_annual)))
    stop("Baseline and analysis years must not overlap.")
  if (!is.null(cfg$analysis_initialization))
    stop("analysis_initialization is removed: analysis always continues from the end of the baseline. Remove this option.")
  if (min(cfg$years_annual) != max(cfg$base_years) + 1L)
    stop("Analysis must begin in the year immediately following base_years; soil and snow state cannot be reset or skipped.")
  if (!is.numeric(cfg$months_out) || !length(cfg$months_out) ||
      any(!is.finite(cfg$months_out)) || any(cfg$months_out != as.integer(cfg$months_out)) ||
      any(!cfg$months_out %in% 1:12)) stop("MONTHS_OUT must contain integers in 1:12.")
  cfg$months_out <- sort(unique(as.integer(cfg$months_out)))
  cfg$month_tag <- wb_month_tag(cfg$months_out)
  valid_annual <- c("cwd", "aet", "pet", "ppt", "soil", "pack", "snow", "snow_peak", "melt", "runoff", "snowfall", "water")
  if (!length(cfg$outputs_to_write) || any(!cfg$outputs_to_write %in% valid_annual))
    stop("OUTPUTS_TO_WRITE must contain: ", paste(valid_annual, collapse = ", "))
  if (is.null(cfg$summarize)) cfg$summarize <- FALSE
  if (is.null(cfg$write_baseline)) cfg$write_baseline <- FALSE
  for (field in c("spinup_max_iter", "wb_cores", "dem_aggregate_fact")) {
    x <- cfg[[field]]
    if (length(x) != 1 || !is.finite(x) || x < 1 || x != as.integer(x)) stop("Invalid ", field)
  }
  for (field in c("soil_tol_mm", "pack_tol_mm"))
    if (length(cfg[[field]]) != 1 || !is.finite(cfg[[field]]) || cfg[[field]] <= 0) stop("Invalid ", field)
  if (!cfg$dem_z_units %in% c("auto", "meters", "feet")) stop("Invalid DEM_Z_UNITS.")
  if (!cfg$awc_units %in% c("mm", "cm", "m", "fraction")) stop("Invalid AWC_UNITS.")
  if (cfg$awc_units == "fraction" && (length(cfg$awc_depth_mm) != 1 ||
      !is.finite(cfg$awc_depth_mm) || cfg$awc_depth_mm <= 0))
    stop("Volumetric AWC requires an explicit positive AWC_DEPTH_MM.")
  if (!cfg$spinup_action %in% c("warn", "stop")) stop("Invalid SPINUP_ACTION.")
  if (!is.finite(cfg$spinup_max_pct) || cfg$spinup_max_pct < 0 || cfg$spinup_max_pct > 100)
    stop("Invalid SPINUP_MAX_PCT.")
  if (length(cfg$lapse_window_km) != 1 || is.na(cfg$lapse_window_km) || cfg$lapse_window_km < 30)
    stop("LAPSE_WINDOW_KM must be at least 30, or Inf for one West-wide fit.")
  for (field in c("lapse_t_limits", "lapse_p_limits"))
    if (length(cfg[[field]]) != 2 || any(!is.finite(cfg[[field]])) || cfg[[field]][1] >= cfg[[field]][2])
      stop("Invalid ", field)
  if (!is.finite(cfg$terra_memfrac) || cfg$terra_memfrac <= 0 || cfg$terra_memfrac > .9)
    stop("TERRA_MEMFRAC must be in (0, 0.9].")
  # ---- PET and radiation (v5) ----
  if (is.null(cfg$pet_method)) cfg$pet_method <- "penman_monteith"
  if (!cfg$pet_method %in% c("penman_monteith", "priestley_taylor", "hamon"))
    stop("PET_METHOD must be 'penman_monteith', 'priestley_taylor' or 'hamon'.")
  one_in <- function(x, lo, hi, label) {
    if (length(x) != 1 || !is.finite(x) || x < lo || x > hi)
      stop(label, " must be a single number in [", lo, ", ", hi, "].")
  }
  one_in(cfg$pt_alpha, 0.5, 2.5, "PT_ALPHA")
  one_in(cfg$surface_albedo, 0, 0.9, "SURFACE_ALBEDO")
  one_in(cfg$rad_diffuse_fraction, 0, 1, "RAD_DIFFUSE_FRACTION")
  one_in(cfg$rad_step_minutes, 1, 60, "RAD_STEP_MINUTES")
  static_ok <- function(x, lo, hi, label) {
    if (identical(x, "prism")) return(invisible(TRUE))
    if (is.numeric(x)) {one_in(x, lo, hi, label); return(invisible(TRUE))}
    stop(label, " must be 'prism' or a number.")
  }
  static_ok(cfg$cloud_transmittance, 0.05, 1, "CLOUD_TRANSMITTANCE")
  static_ok(cfg$humidity_source, 0.05, 1, "HUMIDITY_SOURCE")
  # PRISM normals used as fixed monthly fields every year.
  cfg$static_types <- c(if (identical(cfg$cloud_transmittance, "prism")) "soltrans",
                        if (identical(cfg$humidity_source, "prism")) "tdmean")
  # ---- Penman-Monteith inputs (v6.0) ----
  # Monthly PRISM types needed as series and normals: ppt and tmean always; vpdmax and vpdmin for Penman-Monteith.
  cfg$climate_types <- c("ppt", "tmean", if (cfg$pet_method == "penman_monteith") c("vpdmax", "vpdmin"))
  if (is.null(cfg$wind_path) || identical(cfg$wind_path, NA) || !length(cfg$wind_path) || !nzchar(cfg$wind_path[1])) cfg$wind_path <- NULL
  if (!is.null(cfg$wind_path)) {
    cfg$wind_path <- path.expand(cfg$wind_path)
    if (!file.exists(cfg$wind_path)) stop("Missing WIND_PATH: ", cfg$wind_path)
  }
  if (is.null(cfg$wind_constant)) cfg$wind_constant <- 2
  one_in(cfg$wind_constant, 0.1, 20, "WIND_CONSTANT")
  if (is.null(cfg$aero_scale)) cfg$aero_scale <- 1
  one_in(cfg$aero_scale, 0.1, 2, "AERO_SCALE")
  if (cfg$pet_method == "penman_monteith" && is.null(cfg$wind_path))
    wb_note(sprintf("WIND_PATH is NULL: Penman-Monteith uses a constant 2 m wind of %.2f m/s everywhere.", cfg$wind_constant))
  if (!is.logical(cfg$pet_diagnostics) || length(cfg$pet_diagnostics) != 1 || is.na(cfg$pet_diagnostics))
    stop("PET_DIAGNOSTICS must be TRUE or FALSE.")
  if (is.null(cfg$zonal_path) || !length(cfg$zonal_path) || identical(cfg$zonal_path, NA)) cfg$zonal_path <- NULL
  if (!is.null(cfg$zonal_path)) {
    cfg$zonal_path <- path.expand(cfg$zonal_path)
    if (!file.exists(cfg$zonal_path)) stop("Missing ZONAL_PATH: ", cfg$zonal_path)
  }
  if (is.null(cfg$zonal_monthly)) cfg$zonal_monthly <- character()
  if (any(!cfg$zonal_monthly %in% wb_monthly_vars)) stop("Invalid zonal monthly variable")
  # ---- Snow (v5.3) ----
  for (field in c("snow_t_all", "rain_t_all", "melt_mm_per_deg_day", "melt_radiation_weight"))
    if (length(cfg[[field]]) != 1 || !is.finite(cfg[[field]])) stop("Invalid ", toupper(field))
  if (cfg$rain_t_all <= cfg$snow_t_all) stop("RAIN_T_ALL must be greater than SNOW_T_ALL.")
  if (cfg$melt_mm_per_deg_day <= 0) stop("MELT_MM_PER_DEG_DAY must be positive.")
  if (cfg$melt_radiation_weight < 0 || cfg$melt_radiation_weight > 1) stop("MELT_RADIATION_WEIGHT must be in [0, 1].")
  # ---- Snowmelt method (v5.4) ----
  if (is.null(cfg$melt_method)) cfg$melt_method <- "degree_day"
  if (!cfg$melt_method %in% c("degree_day", "energy")) stop('MELT_METHOD must be "degree_day" or "energy".')
  if (is.null(cfg$melt_mm_per_mj)) cfg$melt_mm_per_mj <- 0
  if (is.null(cfg$snow_albedo)) cfg$snow_albedo <- 0.6
  for (field in c("melt_mm_per_mj", "snow_albedo"))
    if (length(cfg[[field]]) != 1 || !is.finite(cfg[[field]])) stop("Invalid ", toupper(field))
  if (cfg$melt_mm_per_mj < 0) stop("MELT_MM_PER_MJ must be zero or positive.")
  if (cfg$snow_albedo < 0 || cfg$snow_albedo >= 1) stop("SNOW_ALBEDO must be in [0, 1).")
  if (cfg$melt_method == "energy" && cfg$melt_mm_per_mj == 0)
    warning("MELT_METHOD is \"energy\" but MELT_MM_PER_MJ is 0: melt reduces to c_T * T with no radiation term.")
  # Hamon needs no radiation normals unless the energy melt option reads them (v6.0; 5.5 dropped them under Hamon
  # and the energy melt then ran on undefined radiation). Hamon + degree_day reproduces 4.0 exactly.
  if (cfg$pet_method == "hamon" && cfg$melt_method == "degree_day") {
    cfg$static_types <- character(); cfg$climate_types <- c("ppt", "tmean")
  }
  # ---- Monthly output (v5.2) ----
  vars <- wb_monthly_vars
  if (is.null(cfg$monthly_output) || identical(cfg$monthly_output, FALSE) || !length(cfg$monthly_output)) cfg$monthly_output <- character()
  if (identical(cfg$monthly_output, TRUE)) cfg$monthly_output <- c("cwd", "aet")
  if (!is.character(cfg$monthly_output) || any(!cfg$monthly_output %in% vars))
    stop("MONTHLY_OUTPUT must be FALSE, TRUE, or a subset of: ", paste(vars, collapse = ", "))
  cfg$monthly_output <- vars[vars %in% cfg$monthly_output]
  if ("vpd" %in% cfg$monthly_output && cfg$pet_method != "penman_monteith") {
    wb_note("MONTHLY_OUTPUT includes vpd but PET_METHOD is not penman_monteith: the vapour-pressure deficit is not computed, so no vpd rasters are written.")
    cfg$monthly_output <- setdiff(cfg$monthly_output, "vpd")
  }
  cfg$monthly_mask <- sum(2^(match(cfg$monthly_output, vars) - 1L))
  if (is.null(cfg$monthly_years)) cfg$monthly_years <- "analysis"
  # Resolve MONTHLY_YEARS to one integer year vector per monthly variable.
  all_years <- c(cfg$base_years, cfg$years_annual)
  resolve_years <- function(v, label) {
    if (identical(v, "analysis")) return(cfg$years_annual)
    if (identical(v, "all")) return(all_years)
    if (is.numeric(v) && length(v) && all(is.finite(v)) && all(v == as.integer(v))) return(as.integer(v))
    stop(label, " must be 'analysis', 'all', or a vector of years.")
  }
  if (!length(cfg$monthly_output) || !isTRUE(cfg$write_rasters)) {
    # No monthly GeoTIFFs requested: per-variable year schedules are irrelevant.
    cfg$monthly_years_by_var <- list()
  } else if (is.list(cfg$monthly_years)) {
    bad <- setdiff(names(cfg$monthly_years), cfg$monthly_output)
    if (is.null(names(cfg$monthly_years)) || length(bad))
      stop("MONTHLY_YEARS as a list must be named by MONTHLY_OUTPUT variables; unknown: ", paste(bad, collapse = ", "))
    cfg$monthly_years_by_var <- lapply(setNames(cfg$monthly_output, cfg$monthly_output), function(v)
      resolve_years(if (is.null(cfg$monthly_years[[v]])) "analysis" else cfg$monthly_years[[v]], paste0("MONTHLY_YEARS$", v)))
  } else {
    yrs <- resolve_years(cfg$monthly_years, "MONTHLY_YEARS")
    cfg$monthly_years_by_var <- lapply(setNames(cfg$monthly_output, cfg$monthly_output), function(v) yrs)
  }
  if (is.null(cfg$monthly_datatype)) cfg$monthly_datatype <- "FLT4S"
  if (!cfg$monthly_datatype %in% c("FLT4S", "INT2S")) stop("MONTHLY_DATATYPE must be 'FLT4S' or 'INT2S'.")
  if (is.null(cfg$monthly_chunk)) cfg$monthly_chunk <- 4L
  if (!is.numeric(cfg$monthly_chunk) || length(cfg$monthly_chunk) != 1 || !is.finite(cfg$monthly_chunk) || cfg$monthly_chunk < 1)
    stop("MONTHLY_CHUNK must be a positive whole number of monthly variables held in memory at once.")
  cfg$monthly_chunk <- as.integer(cfg$monthly_chunk)
  if (is.null(cfg$monthly_stack)) cfg$monthly_stack <- FALSE
  if (!is.logical(cfg$monthly_stack) || length(cfg$monthly_stack) != 1 || is.na(cfg$monthly_stack))
    stop("MONTHLY_STACK must be TRUE or FALSE.")
  if (is.null(cfg$memory_limit_gb)) cfg$memory_limit_gb <- NA_real_
  if (length(cfg$memory_limit_gb) != 1 || (!is.na(cfg$memory_limit_gb) && (!is.finite(cfg$memory_limit_gb) || cfg$memory_limit_gb <= 0)))
    stop("MEMORY_LIMIT_GB must be NA (use physical memory) or a positive number of GB.")
  # ---- Points and raster switch (v5.5) ----
  if (is.null(cfg$points_path) || identical(cfg$points_path, NA) || !length(cfg$points_path) || !nzchar(cfg$points_path[1])) cfg$points_path <- NULL
  if (!is.null(cfg$points_path)) {
    cfg$points_path <- path.expand(cfg$points_path)
    if (!file.exists(cfg$points_path)) stop("Missing POINTS_PATH: ", cfg$points_path)
  }
  if (is.null(cfg$points_id_field) || identical(cfg$points_id_field, NA) || !length(cfg$points_id_field)) cfg$points_id_field <- NULL
  if (is.null(cfg$points_neighbourhood)) cfg$points_neighbourhood <- 1L
  if (!is.numeric(cfg$points_neighbourhood) || length(cfg$points_neighbourhood) != 1 || !cfg$points_neighbourhood %in% c(1, 3, 5))
    stop("POINTS_NEIGHBOURHOOD must be 1, 3 or 5 (cells per side averaged around each point).")
  cfg$points_neighbourhood <- as.integer(cfg$points_neighbourhood)
  if (is.null(cfg$write_rasters)) cfg$write_rasters <- TRUE
  if (!is.logical(cfg$write_rasters) || length(cfg$write_rasters) != 1 || is.na(cfg$write_rasters)) stop("WRITE_RASTERS must be TRUE or FALSE.")
  if (!cfg$write_rasters && is.null(cfg$points_path) && is.null(cfg$zonal_path))
    stop("No outputs requested: enable write_rasters or supply zonal_path or points_path.")
  if (!cfg$write_rasters && length(cfg$monthly_output)) {
    wb_note("WRITE_RASTERS is FALSE: MONTHLY_OUTPUT rasters are not written (points receive every monthly variable regardless).")
    cfg$monthly_output <- character(); cfg$monthly_mask <- 0L; cfg$monthly_years_by_var <- list()
  }
  if (!cfg$write_rasters && isTRUE(cfg$pet_diagnostics)) {
    wb_note("WRITE_RASTERS is FALSE: PET self-check maps are skipped.")
    cfg$pet_diagnostics <- FALSE
  }
  if (is.null(cfg$prism_mode)) cfg$prism_mode <- "full"
  if (!is.character(cfg$prism_mode) || length(cfg$prism_mode) != 1L ||
      !cfg$prism_mode %in% c("full", "subset", "local"))
    stop("prism_mode must be 'full', 'subset', or 'local'.")
  if (is.null(cfg$prism_buffer_km)) cfg$prism_buffer_km <- 25
  if (!is.numeric(cfg$prism_buffer_km) || length(cfg$prism_buffer_km) != 1L ||
      !is.finite(cfg$prism_buffer_km) || cfg$prism_buffer_km < 0)
    stop("prism_buffer_km must be a nonnegative number.")
  # Local mode never initiates network downloads, regardless of allow_downloads.
  if (cfg$prism_mode == "local") cfg$allow_downloads <- FALSE
  cfg$prism_dir <- file.path(cfg$prism_root, "PRISM_cache")
  cfg$prism_dir_v2 <- file.path(cfg$prism_dir, "api_v2_clean")
  if (cfg$prism_mode == "subset") {
    # Immutable AOI fingerprint prevents overlapping or changing study areas
    # from accidentally selecting an undersized spatial subset.
    geometry_key <- if (!is.null(cfg$aoi_bbox)) paste(cfg$aoi_bbox, collapse = ",") else
      paste(normalizePath(cfg$aoi_path), paste(file.info(cfg$aoi_path)[c("size", "mtime")], collapse = ":"))
    identifier <- paste(geometry_key, cfg$prism_res, cfg$prism_buffer_km,
                        cfg$lapse_window_km, sep = "|")
    keyfile <- tempfile(); writeLines(identifier, keyfile)
    key <- substr(unname(tools::md5sum(keyfile)), 1, 16); unlink(keyfile)
    cfg$prism_subset_dir <- file.path(cfg$prism_dir, "subsets", key)
    cfg$prism_search_dirs <- cfg$prism_subset_dir
  } else cfg$prism_search_dirs <- unique(c(cfg$prism_dir, cfg$prism_dir_v2))
  if (isTRUE(cfg$run_aois_separately))
    stop("The flat output layout requires RUN_AOIS_SEPARATELY <- FALSE to avoid filename collisions.")
  cfg$work_root <- tools::R_user_dir("water_balance_v6", "cache")
  cfg$rcpp_cache_dir <- file.path(cfg$work_root, "cpp")
  cfg
}

wb_setup_runtime <- function(cfg) {
  # Use existing paths only when the expected data files exist.
  # Never point GDAL_DATA at the PROJ directory.
  if (file.exists(file.path(cfg$proj_data_dir, "proj.db")))
    Sys.setenv(PROJ_LIB = cfg$proj_data_dir, PROJ_DATA = cfg$proj_data_dir)
  if (dir.exists(cfg$gdal_data_dir)) Sys.setenv(GDAL_DATA = cfg$gdal_data_dir)
  # Use installed datum grids by default. Opportunistic PROJ network downloads
  # can fail and produce empty geometry in some terra/prism environments.
  terra::projNetwork(isTRUE(cfg$proj_network))
  dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
  td <- tempfile(pattern = "terra_water_balance_v4_", tmpdir = tempdir())
  .wb_log$temp_dir <- td
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  wb_check_space(td)
  terraOptions(progress = 1, memfrac = cfg$terra_memfrac, tempdir = td,
               threads = cfg$wb_cores)
  Sys.setenv(GDAL_NUM_THREADS = as.character(cfg$wb_cores))
  # Fail before processing if the requested target is not projected in meters.
  ref <- rast(nrows = 1, ncols = 1, crs = cfg$metric_crs)
  p <- crs(ref, proj = TRUE)
  if (is.lonlat(ref) || !grepl("\\+units=m( |$)", p))
    stop("METRIC_CRS must describe a projected coordinate system with meter units.")
  wb_compile(cfg)
  invisible(cfg)
}

wb_hash <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 2)
  unname(tools::md5sum(f))
}

wb_file_stamp <- function(paths) {
  paths <- sort(unique(normalizePath(paths, mustWork = TRUE)))
  z <- file.info(paths)
  data.frame(path = paths, size = z$size, mtime = as.numeric(z$mtime), row.names = NULL)
}

wb_input_files <- function(cfg) {
  # Shapefile geometry/projection/attributes all contribute to run identity.
  stem <- tools::file_path_sans_ext(cfg$aoi_path)
  aoi_files <- if (!is.null(cfg$aoi_bbox)) character() else if (tolower(tools::file_ext(cfg$aoi_path)) == "shp") {
    candidates <- list.files(dirname(stem), full.names = TRUE)
    candidates[tolower(tools::file_path_sans_ext(basename(candidates))) == tolower(basename(stem))]
  } else cfg$aoi_path
  unique(c(cfg$dem_path, cfg$awc_path, aoi_files, cfg$wind_path,
           c(paste0(cfg$dem_path, ".aux.xml"), paste0(cfg$awc_path, ".aux.xml"))[
             file.exists(c(paste0(cfg$dem_path, ".aux.xml"), paste0(cfg$awc_path, ".aux.xml")))]))
}

wb_start_run <- function(cfg, index) {
  # Full model settings and code hashes plus input size/mtime prevent accidental
  # reuse across months, source data, model revisions, or initialization choices.
  # Large raster inputs are not checksummed; preserve timestamps only deliberately.
  scientific <- cfg[setdiff(names(cfg), c("resume", "diagnostics", "pet_diagnostics", "memory_limit_gb", "monthly_output", "monthly_years", "monthly_stack", "monthly_mask", "monthly_datatype", "monthly_years_by_var", "monthly_chunk", "stop_after_setup",
      "wb_cores", "terra_memfrac", "gdal_opts", "rcpp_cache_dir", "work_root", "output_root", "allow_downloads"))]
  identity <- list(settings = scientific, inputs = wb_file_stamp(wb_input_files(cfg)),
                   climate = wb_file_stamp(index$file),
                   packages = vapply(c("terra", "prism", "Rcpp", "RcppParallel"),
                                     function(p) as.character(packageVersion(p)), character(1)))
  cfg$run_id <- wb_hash(identity)
  cfg$run_dir <- file.path(cfg$work_root, paste0(cfg$month_tag, "_baseline_",
    min(cfg$base_years), "_", max(cfg$base_years), "_analysis_",
    min(cfg$years_annual), "_", max(cfg$years_annual)), substr(cfg$run_id, 1, 12))
  dir.create(cfg$run_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(identity, file.path(cfg$run_dir, "run_manifest.rds"))
  writeLines(c(paste("Run:", cfg$run_id), paste("Months:", cfg$month_tag),
    paste("Baseline:", paste(range(cfg$base_years), collapse = "-")),
    paste("Analysis:", paste(range(cfg$years_annual), collapse = "-")),
    "Analysis initialization: continuous from final baseline state",
    "Input fingerprints use path, size, and modification time; see run_manifest.rds.",
    capture.output(sessionInfo())), file.path(cfg$run_dir, "run_info.txt"))
  write.csv(index, file.path(cfg$run_dir, "prism_inputs.csv"), row.names = FALSE)
  wb_note(c(paste("Run fingerprint:", cfg$run_id),
    paste("Local restart directory:", cfg$run_dir),
    "Configuration:", capture.output(str(cfg)), "Climate inputs:", capture.output(print(index, row.names = FALSE)),
    "R session:", capture.output(sessionInfo())))
  cat("\nFinal outputs:", cfg$output_root, "\nLocal restart files:", cfg$run_dir, "\n")
  cfg
}

wb_atomic_raster <- function(x, path, cfg, datatype = "FLT4S") {
  old <- wb_enter("Writing and validating raster", files = c(path, sources(x))); on.exit(.wb_log$context <- old, add = TRUE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = "partial_", tmpdir = dirname(path), fileext = ".tif")
  on.exit(unlink(c(tmp, paste0(tmp, ".aux.xml"))), add = TRUE)
  gdal <- cfg$gdal_opts
  if (startsWith(datatype, "INT")) gdal <- sub("^PREDICTOR=3$", "PREDICTOR=2", gdal)
  writeRaster(x, tmp, overwrite = TRUE,
              wopt = list(datatype = datatype, gdal = gdal))
  check <- rast(tmp)
  if (nlyr(check) != nlyr(x) || !compareGeom(check, x, stopOnError = FALSE))
    stop("Written raster failed validation: ", path)
  # Force a complete read and compare the number of valid cells per band.
  # Zero valid cells is permitted only when the corresponding source band is
  # also entirely NA; it is not evidence of a corrupt raster on its own.
  source_counts <- tryCatch(terra::global(x, "notNA")[, 1],
                            error = function(e) stop("Cannot read source raster for ", path, ": ", conditionMessage(e)))
  written_counts <- tryCatch(terra::global(check, "notNA")[, 1],
                             error = function(e) stop("Cannot read written raster for ", path, ": ", conditionMessage(e)))
  if (length(source_counts) != nlyr(x) || length(written_counts) != nlyr(x) ||
      any(!is.finite(source_counts)) || any(!is.finite(written_counts)) ||
      any(source_counts != written_counts))
    stop("Written raster validity differs from source: ", path)
  missing_bands <- which(written_counts == 0)
  if (length(missing_bands)) {
    msg <- paste0("All-NA annual raster band(s) retained: ",
                  paste(names(check)[missing_bands], collapse = ", "), " in ", path)
    warning(msg, call. = FALSE)
  }
  rm(check)
  # Atomic replacement on my macOS filesystem. If unsupported, leave the
  # old file intact and report an error rather than deleting it first.
  if (!file.rename(tmp, path)) stop("Could not atomically replace ", path)
  rast(path)
}



# ---------------- Spatial preparation and unit diagnostics ----------------
wb_require_crs <- function(x, label) {
  s <- crs(x)
  if (is.na(s) || !nzchar(s) || grepl("ENGCRS\\[|LOCAL_CS\\[", s))
    stop(label, " has a missing/local engineering CRS. Supply the verified source CRS; do not assign the target CRS to unknown coordinates.")
  invisible(x)
}

wb_metadata_units <- function(r) {
  u <- tryCatch(units(r), error = function(e) "")
  paste(unique(u[nzchar(u)]), collapse = "; ")
}

wb_unit_kind <- function(u) {
  u <- tolower(trimws(u))
  if (u %in% c("mm", "millimeter", "millimeters", "millimetre", "millimetres")) return("mm")
  if (u %in% c("cm", "centimeter", "centimeters", "centimetre", "centimetres")) return("cm")
  if (u %in% c("m", "meter", "meters", "metre", "metres")) return("m")
  if (u %in% c("ft", "foot", "feet", "us survey foot", "ftus")) return("feet")
  if (u %in% c("m3/m3", "m^3/m^3", "cm3/cm3", "cm^3/cm^3", "fraction", "1")) return("fraction")
  "unknown"
}

wb_check_awc <- function(raw, clip, cfg) {
  if (nlyr(raw) != 1L) stop("AWC must be one layer of available water storage or an explicitly configured volumetric AWC layer.")
  u <- wb_metadata_units(raw)
  cat("AWC raster band units: ", if (nzchar(u)) u else "not declared", "\n", sep = "")
  cat("AWC configured units: ", cfg$awc_units, "; depth: ",
      if (is.finite(cfg$awc_depth_mm)) paste(cfg$awc_depth_mm, "mm") else "not supplied", "\n", sep = "")
  kind <- wb_unit_kind(u)
  if (kind != "unknown" && kind != cfg$awc_units)
    stop("Suspected AWC unit mismatch: metadata says '", u, "', AWC_UNITS says '", cfg$awc_units,
         "'. Verify the source, then correct AWC_UNITS; no automatic conversion was made.")
  if (!nzchar(u)) message("AWC units are assumed from configuration; no metadata proves a mismatch.")
  z <- global(clip, c("min", "max", "mean"), na.rm = TRUE)
  print(z)
  if (!all(is.finite(as.numeric(z[1,])))) stop("No usable AWC values in this AOI.")
  if (cfg$awc_units == "mm" && z[1,"max"] <= 1)
    warning("AWC values are all <= 1 while configured as mm. These could be volumetric fractions. Check source units/depth; values are not rescaled.")
  if (z[1,"min"] <= 0) warning("Nonpositive AWC cells will be excluded from the water balance.")
  scale <- switch(cfg$awc_units, mm = 1, cm = 10, m = 1000, fraction = cfg$awc_depth_mm)
  if (z[1,"max"] * scale > 2000)
    warning("Some AWC storage exceeds 2000 mm. Check units, source depth, and missing-value coding; no additional rescaling is applied.")
  if (cfg$awc_units == "fraction" && (z[1,"min"] < 0 || z[1,"max"] > 1))
    stop("Volumetric AWC must lie between 0 and 1 (apart from NA).")
  if (is.finite(cfg$awc_depth_mm) && z[1,"max"] * scale > cfg$awc_depth_mm)
    warning("AWC storage exceeds the configured profile depth. Check units and depth.")
  scale
}

wb_jobs <- function(cfg) {
  a <- if (!is.null(cfg$aoi_bbox)) {
    b <- cfg$aoi_bbox
    vect(ext(b[1], b[2], b[3], b[4]), crs = "EPSG:4326")
  } else vect(cfg$aoi_path)
  wb_require_crs(a, "AOI")
  if (!nrow(a)) stop("AOI has no features.")
  if (!cfg$run_aois_separately || nrow(a) == 1L) return(list(list(label = "all_aoi", aoi = a)))
  idx <- cfg$aoi_features_to_run
  if (identical(idx, "all")) idx <- seq_len(nrow(a))
  if (!is.numeric(idx) || anyNA(idx) || !length(idx) || any(idx != as.integer(idx)) ||
      any(idx < 1 | idx > nrow(a)) || anyDuplicated(idx)) stop("Invalid AOI_FEATURES_TO_RUN.")
  lapply(idx, function(i) {
    label <- if (cfg$aoi_label_field %in% names(a)) values(a)[[cfg$aoi_label_field]][i] else paste0("aoi_", i)
    list(label = paste0(sprintf("aoi_%02d_", i), wb_sanitize_label(label)), aoi = a[i,])
  })
}

wb_prepare_site_uncached <- function(cfg, job) {
  old <- wb_enter("Preparing terrain and soil", "not applicable", c(cfg$dem_path, cfg$aoi_path, cfg$awc_path)); on.exit(.wb_log$context <- old, add = TRUE)
  cat("\nPreparing ", job$label, "\n", sep = "")
  out_dir <- cfg$run_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw <- wb_step("Opening DEM", cfg$dem_path, rast(cfg$dem_path))
  if (nlyr(raw) != 1) stop("DEM must have one layer.")
  if (!is.null(cfg$dem_source_crs)) {
    warning("Applying user-specified DEM_SOURCE_CRS metadata override: ", cfg$dem_source_crs)
    crs(raw) <- cfg$dem_source_crs
  }
  wb_require_crs(raw, "DEM")
  a <- fillHoles(job$aoi)
  if (!same.crs(a, raw)) a <- project(a, raw)
  # Match v2.2: mask before factor-3 aggregation, then project and mask again.
  raw <- mask(crop(raw, a), a)
  if (cfg$dem_aggregate_fact > 1)
    raw <- aggregate(raw, fact = cfg$dem_aggregate_fact, fun = mean, na.rm = TRUE)
  z <- global(raw, c("min", "max", "mean"), na.rm = TRUE)
  print(z)
  u <- wb_metadata_units(raw)
  chosen <- cfg$dem_z_units
  if (chosen == "auto") {
    chosen <- wb_infer_dem_z_units(z)
    message("DEM elevation heuristic selected ", chosen,
            ". Horizontal CRS units are separate; override DEM_Z_UNITS if source metadata disagrees.")
  }
  uk <- wb_unit_kind(u)
  if (uk %in% c("m", "feet") && uk != if (chosen == "meters") "m" else "feet")
    warning("DEM band units ('", u, "') conflict with the selected elevation units ('", chosen, "'). Verify DEM_Z_UNITS.")
  if (chosen == "feet") raw <- raw * .3048
  z_m <- global(raw, c("min", "max", "mean"), na.rm = TRUE)
  if (z_m[1,"max"] > 4500 || z_m[1,"min"] < -500)
    warning("Converted elevations are unusual for this western-US workflow. Check DEM_Z_UNITS and nodata coding.")
  dem <- project(raw, cfg$metric_crs, method = "bilinear")
  a <- project(a, cfg$metric_crs)
  dem <- mask(crop(dem, a), a)
  if (!any(is.finite(as.numeric(global(dem, "mean", na.rm = TRUE)[1,1])))) stop("DEM/AOI do not overlap.")
  # Transform each DEM cell center; no polygon centroid or feature order enters.
  latitude <- wb_latitude_raster(dem)
  latitude_range <- as.numeric(global(latitude, c("min", "max"), na.rm = TRUE)[1,])
  lat_rad <- latitude * pi / 180
  slope <- terrain(dem, "slope", unit = "degrees") * pi / 180
  aspect <- terrain(dem, "aspect", unit = "degrees")
  # Edge cells (no complete 3x3 neighbourhood) are treated as level ground.
  slope <- ifel(is.na(slope) & !is.na(dem), 0, slope)
  aspect <- ifel(is.na(aspect) & !is.na(dem), 0, aspect)
  names(slope) <- "slope_rad"; names(aspect) <- "aspect_deg"
  Af <- abs(180 - abs(aspect - 225)) * pi / 180
  hl_raw <- clamp(.339 + .808*cos(lat_rad)*cos(slope) - .196*sin(lat_rad)*sin(slope) -
                    .482*cos(Af)*sin(slope), 0, 2, values = TRUE)
  # Divide by the flat-surface value at each cell's own latitude: level ground
  # is exactly 1 everywhere, so only slope and aspect modify PET (day length in
  # the PET equation already carries latitude).
  hl <- clamp(hl_raw / (.339 + .808*cos(lat_rad)), .2, 2, values = TRUE)
  names(hl) <- "HL"
  hl_range <- as.numeric(global(hl, c("min", "max"), na.rm = TRUE)[1,])
  cat("Heat-load modifier (1 = level ground):", paste(round(hl_range, 3), collapse = " to "),
      "; cell latitude range:", paste(round(latitude_range, 5), collapse = " to "), "\n")
  awc <- wb_step("Opening AWC", cfg$awc_path, rast(cfg$awc_path))
  wb_require_crs(awc, "AWC")
  a_awc <- project(a, awc)
  ac <- mask(crop(awc, a_awc, snap = "out"), a_awc)
  scale <- wb_check_awc(awc, ac, cfg)
  soil <- project(ac, dem, method = "bilinear")
  if (scale != 1) soil <- soil * scale
  soil <- mask(soil, dem)
  soil <- ifel(soil > 0, soil, NA)
  names(soil) <- "SOILmax_mm"
  if (!is.finite(global(soil, "mean", na.rm = TRUE)[1,1])) stop("No positive AWC within the DEM/AOI.")
  compareGeom(dem, soil, hl, stopOnError = TRUE)
  metadata <- list(dem_z_units = chosen, awc_configured_units = cfg$awc_units,
    awc_metadata_units = wb_metadata_units(awc), awc_depth_mm = cfg$awc_depth_mm,
    hl_normalization = "flat surface at each cell's latitude", hl_range = hl_range, latitude_method = "individual WGS84 cell centers", latitude_range_deg = latitude_range, resolution = res(dem),
    pet_method = cfg$pet_method, pt_alpha = cfg$pt_alpha, surface_albedo = cfg$surface_albedo,
    rad_diffuse_fraction = cfg$rad_diffuse_fraction, rad_step_minutes = cfg$rad_step_minutes,
    cloud_transmittance = cfg$cloud_transmittance, humidity_source = cfg$humidity_source,
    wind_path = cfg$wind_path, wind_constant = cfg$wind_constant, aero_scale = cfg$aero_scale,
    snow_t_all = cfg$snow_t_all, rain_t_all = cfg$rain_t_all, melt_mm_per_deg_day = cfg$melt_mm_per_deg_day,
    melt_radiation_weight = cfg$melt_radiation_weight,
    radiation = "clear-sky direct beam on the tilted surface integrated over the mid-month day, isotropic diffuse by sky view, terrain-reflected diffuse; no horizon shading")
  saveRDS(metadata, file.path(out_dir, "site_metadata.rds"))
  writeLines(capture.output(str(metadata)), file.path(out_dir, "site_metadata.txt"))
  list(label = job$label, out_dir = out_dir, dem = dem, dem_z_units = chosen,
       SOILmax = soil, HL = hl, slope_rad = slope, aspect_deg = aspect,
       latitude = latitude,
       pet_method = cfg$pet_method, pt_alpha = cfg$pt_alpha, surface_albedo = cfg$surface_albedo,
       aero_scale = cfg$aero_scale,
       snow_t_all = cfg$snow_t_all, rain_t_all = cfg$rain_t_all,
       melt_factor = cfg$melt_mm_per_deg_day, melt_rad_weight = cfg$melt_radiation_weight,
       melt_method = if (cfg$melt_method == "energy") 2L else 1L, melt_rad_factor = cfg$melt_mm_per_mj,
       snow_albedo = cfg$snow_albedo,
       solar_declination = 0.409 * sin(2 * pi * (c(15,45,74,105,135,166,196,227,258,288,319,349) - 80) / 365),
       month_days = c(31,28,31,30,31,30,31,31,30,31,30,31))
}

# Preserve the established adjustment everywhere it is valid. Only repair
# otherwise eligible cells lost when the AOI-masked DEM crosses the PRISM grid.
wb_repair_lapse_edges <- function(site, cfg, dem_low, dz) {
  needs <- is.na(dz) & !is.na(site$dem) & (site$SOILmax > 0)
  n_needed <- global(needs, "sum", na.rm = TRUE)[1,1]
  if (!is.finite(n_needed) || n_needed == 0) return(dz)
  old <- wb_enter("Repairing DEM climate-adjustment edges", "not applicable", cfg$dem_path)
  on.exit(.wb_log$context <- old, add = TRUE)
  raw <- rast(cfg$dem_path)
  if (!is.null(cfg$dem_source_crs)) crs(raw) <- cfg$dem_source_crs
  if (cfg$dem_aggregate_fact > 1)
    raw <- aggregate(raw, fact = cfg$dem_aggregate_fact, fun = mean, na.rm = TRUE)
  if (site$dem_z_units == "feet") raw <- raw * .3048
  # Keep source elevations outside the AOI for this fallback only. Extend the
  # existing grid to the already buffered PRISM footprint; do not choose a new
  # model resolution, recompute terrain/heat load, or refit lapse coefficients.
  footprint <- project(as.polygons(ext(dem_low), crs = crs(dem_low)), site$dem)
  padded_grid <- extend(rast(site$dem), ext(footprint), snap = "out")
  surrounding <- project(raw, padded_grid, method = "bilinear")
  low_surrounding <- project(surrounding, dem_low, method = "bilinear")
  fallback <- site$dem - project(low_surrounding, site$dem, method = "bilinear")
  repair <- ifel(is.na(needs), FALSE, needs) & is.finite(fallback)
  n_fixed <- global(repair, "sum", na.rm = TRUE)[1,1]
  result <- ifel(repair, fallback, dz, wopt = list(datatype = "FLT8S"))
  names(result) <- names(dz)
  note <- sprintf("DEM edge correction: restored climate adjustments at %d of %d eligible cells; existing valid adjustments retained.",
                  as.integer(n_fixed), as.integer(n_needed))
  cat(note, "\n"); wb_note(note)
  result
}

wb_add_lapse <- function(site, cfg, index) {
  old <- wb_enter("Fitting lapse rates", "normals", index$file[index$is_normal]); on.exit(.wb_log$context <- old, add = TRUE)
  cat("\nPreparing PRISM normals for lapse-rate fitting...\n")
  # Crop normals around the AOI bounding box with a multi-cell interpolation
  # margin. No national-size values() allocation is needed for a local AOI.
  # p is snapped to the PRISM grid; the other fields are cropped to it with nearest snapping,
  # which keeps matching cell boundaries without interpolation (wb_climate_year).
  normals <- wb_climate_year(index, site, cfg)
  p <- normals$ppt; tmean <- normals$tmean
  dem_low <- project(site$dem, p[[1]], method = "bilinear")
  dz <- site$dem - project(dem_low, site$dem, method = "bilinear")
  # Deterministic subsampling if the AOI contains >200,000 coarse cells.
  set.seed(cfg$lapse_seed)
  coef <- wb_fit_lapse_from_normals(tmean, p, values(dem_low))
  if (any(!is.finite(unlist(coef)))) stop("Lapse fitting produced non-finite coefficients.")
  site$dz <- wb_repair_lapse_edges(site, cfg, dem_low, dz)
  # Elevation slopes: regional moving-window fits (LAPSE_WINDOW_KM), or the
  # single West-wide fit when the window is Inf.
  site$lapse_t <- coef$lapse_t
  site$beta_p <- coef$beta_p
  report <- data.frame(month = 1:12, lapse_C_per_m = coef$lapse_t, ppt_log_per_m = coef$beta_p)
  if (is.finite(cfg$lapse_window_km)) {
    cat("Fitting regional elevation slopes (", cfg$lapse_window_km, "-km window)...\n", sep = "")
    local <- wb_fit_lapse_local(tmean, p, dem_low, coef, cfg)
    # Interpolate to the DEM; any cell beyond the fitted blocks keeps the West-wide slope.
    to_dem <- function(x, b0) mask(cover(project(x, site$dem, method = "bilinear"),
                                         wb_make_12(site$dem * 0) + b0), site$dem)
    site$lapse_t <- to_dem(local$lapse_t, coef$lapse_t)
    site$beta_p <- to_dem(local$beta_p, coef$beta_p)
    mapped <- mask(c(local$lapse_t, local$beta_p), local$inside, maskvalues = 0)
    spread <- function(x, label) {
      q <- as.matrix(global(x, quantile, probs = c(.05, .5, .95), na.rm = TRUE))
      colnames(q) <- paste0("regional_", label, c("_p05", "_median", "_p95"))
      q
    }
    report <- cbind(report, spread(mapped[[1:12]], "lapse"), spread(mapped[[13:24]], "ppt"))
    names(report)[2:3] <- c("west_wide_lapse_C_per_m", "west_wide_ppt_log_per_m")
    wb_atomic_raster(mapped,
                     file.path(cfg$output_root, "LAPSE_regional_coefficients.tif"), cfg)
  }
  site$normals <- list(ppt = p, tmean = tmean, vpd = normals$vpd)
  write.csv(report, file.path(site$out_dir, "lapse_coefficients.csv"), row.names = FALSE)
  site$keep <- integer(12)
  site$keep[cfg$months_out] <- 1L
  site$static <- wb_static_normals(site, cfg, index, p)
  site$wind <- wb_read_wind(site, cfg)
  site
}

# Static 2 m wind climatology for Penman-Monteith (v6.0): a 12-layer raster on any
# grid (e.g. gridMET monthly means from calibrate_pt_alpha_vpd.R), cropped around
# the site. Kept on its own grid and interpolated per cell in the engine, like the
# cloud-transmittance and dew-point normals. NULL when a constant wind is used.
wb_read_wind <- function(site, cfg) {
  if (cfg$pet_method != "penman_monteith" || is.null(cfg$wind_path)) return(NULL)
  old <- wb_enter("Opening wind climatology", files = cfg$wind_path); on.exit(.wb_log$context <- old, add = TRUE)
  w <- wb_step("Opening wind climatology", cfg$wind_path, rast(cfg$wind_path))
  wb_require_crs(w, "wind climatology")
  if (nlyr(w) != 12) stop("WIND_PATH must hold 12 monthly layers of 2 m wind speed (m/s); found ", nlyr(w))
  w <- wb_crop_climate(w, site)
  names(w) <- paste0("u2_", sprintf("%02d", 1:12))
  rng <- as.numeric(as.matrix(global(w, c("min", "max"), na.rm = TRUE)))
  if (!any(is.finite(rng))) stop("The wind climatology has no values over the study area: ", cfg$wind_path)
  if (min(rng, na.rm = TRUE) < 0 || max(rng, na.rm = TRUE) > 25)
    stop("WIND_PATH values are expected to be 2 m wind speeds in m/s (0-25); found ", paste(round(range(rng, na.rm = TRUE), 2), collapse = " to "))
  cat("Wind climatology (2 m, m/s) over the site: ", paste(round(range(rng, na.rm = TRUE), 2), collapse = " to "), "\n", sep = "")
  w
}

# Static PRISM normals (cloud transmittance, dew point) on the cropped
# precipitation-normal window. A normal on a different grid is resampled to it.
wb_static_normals <- function(site, cfg, index, p) {
  if (!length(cfg$static_types)) return(NULL)
  old <- wb_enter("Opening static PRISM normals", "normals", index$file[index$is_normal & index$type %in% cfg$static_types])
  on.exit(.wb_log$context <- old, add = TRUE)
  out <- list()
  for (type in cfg$static_types) {
    x <- wb_climate_stack(index, type)
    x <- if (compareGeom(x, p, stopOnError = FALSE, ext = FALSE, rowcol = FALSE)) crop(x, p, snap = "near") else {
      warning("PRISM ", type, " normals are on a different grid than the precipitation normals; resampling bilinearly.")
      project(wb_crop_climate(x, site), p, method = "bilinear")
    }
    if (!compareGeom(x, p, stopOnError = FALSE)) stop("Could not align PRISM ", type, " normals with the climate window.")
    rng <- as.numeric(as.matrix(global(x, c("min", "max"), na.rm = TRUE)))
    if (type == "soltrans" && (min(rng, na.rm = TRUE) < 0 || max(rng, na.rm = TRUE) > 1.0001))
      stop("soltrans normals are expected to be a 0-1 transmittance; found ", paste(round(range(rng), 3), collapse = " to "))
    if (type == "tdmean" && (min(rng, na.rm = TRUE) < -60 || max(rng, na.rm = TRUE) > 40))
      stop("tdmean normals are expected in degrees C; found ", paste(round(range(rng), 1), collapse = " to "))
    cat("PRISM ", type, " normals: ", paste(round(range(rng), 3), collapse = " to "), "\n", sep = "")
    out[[type]] <- x
  }
  out
}

# ---------------- Non-destructive PRISM discovery and acquisition ----------------
wb_parse_prism <- function(file) {
  b <- tolower(basename(file))
  txt <- paste(tolower(basename(dirname(file))), b)
  grab <- function(pattern, text = b) {
    h <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
    if (length(h) > 1) h[-1] else character()
  }
  typ <- grab("prism_([^_[:space:]]+)_")
  if (!length(typ)) typ <- grab("prism_([^_[:space:]]+)_", txt)
  typ <- if (length(typ)) typ[1] else NA_character_
  resolution <- if (grepl("_30s_|_800m", txt)) "800m" else
    if (grepl("_25m_|_4km", txt)) "4km" else
      if (grepl("_15s_|_400m", txt)) "400m" else NA_character_
  region <- grab("prism_[^_]+_(us|ak|hi)_", txt)
  region <- if (length(region)) region[1] else "us" # legacy CONUS names
  normal <- grepl("avg_30y|30yr|normal", txt)
  year <- month <- normal_end <- NA_integer_
  h <- grab("_([12][0-9]{3})([01][0-9])_avg_30y", txt)
  if (length(h)) {
    normal <- TRUE; normal_end <- as.integer(h[1]); month <- as.integer(h[2])
  } else {
    h <- grab("_([12][0-9]{3})([01][0-9])(?:_|\\.|$)")
    if (length(h) && !normal) {
      year <- as.integer(h[1]); month <- as.integer(h[2])
    } else {
      h <- grab("_([01][0-9])(?:_bil)?\\.(tif|bil|zip)$")
      if (length(h)) {normal <- TRUE; month <- as.integer(h[1])}
    }
  }
  data.frame(file = file, type = typ, resolution = resolution, region = region,
             year = year, month = month, is_normal = normal, normal_end = normal_end)
}

wb_scan_prism <- function(cfg, zip = FALSE) {
  files <- unique(unlist(lapply(cfg$prism_search_dirs, function(d)
    list.files(d, pattern = if (zip) "\\.zip$" else "\\.(tif|bil)$",
               recursive = TRUE, full.names = TRUE, ignore.case = TRUE)), use.names = FALSE))
  if (!length(files)) return(wb_parse_prism("")[FALSE,])
  d <- do.call(rbind, lapply(files, wb_parse_prism))
  d <- d[which(d$type %in% c(cfg$climate_types, cfg$static_types) & d$resolution == cfg$prism_res &
                 d$region == "us" & d$month %in% 1:12),,drop = FALSE]
  if (!nrow(d)) return(d)
  # Static types (soltrans, tdmean) are used only as normals; climate types (ppt, tmean and, for
  # Penman-Monteith, vpdmax and vpdmin) as normals and as monthly series (v6.0).
  d <- d[which((d$is_normal & (is.na(d$normal_end) | d$normal_end == cfg$prism_normal_end_year)) |
      (!d$is_normal & d$type %in% cfg$climate_types & d$year %in% c(cfg$base_years, cfg$years_annual))),,drop = FALSE]
  d
}

wb_select_prism <- function(cfg, index) {
  if (!nrow(index)) return(index)
  index$mtime <- as.numeric(file.info(index$file)$mtime)
  # Prefer known requested normal period to unknown legacy period; then newest
  # file, then TIFF over BIL, then path for reproducible tie-breaking.
  index <- index[order(is.na(index$normal_end), -index$mtime,
                       !grepl("\\.tif$", index$file, ignore.case = TRUE), index$file),]
  key <- paste(index$type, index$is_normal, ifelse(index$is_normal, 0, index$year), index$month)
  index <- index[!duplicated(key),,drop = FALSE]
  index[order(index$type, !index$is_normal, index$year, index$month),,drop = FALSE]
}

wb_missing_prism <- function(cfg, index) {
  years <- sort(unique(c(cfg$base_years, cfg$years_annual)))
  needed <- expand.grid(type = cfg$climate_types, year = c(NA_integer_, years),
                        month = 1:12, stringsAsFactors = FALSE)
  if (length(cfg$static_types))
    needed <- rbind(needed, expand.grid(type = cfg$static_types, year = NA_integer_,
                                        month = 1:12, stringsAsFactors = FALSE))
  key <- function(d) paste(d$type, ifelse(is.na(d$year), "normal", d$year), d$month)
  have <- index
  have$year[have$is_normal] <- NA_integer_
  needed[!key(needed) %in% key(have),,drop = FALSE]
}

wb_extract_prism_zips <- function(cfg, missing) {
  z <- wb_scan_prism(cfg, zip = TRUE)
  if (!nrow(z) || !nrow(missing)) return(invisible(NULL))
  key <- function(d) paste(d$type, ifelse(is.na(d$year), "normal", d$year), d$month)
  z$year[z$is_normal] <- NA_integer_
  z <- z[key(z) %in% key(missing),,drop = FALSE]
  for (f in z$file) {
    entries <- tryCatch(utils::unzip(f, list = TRUE), error = function(e) NULL)
    if (is.null(entries) || !nrow(entries)) {warning("Unreadable ZIP retained: ", f); next}
    if (any(grepl("(^[/\\\\])|(^[A-Za-z]:)|(^|[/\\\\])\\.\\.([/\\\\]|$)", entries$Name)))
      stop("Unsafe ZIP entry path in ", f)
    dest <- file.path(cfg$prism_dir_v2, "_v3_extracted", tools::file_path_sans_ext(basename(f)))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(f, exdir = dest, overwrite = FALSE)
  }
  invisible(NULL)
}

# Optional compact PRISM archive. PRISM's public web service returns ZIPs
# containing COGs, rather than directly addressable COG URLs. This mode saves
# persistent disk, NOT necessarily download bandwidth. Full archives are
# downloaded into an ephemeral staging directory, cropped, and discarded.
wb_prism_subset_extent <- function(cfg, raster) {
  # Expand the AOI in a metric CRS before projection, including the entire
  # lapse-fitting moving window and four coarse interpolation cells.
  if (!is.null(cfg$aoi_bbox)) {
    a <- terra::as.polygons(terra::ext(cfg$aoi_bbox), crs = "EPSG:4326")
  } else a <- terra::vect(cfg$aoi_path)
  if (!nzchar(terra::crs(a))) stop("AOI CRS is missing: cannot subset PRISM safely.")
  a <- terra::project(a, "EPSG:5070")
  # 4 extra PRISM pixels ensures that boundary interpolation has neighbours.
  margin <- max(cfg$prism_buffer_km, cfg$lapse_window_km) * 1000 +
    4 * if (cfg$prism_res == "4km") 4000 else if (cfg$prism_res == "800m") 800 else 400
  area <- terra::as.polygons(terra::extend(terra::ext(a), margin), crs = "EPSG:5070")
  terra::ext(terra::project(area, terra::crs(raster)))
}

wb_ensure_prism_subset <- function(cfg) {
  dir.create(cfg$prism_subset_dir, recursive = TRUE, showWarnings = FALSE)
  index <- wb_select_prism(cfg, wb_scan_prism(cfg))
  missing <- wb_missing_prism(cfg, index)
  if (!nrow(missing)) return(index)
  if (!cfg$allow_downloads)
    stop("PRISM subset is missing ", nrow(missing), " variable/month combinations; downloads disabled.")
  # One year/type group at a time to limit peak temporary disk usage.
  groups <- split(missing, paste(missing$type, ifelse(is.na(missing$year), "normals", missing$year)))
  for (part in groups) {
    # Recheck to skip groups already completed by a prior interrupted run.
    outstanding <- wb_missing_prism(cfg, wb_select_prism(cfg, wb_scan_prism(cfg)))
    key <- function(d) paste(d$type, d$year, d$month)
    part <- part[key(part) %in% key(outstanding), , drop = FALSE]
    if (!nrow(part)) next
    stage <- tempfile("prism_subset_staging_")
    dir.create(stage, recursive = TRUE)
    tryCatch({
      prism::prism_set_dl_dir(stage)
      type <- part$type[1]
      if (is.na(part$year[1])) {
        prism::get_prism_normals(type = type, resolution = cfg$prism_res,
                                mon = unique(part$month), keepZip = TRUE)
      } else {
        prism::get_prism_monthlys(type = type, years = part$year[1],
                                 mon = unique(part$month), resolution = cfg$prism_res,
                                 keepZip = TRUE)
      }
      tmp_cfg <- cfg
      tmp_cfg$prism_search_dirs <- stage
      tmp_cfg$prism_dir_v2 <- stage
      wb_extract_prism_zips(tmp_cfg, part)
      sources <- wb_scan_prism(tmp_cfg)
      key_part <- paste(part$type, ifelse(is.na(part$year), "normal", part$year), part$month)
      key_src <- paste(sources$type, ifelse(sources$is_normal, "normal", sources$year), sources$month)
      sources <- sources[key_src %in% key_part, , drop = FALSE]
      if (!nrow(sources)) stop("PRISM archive contained no usable TIFFs for ", type)
      for (i in seq_len(nrow(sources))) {
        src <- sources$file[i]
        x <- terra::rast(src)
        window <- wb_prism_subset_extent(cfg, x)
        clipped <- terra::crop(x, window, snap = "out")
        if (is.null(clipped) || !terra::ncell(clipped))
          stop("PRISM spatial subset does not overlap AOI: ", src)
        dst <- file.path(cfg$prism_subset_dir, paste0(tools::file_path_sans_ext(basename(src)), "_subset.tif"))
        tempdst <- paste0(dst, ".partial.tif")
        terra::writeRaster(clipped, tempdst, overwrite = TRUE,
                           wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES")))
        if (!file.rename(tempdst, dst)) stop("Cannot finalize PRISM subset: ", dst)
      }
    }, finally = unlink(stage, recursive = TRUE, force = TRUE))
  }
  index <- wb_select_prism(cfg, wb_scan_prism(cfg))
  missing <- wb_missing_prism(cfg, index)
  if (nrow(missing)) stop("PRISM subset cache incomplete; ", nrow(missing), " combinations missing.")
  cat("PRISM spatial subset complete: ", nrow(index), " rasters at ", cfg$prism_subset_dir, "\n", sep = "")
  index
}

wb_ensure_prism <- function(cfg) {
  if (identical(cfg$prism_mode, "subset")) return(wb_ensure_prism_subset(cfg))
  index <- wb_select_prism(cfg, wb_scan_prism(cfg))
  missing <- wb_missing_prism(cfg, index)
  if (nrow(missing)) {
    # Existing archives are retained; extraction creates only missing files.
    wb_extract_prism_zips(cfg, missing)
    index <- wb_select_prism(cfg, wb_scan_prism(cfg))
    missing <- wb_missing_prism(cfg, index)
  }
  if (nrow(missing) && !cfg$allow_downloads)
    stop("PRISM cache is missing ", nrow(missing), " month/variable combinations. Enable ALLOW_DOWNLOADS or supply the missing data.")
  if (nrow(missing)) {
    dir.create(cfg$prism_dir_v2, recursive = TRUE, showWarnings = FALSE)
    prism::prism_set_dl_dir(cfg$prism_dir_v2)
    for (type in unique(missing$type)) {
      m <- missing[missing$type == type,,drop = FALSE]
      if (anyNA(m$year))
        prism::get_prism_normals(type = type, resolution = cfg$prism_res,
                                mon = m$month[is.na(m$year)], keepZip = TRUE)
      annual <- m[!is.na(m$year),,drop = FALSE]
      for (yr in unique(annual$year))
        prism::get_prism_monthlys(type = type, years = yr, mon = annual$month[annual$year == yr],
                                 resolution = cfg$prism_res, keepZip = TRUE)
    }
    # Support package versions that retain ZIPs without extracting them.
    index <- wb_select_prism(cfg, wb_scan_prism(cfg))
    wb_extract_prism_zips(cfg, wb_missing_prism(cfg, index))
    index <- wb_select_prism(cfg, wb_scan_prism(cfg))
    missing <- wb_missing_prism(cfg, index)
  }
  if (nrow(missing)) {
    hint <- paste0(if (any(missing$type %in% c("soltrans", "tdmean")))
      "\nThe soltrans/tdmean normals need a prism package version that offers these types (>= 0.2.1); alternatively set CLOUD_TRANSMITTANCE / HUMIDITY_SOURCE to numbers." else "",
      if (any(missing$type %in% c("vpdmax", "vpdmin")))
      "\nThe vpdmax/vpdmin series are needed by PET_METHOD = \"penman_monteith\" (prism package >= 0.2.1 offers them); check the download folder is writable, or choose another PET method." else "")
    stop("PRISM archive remains incomplete:\n", paste(capture.output(head(missing, 24)), collapse = "\n"), hint)
  }
  if (any(index$is_normal & is.na(index$normal_end)))
    warning("Some legacy normals do not encode their reference-period end year. Confirm their period in prism_inputs.csv and the source metadata.")
  cat("PRISM input rasters selected:", nrow(index), "; cache preserved.\n")
  index
}

wb_climate_stack <- function(index, type, year = NULL) {
  old <- wb_enter(paste("Opening PRISM", type), if (is.null(year)) "normals" else year); on.exit(.wb_log$context <- old, add = TRUE)
  use <- index$type == type
  use <- use & if (is.null(year)) index$is_normal else !index$is_normal & !is.na(index$year) & index$year == year
  d <- index[which(use),,drop = FALSE]
  if (nrow(d) != 12L || !setequal(d$month, 1:12)) stop("Expected exactly 12 months for ", type, " ", year)
  d <- d[order(d$month),]
  # Open each month individually so header warnings identify a single file.
  x <- rast(lapply(d$file, function(f) wb_step("Opening PRISM month", f, rast(f))))
  names(x) <- paste0(type, "_", sprintf("%02d", 1:12))
  wb_require_crs(x, paste("PRISM", type))
  x
}



# PRISM vpdmax and vpdmin are supplied in hPa; convert their mean to kPa
# for Penman-Monteith (1 kPa = 10 hPa). Normals when year is NULL.
# Leave cached PRISM rasters in their original units.
wb_vpd_stack <- function(index, year = NULL) {
  hi <- wb_climate_stack(index, "vpdmax", year); lo <- wb_climate_stack(index, "vpdmin", year)
  if (!compareGeom(hi, lo, stopOnError = FALSE)) stop("PRISM vpdmax and vpdmin grids differ (", if (is.null(year)) "normals" else year, ").")
  x <- (hi + lo) / 20  # mean of two hPa fields, converted to kPa
  names(x) <- paste0("vpd_", sprintf("%02d", 1:12))
  x
}

# Open a monthly climate set (ppt, tmean and, for Penman-Monteith, vpd) cropped to the site window.
wb_climate_year <- function(index, site, cfg, year = NULL) {
  p <- wb_crop_climate(wb_climate_stack(index, "ppt", year), site)
  tmean <- crop(wb_climate_stack(index, "tmean", year), p, snap = "near")
  if (!compareGeom(p, tmean, stopOnError = FALSE))
    stop("PRISM precipitation and temperature grids differ (", if (is.null(year)) "normals" else year, ").")
  vpd <- NULL
  if ("vpdmax" %in% cfg$climate_types) {
    vpd <- crop(wb_vpd_stack(index, year), p, snap = "near")
    if (!compareGeom(p, vpd, stopOnError = FALSE))
      stop("PRISM vapour-pressure deficit grid differs from the precipitation grid (", if (is.null(year)) "normals" else year, ").")
    rng <- as.numeric(as.matrix(global(vpd, c("min", "max"), na.rm = TRUE)))
    if (min(rng, na.rm = TRUE) < -0.01 || max(rng, na.rm = TRUE) > 12)
      stop("Converted PRISM vapour-pressure deficit is expected in kPa (0-12); found ", paste(round(range(rng), 2), collapse = " to "))
  }
  list(ppt = p, tmean = tmean, vpd = vpd)
}

# Climate-grid positions are prepared once per distinct grid and reused.
wb_grid_for <- function(ppt, tmean, site, mem, grids) {
  key <- wb_grid_key(ppt)
  if (is.null(grids[[key]])) grids[[key]] <- wb_memory_grid(mem, site, c(ppt[[1]], tmean[[1]]), site$dz)
  grids[[key]]
}

wb_spinup <- function(site, cfg, mem, grids) {
  old <- wb_enter("Spin-up", "normals"); on.exit(.wb_log$context <- old, add = TRUE)
  path <- file.path(site$out_dir, "state", paste0("SPINUP_", cfg$month_tag, ".rds"))
  key <- list(run = cfg$run_id, site = site$label, stage = "spinup")
  cached <- if (cfg$resume) wb_read_state(path, key) else NULL
  if (!is.null(cached)) {
    cat("Using saved spin-up state.\n")
    if (!isTRUE(cached$info$converged)) warning("Resuming from a spin-up that did not meet the configured convergence criteria; see spinup_diagnostics.csv.")
    return(cached)
  }
  grid <- wb_grid_for(site$normals$ppt, site$normals$tmean, site, mem, grids)
  soil <- mem$awc
  pack <- numeric(mem$n)
  records <- list()
  gain_streak <- numeric(mem$n)
  converged <- FALSE
  for (i in seq_len(cfg$spinup_max_iter)) {
    x <- wb_year(mem, grid, site$normals$ppt, site$normals$tmean, soil, pack, site, vpd = site$normals$vpd)
    ds <- wb_change_stats(soil, x[, "soil_end"], cfg$soil_tol_mm)
    dp <- wb_change_stats(pack, x[, "pack_end"], cfg$pack_tol_mm)
    records[[i]] <- data.frame(iteration = i, soil_mean_mm = ds$mean, soil_max_mm = ds$max,
      soil_pct_above_tol = ds$pct, pack_mean_mm = dp$mean, pack_max_mm = dp$max,
      pack_pct_above_tol = dp$pct, compared_cells = min(ds$count, dp$count),
      newly_invalid_cells = max(ds$lost, dp$lost))
    cat(sprintf("Spin-up %d: soil mean %.4f/max %.3f mm (%.4f%% > tolerance); snow mean %.4f/max %.3f mm (%.4f%% > tolerance)\n",
                i, ds$mean, ds$max, ds$pct, dp$mean, dp$max, dp$pct))
    snow_change <- x[, "pack_end"] - pack
    gain_streak <- ifelse(snow_change > cfg$pack_tol_mm, gain_streak + 1, 0)
    soil <- x[, "soil_end"]; pack <- x[, "pack_end"]
    converged <- all(is.finite(c(ds$mean, dp$mean, ds$pct, dp$pct))) &&
      ds$count > 0 && dp$count > 0 && ds$lost == 0 && dp$lost == 0 &&
      ds$mean <= cfg$soil_tol_mm && dp$mean <= cfg$pack_tol_mm &&
      ds$pct <= cfg$spinup_max_pct && dp$pct <= cfg$spinup_max_pct
    if (converged) break
  }
  write.csv(do.call(rbind, records), file.path(site$out_dir, "spinup_diagnostics.csv"), row.names = FALSE)
  if (cfg$write_rasters) wb_atomic_raster(wb_to_raster(mem, cbind(ds$delta, dp$delta), c("soil_change_mm", "pack_change_mm")),
                   file.path(site$out_dir, "spinup_final_change_mm.tif"), cfg)
  status <- if (converged) "Spin-up met configured mean and cell-fraction tolerances." else
    "Spin-up did NOT meet configured tolerances. Inspect residual maps and diagnostics; persistent snow accumulation may not stabilize with more iterations."
  writeLines(c(status, paste("Iterations:", i), paste("Allowed percent above each 1-cell tolerance:", cfg$spinup_max_pct)),
             file.path(site$out_dir, "spinup_status.txt"))
  wb_note(c(status,
    paste("Spin-up diagnostics directory:", site$out_dir), capture.output(print(do.call(rbind, records)))))
  # Snow anomalies remain visible without becoming an additional stopping rule.
  snow_path <- file.path(cfg$output_root, paste0("SNOW_spinup_baseline_",
    min(cfg$base_years), "_", max(cfg$base_years), ".tif"))
  if (cfg$write_rasters) wb_atomic_raster(wb_to_raster(mem, cbind(pack, snow_change),
    c("spinup_end_snowpack_SWE_mm", "last_spinup_year_change_SWE_mm")), snow_path, cfg)
  change_range <- range(snow_change, na.rm = TRUE)
  persistent_pct <- 100 * mean(gain_streak >= 3, na.rm = TRUE)
  snow_summary <- c(sprintf("Snow diagnostic: final annual change %.3f to %.3f mm SWE; %.4f%% of valid cells gained >%.3f mm for at least 3 consecutive spin-up years.",
    change_range[1], change_range[2], persistent_pct, cfg$pack_tol_mm),
    if (cfg$write_rasters) paste("Snow map (end pack and signed annual change):", snow_path) else "Snow map not written (WRITE_RASTERS is FALSE).",
    "Persistent gains are diagnostic, not proof of an error or a glacier. Spin-up-year count affects accumulated storage.")
  cat(paste(snow_summary, collapse = "\n"), "\n")
  wb_note(snow_summary)
  if (!converged) {
    if (cfg$spinup_action == "stop") stop(status) else warning(status)
  }
  cat(status, "\n")
  wb_save_state(path, list(key = key, soil = soil, pack = pack,
                           info = list(converged = converged, iterations = i)))
}

# PET self-check over the climate normals: monthly Hamon x heat load (v4),
# Priestley-Taylor and, when the vapour-pressure deficit normals are present,
# Penman-Monteith PET from the spun-up state, plus the radiation table. Written
# whatever PET_METHOD drives the run, so a new area can be checked at once.
# Monthly tables come from a subsample of cells; annual maps cover every cell.
wb_pet_check <- function(site, cfg, mem, grids, spin) {
  old <- wb_enter("PET self-check", "normals"); on.exit(.wb_log$context <- old, add = TRUE)
  grid <- wb_grid_for(site$normals$ppt, site$normals$tmean, site, mem, grids)
  tag <- paste0("normals_", cfg$prism_normal_end_year - 29L, "_", cfg$prism_normal_end_year)
  # Monthly PET on a subsample (all cells when the area is small).
  set.seed(cfg$lapse_seed)
  sub <- if (mem$n > 500000L) sort(sample.int(mem$n, 500000L)) else seq_len(mem$n)
  take <- function(v) v[sub]; takem <- function(m) m[sub, , drop = FALSE]
  s <- mem$static
  mem_sub <- list(n = length(sub), awc = take(mem$awc), hl = take(mem$hl), latitude = take(mem$latitude),
    elevation = take(mem$elevation), rad = takem(mem$rad), lapse_t = takem(mem$lapse_t), beta_p = takem(mem$beta_p),
    static = wb_static_subset(s, sub))
  grid_sub <- list(key = grid$key, nr = grid$nr, nc = grid$nc, rowf = take(grid$rowf), colf = take(grid$colf),
                   covered = take(grid$covered), dz = take(grid$dz))
  methods <- c("hamon", "priestley_taylor", if (!is.null(site$normals$vpd)) "penman_monteith")
  monthly <- function(method) wb_year(mem_sub, grid_sub, site$normals$ppt, site$normals$tmean,
                                      take(spin$soil), take(spin$pack), site, method = method, diag = 2L,
                                      vpd = site$normals$vpd)[, 5:16, drop = FALSE]
  pet_h <- colMeans(monthly("hamon"), na.rm = TRUE); pet_p <- colMeans(monthly("priestley_taylor"), na.rm = TRUE)
  report <- cbind(mem$rad_summary, transmittance_fill = s$trans_fill,
    pet_hamon_mm = pet_h, pet_priestley_taylor_mm = pet_p, pet_ratio_pt_over_hamon = pet_p / pet_h)
  if ("penman_monteith" %in% methods) {
    pet_m <- colMeans(monthly("penman_monteith"), na.rm = TRUE)
    report <- cbind(report, pet_penman_monteith_mm = pet_m, pet_ratio_pm_over_pt = pet_m / pet_p)
  }
  csv <- file.path(cfg$output_root, paste0("PET_check_", tag, ".csv"))
  write.csv(report, csv, row.names = FALSE)
  # Annual PET, AET and CWD for every cell, one method at a time.
  annual <- NULL; maps <- list()
  for (method in methods) {
    x <- wb_year(mem, grid, site$normals$ppt, site$normals$tmean, spin$soil, spin$pack, site, method = method, diag = 1L,
                 vpd = site$normals$vpd)
    annual <- rbind(annual, data.frame(method = method, PET_mm = mean(x[, "pet"], na.rm = TRUE),
      AET_mm = mean(x[, "aet"], na.rm = TRUE), CWD_mm = mean(x[, "cwd"], na.rm = TRUE)))
    maps[[method]] <- x[, "pet"]
    rm(x); gc(verbose = FALSE)
  }
  write.csv(annual, file.path(cfg$output_root, paste0("PET_check_annual_", tag, ".csv")), row.names = FALSE)
  layers <- cbind(maps$hamon, maps$priestley_taylor, maps$priestley_taylor / maps$hamon)
  layer_names <- c("PET_hamon_mm", "PET_priestley_taylor_mm", "ratio_pt_over_hamon")
  if (!is.null(maps$penman_monteith)) {
    layers <- cbind(layers, maps$penman_monteith, maps$penman_monteith / maps$priestley_taylor)
    layer_names <- c(layer_names, "PET_penman_monteith_mm", "ratio_pm_over_pt")
  }
  wb_atomic_raster(wb_to_raster(mem, layers, layer_names),
    file.path(cfg$output_root, paste0("PET_check_annual_", tag, ".tif")), cfg)
  rm(layers)
  rm(maps); gc(verbose = FALSE)
  wb_write_terrain_modifier(site, cfg, file.path(cfg$output_root, "RAD_terrain_modifier_by_month.tif"))
  lines <- c(sprintf("PET self-check over the normals (monthly table from %s cells; annual means over all cells and MONTHS_OUT):",
                     format(length(sub), big.mark = ",")),
             capture.output(print(report, row.names = FALSE, digits = 4)),
             capture.output(print(annual, row.names = FALSE, digits = 5)),
             paste("Written:", csv))
  cat(paste(lines, collapse = "\n"), "\n"); wb_note(lines)
  invisible(report)
}

# Monthly output files: <output_root>/monthly/<VAR>_monthly_<unit>_<year>.tif,
# 12 layers each (mm; soil and pack are end-of-month states; ppt is the
# downscaled monthly precipitation and snowfall its snow share; vpd is the
# downscaled vapour-pressure deficit in kPa, Penman-Monteith runs only).
wb_monthly_vars <- c("cwd", "aet", "pet", "water", "soil", "pack", "ppt", "snowfall", "tmean", "rs", "rnet", "melt", "vpd", "runoff")
wb_monthly_units <- c(cwd = "mm", aet = "mm", pet = "mm", water = "mm", soil = "mm", pack = "mm", ppt = "mm", snowfall = "mm",
                      tmean = "C", rs = "MJm2d", rnet = "MJm2d", melt = "mm", vpd = "kPa", runoff = "mm")
wb_monthly_paths <- function(cfg, vars, year) {
  file.path(cfg$output_root, "monthly", paste0(toupper(vars), "_monthly_", wb_monthly_units[vars], "_", year, ".tif"))
}

wb_write_monthly <- function(mem, x, cfg, year, vars = cfg$monthly_output) {
  old <- wb_enter("Writing monthly rasters", year); on.exit(.wb_log$context <- old, add = TRUE)
  for (v in vars) {
    cols <- paste0(v, "_m", sprintf("%02d", 1:12))
    if (!all(cols %in% colnames(x))) stop("Monthly columns for ", v, " are missing from the engine output.")
    vals <- x[, cols]
    # Whole millimetres for water variables when INT2S is chosen; temperature and radiation stay Float32.
    dtype <- if (wb_monthly_units[v] == "mm") cfg$monthly_datatype else "FLT4S"
    if (dtype == "INT2S") vals <- round(pmin(vals, 32767))   # Int16 halves the size and compresses well
    wb_atomic_raster(wb_to_raster(mem, vals, paste0(toupper(v), "_", year, "_", sprintf("%02d", 1:12))),
                     wb_monthly_paths(cfg, v, year), cfg, dtype)
  }
  invisible(NULL)
}

# Optional stacked monthly product per variable: 12 x N layers over the years
# that have monthly files. Large; off unless MONTHLY_STACK is TRUE.
wb_monthly_stacks <- function(cfg) {
  if (!length(cfg$monthly_output) || !cfg$monthly_stack) return(invisible(NULL))
  old <- wb_enter("Stacking monthly rasters"); on.exit(.wb_log$context <- old, add = TRUE)
  for (v in cfg$monthly_output) {
    files <- wb_monthly_paths(cfg, v, sort(cfg$monthly_years_by_var[[v]]))
    files <- files[file.exists(files)]
    if (!length(files)) next
    x <- rast(lapply(files, function(f) wb_step("Reading monthly raster for stack", f, rast(f))))
    yrs <- range(as.integer(sub(".*_([0-9]{4})\\.tif$", "\\1", basename(files))))
    path <- file.path(cfg$output_root, paste0(toupper(v), "_monthly_", wb_monthly_units[v], "_", yrs[1], "_", yrs[2], ".tif"))
    wb_atomic_raster(x, path, cfg)
    wb_note(paste("Monthly stack:", path)); cat("Monthly stack:", path, "\n")
  }
  invisible(NULL)
}

# ---------------- Points (v5.5) ----------------
# Points come from a CSV with longitude/latitude columns (WGS84) or any vector
# file terra can read. Each point is matched to its DEM cell; with a
# neighbourhood of 3 or 5 the cells around it are averaged. Every monthly
# variable the engine knows and the annual CWD/AET sums are written to
# <output_root>/points/points_monthly.csv and points_annual.csv.
wb_read_points <- function(cfg, site) {
  old <- wb_enter("Reading points", files = cfg$points_path); on.exit(.wb_log$context <- old, add = TRUE)
  ext_ <- tolower(tools::file_ext(cfg$points_path))
  if (ext_ %in% c("csv", "txt")) {
    d <- read.csv(cfg$points_path, stringsAsFactors = FALSE, check.names = FALSE)
    nm <- tolower(trimws(names(d)))
    lon_col <- which(nm %in% c("lon", "long", "longitude", "x", "lng"))[1]
    lat_col <- which(nm %in% c("lat", "latitude", "y"))[1]
    if (is.na(lon_col) || is.na(lat_col))
      stop("POINTS_PATH CSV needs longitude and latitude columns (names like lon/longitude/x and lat/latitude/y).")
    lon <- suppressWarnings(as.numeric(d[[lon_col]])); lat <- suppressWarnings(as.numeric(d[[lat_col]]))
    id_col <- if (!is.null(cfg$points_id_field)) which(names(d) == cfg$points_id_field)[1] else setdiff(seq_along(d), c(lon_col, lat_col))[1]
    if (!is.null(cfg$points_id_field) && is.na(id_col)) stop("POINTS_ID_FIELD '", cfg$points_id_field, "' is not a column of ", cfg$points_path)
    ids <- if (is.na(id_col)) seq_len(nrow(d)) else d[[id_col]]
    attrs <- d[, setdiff(seq_along(d), c(lon_col, lat_col, id_col)), drop = FALSE]
    v <- vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
  } else {
    v <- vect(cfg$points_path)
    wb_require_crs(v, "points")
    if (!all(geomtype(v) == "points")) v <- centroids(v)
    d <- as.data.frame(v)
    id_col <- if (!is.null(cfg$points_id_field)) which(names(d) == cfg$points_id_field)[1] else if (ncol(d)) 1L else NA
    if (!is.null(cfg$points_id_field) && is.na(id_col)) stop("POINTS_ID_FIELD '", cfg$points_id_field, "' is not a field of ", cfg$points_path)
    ids <- if (is.na(id_col)) seq_len(nrow(v)) else d[[id_col]]
    attrs <- if (ncol(d)) d[, setdiff(seq_along(d), id_col), drop = FALSE] else data.frame(row.names = seq_len(nrow(v)))
    ll <- crds(project(v, "EPSG:4326")); lon <- ll[, 1]; lat <- ll[, 2]
  }
  if (anyDuplicated(ids)) stop("Point ids must be unique (", sum(duplicated(ids)), " duplicates).")
  xy <- crds(project(v, crs(site$dem)))
  cell <- cellFromXY(site$dem, xy)
  pts <- data.frame(point_id = as.character(ids), longitude = lon, latitude = lat, cell = cell, stringsAsFactors = FALSE)
  pts$elevation_m <- NA_real_; ok <- !is.na(cell); pts$elevation_m[ok] <- values(site$dem, mat = FALSE)[cell[ok]]
  half <- (cfg$points_neighbourhood - 1L) %/% 2L
  rc <- rowColFromCell(site$dem, ifelse(is.na(cell), 1L, cell))
  nb <- do.call(cbind, lapply(-half:half, function(dr) lapply(-half:half, function(dc) {
    r <- rc[, 1] + dr; cc <- rc[, 2] + dc
    good <- !is.na(cell) & r >= 1 & r <= nrow(site$dem) & cc >= 1 & cc <= ncol(site$dem)
    out <- rep(NA_integer_, nrow(pts)); out[good] <- cellFromRowCol(site$dem, r[good], cc[good]); out })))
  nb <- matrix(as.integer(unlist(nb)), nrow = nrow(pts))
  cat("Points read:", nrow(pts), " (", sum(!is.na(cell)), "inside the DEM extent)\n")
  list(pts = pts, attrs = attrs, nb_cells = nb)
}

# Match the neighbourhood cells to engine rows; points with no modelled cell are reported and dropped.
wb_points_rows <- function(points, mem, cfg) {
  rows <- matrix(match(points$nb_cells, mem$cells), nrow = nrow(points$pts))
  has <- rowSums(!is.na(rows)) > 0
  dropped <- points$pts[!has, c("point_id", "longitude", "latitude", "elevation_m")]
  if (nrow(dropped)) {
    dropped$reason <- ifelse(is.na(points$pts$cell[!has]), "outside the DEM extent", "no modelled cell (outside the AOI or no soil data)")
    dir.create(file.path(cfg$output_root, "points"), recursive = TRUE, showWarnings = FALSE)
    write.csv(dropped, file.path(cfg$output_root, "points", "points_dropped.csv"), row.names = FALSE)
    msg <- sprintf("%d of %d points have no modelled cell and are listed in points/points_dropped.csv.", nrow(dropped), nrow(points$pts))
    cat(msg, "\n"); wb_note(msg)
  }
  if (!any(has)) stop("No point falls on a modelled cell.")
  points$pts <- points$pts[has, ]; points$attrs <- points$attrs[has, , drop = FALSE]
  points$nb_cells <- points$nb_cells[has, , drop = FALSE]
  points$rows <- rows[has, , drop = FALSE]                   # rows of mem (NA where a neighbour cell is not modelled)
  points$sub <- sort(unique(as.vector(points$rows[!is.na(points$rows)])))  # distinct engine rows the points need
  points$sub_rows <- matrix(match(points$rows, points$sub), nrow = nrow(points$rows))  # rows within the subset
  points
}

# The in-memory inputs for a subset of rows, so the engine can be run on the
# point cells alone (full-raster runs) at negligible cost.
wb_mem_subset <- function(mem, idx) {
  list(n = length(idx), cells = mem$cells[idx], template = mem$template, awc = mem$awc[idx], hl = mem$hl[idx],
       latitude = mem$latitude[idx], elevation = mem$elevation[idx], rad = mem$rad[idx, , drop = FALSE],
       lapse_t = mem$lapse_t[idx, , drop = FALSE], beta_p = mem$beta_p[idx, , drop = FALSE],
       static = wb_static_subset(mem$static, idx))
}
wb_grid_subset <- function(grid, idx) list(key = grid$key, nr = grid$nr, nc = grid$nc, rowf = grid$rowf[idx],
                                           colf = grid$colf[idx], covered = grid$covered[idx], dz = grid$dz[idx])

# Point values for one year from an engine result that holds every monthly
# variable (rows = points$sub when `subset` is TRUE, otherwise all of mem).
wb_points_year <- function(x, points, year, stage, cfg, subset = TRUE) {
  rows <- if (subset) points$sub_rows else points$rows
  avg <- function(col) { v <- matrix(x[as.vector(rows), col], nrow = nrow(rows)); m <- rowMeans(v, na.rm = TRUE); m[is.nan(m)] <- NA; m }
  monthly <- do.call(rbind, lapply(1:12, function(m) {
    out <- data.frame(point_id = points$pts$point_id, stage = stage, year = year, month = m, stringsAsFactors = FALSE)
    for (v in wb_monthly_vars) out[[v]] <- avg(paste0(v, "_m", sprintf("%02d", m)))
    out }))
  annual <- data.frame(point_id = points$pts$point_id, stage = stage, year = year,
                       cwd = avg("cwd"), aet = avg("aet"), stringsAsFactors = FALSE)
  list(monthly = monthly, annual = annual)
}

# Assemble every year's point results into the two CSVs.
wb_points_write <- function(cfg, site, points, stage_dirs) {
  old <- wb_enter("Writing point results"); on.exit(.wb_log$context <- old, add = TRUE)
  files <- unlist(lapply(stage_dirs, function(d) list.files(d, pattern = "^points_[0-9]{4}\\.rds$", full.names = TRUE)))
  res <- lapply(files, readRDS)
  monthly <- do.call(rbind, lapply(res, `[[`, "monthly")); annual <- do.call(rbind, lapply(res, `[[`, "annual"))
  meta <- cbind(points$pts[, c("point_id", "longitude", "latitude", "elevation_m")], points$attrs)
  dir.create(file.path(cfg$output_root, "points"), recursive = TRUE, showWarnings = FALSE)
  monthly <- merge(meta, monthly, by = "point_id"); annual <- merge(meta, annual, by = "point_id")
  monthly <- monthly[order(monthly$point_id, monthly$year, monthly$month), ]
  annual <- annual[order(annual$point_id, annual$year), ]
  names(annual)[names(annual) %in% c("cwd", "aet")] <- paste0(c("cwd", "aet"), "_", cfg$month_tag, "_mm")
  fm <- file.path(cfg$output_root, "points", "points_monthly.csv"); fa <- file.path(cfg$output_root, "points", "points_annual.csv")
  write.csv(monthly, fm, row.names = FALSE); write.csv(annual, fa, row.names = FALSE)
  writeLines(c("points_monthly.csv: one row per point, year and month.",
    "  cwd, aet, pet: monthly deficit, actual and potential evapotranspiration (mm)",
    "  water: rain + snowmelt reaching the soil (mm); soil, pack: end-of-month soil water and snowpack SWE (mm)",
    "  ppt, snowfall: downscaled precipitation and its snow share (mm); melt: snowmelt (mm)",
    "  tmean: downscaled mean temperature (C); rs: incoming shortwave on the cell's slope under the month's cloudiness (MJ m-2 d-1)",
    "  rnet: net radiation at a snow surface, (1 - SNOW_ALBEDO) x rs - net longwave (MJ m-2 d-1)",
    "  vpd: downscaled monthly vapour-pressure deficit (kPa); NA unless PET_METHOD is penman_monteith",
    sprintf("points_annual.csv: CWD and AET summed over MONTHS_OUT (%s), matching the annual rasters.", cfg$month_tag),
    sprintf("Values are the mean of a %d x %d cell neighbourhood around each point (POINTS_NEIGHBOURHOOD); stage is baseline or analysis.",
            cfg$points_neighbourhood, cfg$points_neighbourhood),
    "Points without a modelled cell are listed in points_dropped.csv."), file.path(cfg$output_root, "points", "README.txt"))
  wb_note(c(paste("Point results:", fm), paste("Point results:", fa)))
  cat("Point results:", fm, "\n")
  invisible(NULL)
}

# Persist the converged end-of-baseline state independently of the requested
# output variables, output directory, and zonal geography. Changes to the
# scientific inputs invalidate the cache. No intermediate spin-up rasters are
# retained here.
wb_initialization_cache <- function(cfg, index, site, mem) {
  science_fields <- c(
    "base_years", "dem_aggregate_fact", "dem_z_units", "dem_source_crs",
    "metric_crs", "awc_units", "awc_depth_mm", "prism_res",
    "prism_normal_end_year", "lapse_window_km", "lapse_t_limits",
    "lapse_p_limits", "lapse_seed", "spinup_max_iter", "soil_tol_mm",
    "pack_tol_mm", "spinup_max_pct", "spinup_action", "pet_method",
    "pt_alpha", "wind_constant", "aero_scale", "surface_albedo",
    "rad_diffuse_fraction", "rad_step_minutes", "cloud_transmittance",
    "humidity_source", "snow_t_all", "rain_t_all", "melt_method",
    "melt_mm_per_deg_day", "melt_mm_per_mj", "snow_albedo",
    "melt_radiation_weight", "code_md5"
  )
  historical <- index[is.na(index$year) | index$is_normal |
                        (!is.na(index$year) & index$year %in% cfg$base_years), , drop = FALSE]
  # Source climate files, not output settings, determine initialization identity.
  identity <- list(version = 1L, environment = list(inputs = wb_file_stamp(wb_input_files(cfg)),
                                        n = mem$n, grid = c(terra::nrow(site$dem), terra::ncol(site$dem),
                                         terra::res(site$dem), as.vector(terra::ext(site$dem)))),
                   science = cfg[intersect(science_fields, names(cfg))],
                   baseline_climate = wb_file_stamp(historical$file),
                   packages = vapply(c("terra", "prism", "Rcpp", "RcppParallel"),
                                     function(x) as.character(utils::packageVersion(x)), character(1)))
  key <- wb_hash(identity)
  root <- cfg$environment_cache
  if (is.null(root) || !nzchar(root))
    root <- file.path(dirname(normalizePath(cfg$dem_path, mustWork = TRUE)), ".wbr_prepared_environment")
  path <- file.path(path.expand(root), "initialization_states", key, "baseline_end.rds")
  list(key = key, path = path, identity = identity)
}

wb_cached_baseline <- function(path, key, n, last_year) {
  x <- wb_read_state(path, key)
  if (is.null(x) || !identical(x$year, as.integer(last_year))) return(NULL)
  for (name in c("soil", "pack", "snow_autumn", "snow_autumn_peak")) {
    if (!is.numeric(x[[name]]) || length(x[[name]]) != n) return(NULL)
  }
  if (any(!is.finite(x$soil)) || any(!is.finite(x$pack))) return(NULL)
  x
}

wb_run_sequence <- function(site, cfg, index, years, initial, stage, mem, grids, points = NULL, zonal = NULL) {
  emit_annual <- isTRUE(cfg$write_rasters) && (stage != "baseline" || isTRUE(cfg$write_baseline))
  if (stage == "baseline" && emit_annual && any(cfg$outputs_to_write %in% c("snow", "pack", "snow_peak")) && is.null(initial$snow_autumn)) stop("Baseline seasonal snow requires preceding autumn state; disable write_baseline or supply prior autumn values.")

  old <- wb_enter(stage); on.exit(.wb_log$context <- old, add = TRUE)
  dir <- file.path(site$out_dir, paste0(stage, "_years"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, paste0("WB_", cfg$month_tag, "_mm_", years, ".tif"))
  # Monthly rasters (one file per variable and year) go straight to the output
  # folder; they are part of each year's completion record for resume checks.
  vars_for_year <- function(y) cfg$monthly_output[vapply(cfg$monthly_output, function(v) y %in% cfg$monthly_years_by_var[[v]], logical(1))]
  monthly_paths <- function(j) { v <- vars_for_year(years[j]); if (length(v)) wb_monthly_paths(cfg, v, years[j]) else character() }
  zonal_paths <- function(j) if (is.null(zonal)) character() else file.path(dir, paste0("zonal_", years[seq_len(j)], ".rds"))
  points_paths <- if (is.null(points)) function(j) character() else function(j) file.path(dir, paste0("points_", years[seq_len(j)], ".rds"))
  stamp_paths <- function(j) c(if (emit_annual) paths[seq_len(j)], unlist(lapply(seq_len(j), monthly_paths)), points_paths(j), zonal_paths(j))
  all_monthly <- as.integer(sum(2^(seq_along(wb_monthly_vars) - 1L)))
  # One restart file per stage holds the soil and snow carried out of the last
  # finished year. It is trusted only if this run wrote it and every annual
  # raster it lists is unchanged.
  state_path <- file.path(dir, "state.rds")
  key <- list(run = cfg$run_id, site = site$label, stage = stage)
  state <- if (cfg$resume) wb_read_state(state_path, key) else NULL
  if (!is.null(state) && !isTRUE(tryCatch(state$year %in% years &&
      identical(state$stamp, wb_file_stamp(stamp_paths(sum(years <= state$year)))), error = function(e) FALSE))) state <- NULL
  if (is.null(state)) state <- list(key = key, year = min(years) - 1L, soil = initial$soil, pack = initial$pack)
  else cat("Resuming ", stage, " after December ", state$year, "\n", sep = "")
  for (j in which(years > state$year)) {
    wb_check_space(.wb_log$temp_dir)
    yr <- years[j]
    wb_enter(paste("Calculating", stage), yr,
      index$file[!index$is_normal & !is.na(index$year) & index$year == yr])
    cat("Calculating ", stage, " ", yr, "\n", sep = "")
    clim <- wb_climate_year(index, site, cfg, yr)
    p <- clim$ppt; tmean <- clim$tmean; vpd <- clim$vpd
    mvars <- vars_for_year(yr)
    zvars <- if (is.null(zonal)) character() else cfg$zonal_monthly
    annual_vars <- unique(c(if (emit_annual) ifelse(cfg$outputs_to_write %in% c("snow", "snow_peak"), "pack", cfg$outputs_to_write) else character(),
                            if (stage == "baseline") "pack" else character()))
    needed <- unique(c(mvars, zvars, annual_vars))
    grid <- wb_grid_for(p, tmean, site, mem, grids)
    chunks <- if (!cfg$write_rasters) list(needed) else if (length(needed))
      split(needed, ceiling(seq_along(needed) / cfg$monthly_chunk)) else list(character())
    mask_of <- function(v) as.integer(sum(2^(match(v, wb_monthly_vars) - 1L)))
    x <- wb_year(mem, grid, p, tmean, state$soil, state$pack, site,
                 monthly = mask_of(chunks[[1]]), vpd = vpd)
    annual_values <- list()
    prev_snow <- state$snow_autumn
    next_snow <- NULL
    # Preserve original annual CWD/AET values exactly as calculated by the engine.
    if (emit_annual && "cwd" %in% cfg$outputs_to_write) annual_values$cwd <- x[, "cwd"]
    if (emit_annual && "aet" %in% cfg$outputs_to_write) annual_values$aet <- x[, "aet"]
    zrows <- list()
    next_peak <- NULL
    collect <- function(dat, chunk) {
      requested <- intersect(chunk, annual_vars)
      for (v in requested) {
        cols <- paste0(v, "_m", sprintf("%02d", 1:12))
        if (!all(cols %in% colnames(dat))) stop("Missing annual monthly columns: ", v)
        a <- dat[, cols, drop = FALSE]
        if (v == "pack") {
          next_snow <<- rowSums(a[, 9:12, drop = FALSE])
          if ("snow" %in% cfg$outputs_to_write || "pack" %in% cfg$outputs_to_write) {
            val <- if (is.null(prev_snow)) rep(NA_real_, nrow(a)) else
              (prev_snow + rowSums(a[, 1:5, drop = FALSE])) / 9
            for (nm in intersect(c("snow", "pack"), cfg$outputs_to_write)) annual_values[[nm]] <<- val
          }
          if ("snow_peak" %in% cfg$outputs_to_write) {
            # Peak from preceding autumn requires the previous months' individual states.
            val <- if (is.null(state$snow_autumn_peak)) rep(NA_real_, nrow(a)) else
              pmax(state$snow_autumn_peak, apply(a[, 1:5, drop = FALSE], 1, max))
            annual_values$snow_peak <<- val
          }
          next_peak <<- apply(a[, 9:12, drop = FALSE], 1, max)
        } else if (v == "soil") {
          annual_values$soil <<- rowMeans(a[, 3:10, drop = FALSE])
        } else if (!v %in% c("cwd", "aet")) {
          annual_values[[v]] <<- rowSums(a[, cfg$months_out, drop = FALSE])
        }
      }
      output <- intersect(chunk, mvars)
      if (length(output)) wb_write_monthly(mem, dat, cfg, yr, output)
      selected <- intersect(chunk, zvars)
      if (length(selected)) zrows[[length(zrows) + 1L]] <<- wb_zonal_chunk(dat, selected, zonal, yr, stage)
    }
    collect(x, chunks[[1]])
    for (ch in chunks[-1]) {
      xc <- wb_year(mem, grid, p, tmean, state$soil, state$pack, site, monthly = mask_of(ch), vpd = vpd)
      collect(xc, ch); rm(xc); gc(verbose = FALSE)
    }
    if (emit_annual) {
      missing <- setdiff(cfg$outputs_to_write, names(annual_values))
      if (length(missing)) stop("Annual outputs not collected: ", paste(missing, collapse = ", "))
      wb_atomic_raster(wb_to_raster(mem, do.call(cbind, annual_values[cfg$outputs_to_write]),
                                      cfg$outputs_to_write), paths[j], cfg)
    }
    if (!is.null(zonal)) {
      if (!length(zrows)) stop("Zonal analysis has no variables to collect.")
      wb_save_state(zonal_paths(j)[j], do.call(rbind, zrows))
    }
    if (!is.null(points)) {
      if (cfg$write_rasters || !is.null(zonal)) {
        # Full run: repeat the year on the point cells alone with every monthly variable (same inputs, same state).
        xp <- wb_year(wb_mem_subset(mem, points$sub), wb_grid_subset(grid, points$sub), p, tmean,
                      state$soil[points$sub], state$pack[points$sub], site, awc = mem$awc[points$sub], monthly = all_monthly, vpd = vpd)
        saveRDS(wb_points_year(xp, points, yr, stage, cfg, subset = TRUE), points_paths(j)[j])
      } else saveRDS(wb_points_year(x, points, yr, stage, cfg, subset = FALSE), points_paths(j)[j])
    }
    if (cfg$diagnostics) print(apply(x, 2, function(v) c(min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE), mean = mean(v, na.rm = TRUE))))
    state$year <- yr; state$soil <- x[, "soil_end"]; state$pack <- x[, "pack_end"]
    state$snow_autumn <- next_snow
    state$snow_autumn_peak <- next_peak
    state$stamp <- wb_file_stamp(stamp_paths(j))
    wb_save_state(state_path, state)
    wb_cleanup_temp(site)
  }
  list(final = state, paths = paths)
}

wb_summaries <- function(site, cfg, baseline, analysis) {
  old <- wb_enter("Writing final summaries"); on.exit(.wb_log$context <- old, add = TRUE)
  summary_dir <- file.path(cfg$output_root, "summaries")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  for (output in cfg$outputs_to_write) {
    read_band <- function(f) wb_step("Reading annual result for summary", f, terra::rast(f)[[output]])
    for (stage in if (isTRUE(cfg$write_baseline)) c("baseline", "analysis") else "analysis") {
      files <- if (stage == "baseline") baseline$paths else analysis$paths
      years <- if (stage == "baseline") cfg$base_years else cfg$years_annual
      if (!all(file.exists(files))) stop("Missing ", stage, " yearly raster files: ", paste(files[!file.exists(files)], collapse = ", "))
      layers <- do.call(c, lapply(files, read_band))
      names(layers) <- paste0(toupper(output), "_", years)
      # Retain the original analysis annual brick, unchanged in its usual location.
      if (stage == "analysis") {
        path <- file.path(cfg$output_root, paste0(toupper(output), "_", cfg$month_tag, "_mm_", min(years), "_", max(years), ".tif"))
        wb_atomic_raster(layers, path, cfg)
        cat("Final stack:", path, "\n")
      }
      if (isTRUE(cfg$summarize)) {
        # Reporting-period means are independent of the spin-up/baseline purpose.
        # An incomplete first snow year is NA and excluded from the mean;
        # zero SWE is a valid observation and is included.
        mean_r <- terra::mean(layers, na.rm = TRUE)
        names(mean_r) <- paste0(output, "_", stage, "_mean_mm")
        path <- file.path(summary_dir, paste0(toupper(output), "_mean_", stage, "_", min(years), "_", max(years), ".tif"))
        wb_atomic_raster(mean_r, path, cfg)
        wb_note(paste("Period mean:", path))
      }
    }
  }
  invisible(TRUE)
}

wb_main_impl <- function(cfg) {
  wb_load_packages()
  cfg <- wb_validate(cfg)
  wb_setup_runtime(cfg)
  if (cfg$stop_after_setup) {
    cfg$run_dir <- file.path(cfg$work_root, "setup_check")
    for (job in wb_jobs(cfg)) wb_prepare_site(cfg, job)
    cat("Setup check finished; no climate download or annual simulation was run.\n")
    return(invisible(cfg))
  }
  index <- wb_ensure_prism(cfg)
  cfg <- wb_start_run(cfg, index)
  for (job in wb_jobs(cfg)) {
    site <- wb_prepare_site(cfg, job)
    wb_cleanup_temp(site)
    site <- wb_add_lapse(site, cfg, index)
    wb_cleanup_temp(site)
    points <- if (!is.null(cfg$points_path)) wb_read_points(cfg, site) else NULL
    mem <- wb_memory_site(site, cfg, keep_cells = if (!cfg$write_rasters && is.null(cfg$zonal_path) && !is.null(points)) unique(as.vector(points$nb_cells)) else NULL)
    mem$static <- wb_memory_static(mem, site, cfg)
    if (!is.null(points)) {
      points <- wb_points_rows(points, mem, cfg)
      cat("Points with modelled cells:", nrow(points$pts), if (!cfg$write_rasters) "(points-only run: no rasters will be written)" else "", "\n")
    }
    zonal <- if (is.null(cfg$zonal_path)) NULL else wb_zonal_prepare(cfg$zonal_path, cfg$zonal_id_field, mem)
    write.csv(mem$rad_summary, file.path(site$out_dir, "radiation_summary.csv"), row.names = FALSE)
    cat("Clear-sky radiation on slopes, monthly mean (MJ m-2 d-1):",
        paste(round(mem$rad_summary$clear_sky_slope_MJ_m2_d, 1), collapse = " "), "\n")
    cat("Terrain modifier median (p05-p95) by month:",
        paste0(sprintf("%.2f", mem$rad_summary$terrain_modifier_median), "(",
               sprintf("%.2f", mem$rad_summary$terrain_modifier_p05), "-",
               sprintf("%.2f", mem$rad_summary$terrain_modifier_p95), ")", collapse = " "), "\n")
    cat("PET method:", cfg$pet_method, "\n")
    grids <- new.env()
    init_cache <- wb_initialization_cache(cfg, index, site, mem)
    initial <- wb_cached_baseline(init_cache$path, init_cache$key, mem$n, max(cfg$base_years))
    # Baseline raster output is an explicit opt-in; cached states cannot replace
    # those requested historical output products.
    if (isTRUE(cfg$write_baseline)) initial <- NULL
    if (is.null(initial)) {
      spin <- wb_spinup(site, cfg, mem, grids)
      wb_cleanup_temp(site)
      if (cfg$pet_diagnostics) wb_pet_check(site, cfg, mem, grids, spin)
      wb_cleanup_temp(site)
      baseline <- wb_run_sequence(site, cfg, index, cfg$base_years, spin, "baseline", mem, grids, points, zonal)
      initial <- baseline$final
      # Seasonal outputs require preceding autumn snow information. Save it with
      # the initialized state rather than allowing output selection to erase it.
      if (!is.null(initial$snow_autumn) && !is.null(initial$snow_autumn_peak)) {
        cached <- initial
        cached$key <- init_cache$key
        cached$stamp <- NULL
        wb_save_state(init_cache$path, cached)
        cat("Saved reusable baseline-end state: ", init_cache$path, "\n", sep = "")
      } else {
        warning("Initialization state lacks preceding autumn snow history; cache not saved.")
      }
    } else {
      cat("Reusing initialized baseline state from ", init_cache$path, "\n", sep = "")
      baseline <- list(final = initial, paths = character())
    }
    cat("Analysis initialization: continuous from final baseline state\n")
    analysis <- wb_run_sequence(site, cfg, index, cfg$years_annual, initial, "analysis", mem, grids, points, zonal)
    if (cfg$write_rasters) { wb_summaries(site, cfg, baseline, analysis); wb_monthly_stacks(cfg) }
    if (!is.null(points)) wb_points_write(cfg, site, points, file.path(site$out_dir, c("baseline_years", "analysis_years")))
    if (!is.null(zonal)) wb_zonal_export(cfg, site,
      if (length(baseline$paths)) c("baseline", "analysis") else "analysis", zonal)
    for (f in c("site_metadata.txt", "lapse_coefficients.csv", "radiation_summary.csv", "spinup_status.txt", "spinup_diagnostics.csv")) {
      path <- file.path(site$out_dir, f)
      if (file.exists(path)) wb_note(c(paste("Diagnostics:", f), readLines(path)))
    }
    # Retain ending states locally for provenance/future extension; remove only
    # this run's intermediate checkpoint directories after final validation.
    if (cfg$write_rasters) {
      wb_atomic_raster(wb_to_raster(mem, cbind(baseline$final$soil, baseline$final$pack), c("soil_end", "pack_end")),
        file.path(cfg$run_dir, paste0("state_Dec_", max(cfg$base_years), ".tif")), cfg, "FLT8S")
      wb_atomic_raster(wb_to_raster(mem, cbind(analysis$final$soil, analysis$final$pack), c("soil_end", "pack_end")),
        file.path(cfg$run_dir, paste0("state_Dec_", max(cfg$years_annual), ".tif")), cfg, "FLT8S")
    }
    if (!cfg$diagnostics) {
      unlink(file.path(site$out_dir, c("baseline_years", "analysis_years", "state")), recursive = TRUE)
      unlink(file.path(site$out_dir, "spinup_final_change_mm.tif"))
    }
  }
  wb_note(paste("Completed:", Sys.time()))
  cat("\nCompleted. Outputs:", cfg$output_root, "\n")
  invisible(cfg)
}

# A single calling handler logs while the nested file/year context is still active.
wb_main <- function(cfg) {
  dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
  .wb_log$report <- file.path(cfg$output_root, "run_report.txt")
  .wb_log$context <- list(step = "startup", year = "not applicable", files = "not applicable")
  .wb_log$temp_dir <- NULL
  on.exit({
    if (!is.null(.wb_log$temp_dir)) unlink(.wb_log$temp_dir, recursive = TRUE)
    .wb_log$temp_dir <- NULL
    .wb_log$report <- NULL
  }, add = TRUE)
  wb_note(c("", paste("===== Run started:", Sys.time(), "====="),
    "Warnings are logged immediately. Multiple file paths indicate inputs involved in a combined operation.",
    "Raster read/write failures and logging timeouts abort processing. Completed-run checkpoints are pruned unless DIAGNOSTICS=TRUE."))
  withCallingHandlers(wb_main_impl(cfg),
    warning = function(w) {
      wb_condition(w, "WARNING")
      if (grepl("GDAL error|cannot read|cannot write|read failed|write failed|IReadBlock|TIFFRead|TIFFWrite|No space left|disk full|Operation timed out|Problem closing connection", conditionMessage(w), ignore.case = TRUE))
        stop(structure(list(message = paste("Raster I/O failure:", conditionMessage(w)), call = NULL),
                       class = c("wb_io_error", "error", "condition")))
      invokeRestart("muffleWarning")
    },
    error = function(e) wb_condition(e, "ERROR"))
}

wb_environment_input_signature <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  files <- path
  if (tolower(tools::file_ext(path)) == "shp") {
    stem <- tools::file_path_sans_ext(path)
    files <- Sys.glob(paste0(stem, ".*"))
  }
  info <- file.info(files)
  if (anyNA(info$size) || anyNA(info$mtime)) stop("Cannot fingerprint input: ", path)
  data.frame(path = normalizePath(files), bytes = as.numeric(info$size),
             modified = as.numeric(info$mtime), stringsAsFactors = FALSE)
}

wb_environment_key <- function(cfg, job) {
  # Keep this independent of output_root, run years, PRISM, and analysis settings.
  settings <- list(version = 1L,
    dem = wb_environment_input_signature(cfg$dem_path),
    awc = wb_environment_input_signature(cfg$awc_path),
    aoi = if (!is.null(cfg$aoi_bbox)) cfg$aoi_bbox else wb_environment_input_signature(cfg$aoi_path),
    job_label = job$label,
    aoi_features_to_run = cfg$aoi_features_to_run,
    run_aois_separately = cfg$run_aois_separately,
    aoi_label_field = cfg$aoi_label_field,
    aggregate = cfg$dem_aggregate_fact, z_units = cfg$dem_z_units,
    source_crs = cfg$dem_source_crs, metric_crs = cfg$metric_crs,
    awc_units = cfg$awc_units, awc_depth_mm = cfg$awc_depth_mm)
  temp <- tempfile(fileext = ".rds")
  on.exit(unlink(temp), add = TRUE)
  saveRDS(settings, temp, version = 2)
  unname(tools::md5sum(temp))
}

wb_prepare_site <- function(cfg, job) {
  if (identical(cfg$prepared_environment, FALSE))
    return(wb_prepare_site_uncached(cfg, job))
  cache_root <- cfg$environment_cache
  if (is.null(cache_root) || !nzchar(cache_root))
    cache_root <- file.path(dirname(normalizePath(cfg$dem_path, mustWork = TRUE)),
                            ".wbr_prepared_environment")
  key <- wb_environment_key(cfg, job)
  dest <- file.path(path.expand(cache_root), key)
  manifest <- file.path(dest, "manifest.rds")
  fields <- c("dem", "SOILmax", "HL", "slope_rad", "aspect_deg", "latitude")
  raster_paths <- file.path(dest, paste0(fields, ".tif"))
  if (file.exists(manifest) && all(file.exists(raster_paths))) {
    meta <- tryCatch(readRDS(manifest), error = function(e) NULL)
    if (is.list(meta) && identical(meta$key, key) &&
        identical(meta$fields, fields)) {
      site <- meta$site
      loaded <- tryCatch(lapply(raster_paths, terra::rast), error = function(e) NULL)
      if (!is.null(loaded) && all(vapply(loaded, function(x) terra::nlyr(x) == 1L, logical(1)))) {
        names(loaded) <- fields
        site[fields] <- loaded
        site$out_dir <- cfg$run_dir
        dir.create(site$out_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines(capture.output(str(meta$metadata)), file.path(site$out_dir, "site_metadata.txt"))
        saveRDS(meta$metadata, file.path(site$out_dir, "site_metadata.rds"))
        cat("Reusing prepared geography: ", dest, "\n", sep = "")
        return(site)
      }
    }
  }
  site <- wb_prepare_site_uncached(cfg, job)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  # Write files before publishing the manifest. An interrupted write is never
  # treated as a completed cache entry.
  for (i in seq_along(fields)) {
    tmp <- tempfile(pattern = paste0(fields[i], "_"), tmpdir = dest, fileext = ".tif")
    terra::writeRaster(site[[fields[i]]], tmp, overwrite = TRUE,
                       wopt = list(datatype = "FLT8S", gdal = cfg$gdal_opts))
    if (file.exists(raster_paths[i])) unlink(raster_paths[i])
    if (!file.rename(tmp, raster_paths[i])) stop("Cannot finalize prepared raster: ", raster_paths[i])
  }
  metadata <- readRDS(file.path(site$out_dir, "site_metadata.rds"))
  saved <- site
  saved[fields] <- NULL
  meta <- list(key = key, fields = fields, site = saved, metadata = metadata)
  tmp_manifest <- tempfile(pattern = "manifest_", tmpdir = dest, fileext = ".rds")
  saveRDS(meta, tmp_manifest)
  if (file.exists(manifest)) unlink(manifest)
  if (!file.rename(tmp_manifest, manifest)) stop("Cannot finalize prepared environment manifest")
  cat("Saved reusable prepared geography: ", dest, "\n", sep = "")
  site
}

wb_zonal_prepare <- function(zones, id_field, mem) {
  if (is.character(zones)) zones <- if (tolower(tools::file_ext(zones)) %in%
      c("tif", "tiff", "vrt", "grd", "img", "nc")) terra::rast(zones) else terra::vect(zones)
  if (!inherits(zones, c("SpatRaster", "SpatVector"))) stop("ZONAL_PATH must be a raster or vector.")
  template <- mem$template
  if (inherits(zones, "SpatRaster")) {
    if (terra::nlyr(zones) != 1L) stop("Zonal raster must have exactly one layer.")
    categories <- tryCatch(terra::cats(zones)[[1]], error = function(e) NULL)
    zid <- terra::project(zones, template, method = "near")
    values <- terra::values(zid, mat = FALSE)[mem$cells]
    if (any(is.finite(values) & values != floor(values))) stop("Zonal rasters require integer category codes.")
    id <- as.character(values)
    label_lookup <- if (!is.null(categories) && ncol(categories) >= 2L)
      setNames(as.character(categories[[2]]), as.character(categories[[1]])) else NULL
  } else {
    if (!terra::same.crs(zones, template)) zones <- terra::project(zones, template)
    typ <- terra::geomtype(zones)
    ids <- if (is.null(id_field)) as.character(seq_len(nrow(zones))) else {
      if (!id_field %in% names(zones)) stop("ZONAL_ID_FIELD not found: ", id_field)
      as.character(terra::values(zones)[[id_field]])
    }
    if (anyNA(ids) || any(!nzchar(ids))) stop("Missing or blank zone identifiers.")
    if (identical(typ, "polygons")) {
      # Overlapping polygons: the last feature wins, matching terra::rasterize.
      zones$wb_zone_code <- match(ids, unique(ids))
      zid <- terra::rasterize(zones, template, field = "wb_zone_code", touches = FALSE)
      codes <- terra::values(zid, mat = FALSE)[mem$cells]
      id <- unique(ids)[codes]
    } else if (identical(typ, "points")) {
      cells <- terra::cellFromXY(template, terra::crds(zones))
      ids_by_cell <- rep(NA_character_, terra::ncell(template))
      # Multiple points in a model cell are different sampling units; disallow ambiguity.
      if (anyDuplicated(cells[!is.na(cells)])) stop("Multiple zonal points share a model cell; use POINTS_PATH for overlapping point samples.")
      ids_by_cell[cells[!is.na(cells)]] <- ids[!is.na(cells)]
      id <- ids_by_cell[mem$cells]
    } else stop("Zonal vectors must contain points or polygons.")
    label_lookup <- NULL
  }
  good <- !is.na(id) & nzchar(id) & id != "NA"
  labels <- unique(id[good])
  if (!length(labels)) stop("No modeled cells intersect ZONAL_PATH.")
  group <- match(id, labels)
  meta <- data.frame(zone_id = labels, zone_label = labels,
                     n_zone_cells = tabulate(group[good], nbins = length(labels)),
                     stringsAsFactors = FALSE)
  if (!is.null(label_lookup)) {
    l <- label_lookup[meta$zone_id]
    meta$zone_label[!is.na(l)] <- l[!is.na(l)]
  }
  list(group = group, meta = meta)
}

wb_zonal_chunk <- function(x, vars, zonal, year, stage) {
  out <- vector("list", length(vars) * 12L); k <- 0L
  group <- zonal$group; ng <- nrow(zonal$meta)
  for (v in vars) for (m in 1:12) {
    col <- paste0(v, "_m", sprintf("%02d", m))
    if (!col %in% colnames(x)) stop("Missing zonal monthly result: ", col)
    value <- x[, col]
    valid <- which(!is.na(group) & is.finite(value))
    sums <- numeric(ng); sumsq <- numeric(ng); counts <- integer(ng)
    if (length(valid)) {
      z <- rowsum(value[valid], group[valid], reorder = TRUE)
      n <- tabulate(group[valid], nbins = ng)
      sums[as.integer(rownames(z))] <- z[, 1]
      z2 <- rowsum(value[valid]^2, group[valid], reorder = TRUE)
      sumsq[as.integer(rownames(z2))] <- z2[, 1]
      counts <- n
    }
    # Sample spatial variance (n - 1), computed from the cell-level
    # first and second moments. Never confuse temporal variability with
    # within-zone spatial variability; zero-valued snow is valid.
    mn <- ifelse(counts > 0L, sums / pmax(counts, 1L), NA_real_)
    spatial_var <- rep(NA_real_, ng)
    many <- which(counts > 1L)
    spatial_var[many] <- pmax(0, (sumsq[many] - sums[many]^2 / counts[many]) / (counts[many] - 1L))
    spatial_sd <- sqrt(spatial_var)
    spatial_cv <- ifelse(is.finite(mn) & abs(mn) > 0, spatial_sd / abs(mn), NA_real_)
    k <- k + 1L
    out[[k]] <- data.frame(zone_id = zonal$meta$zone_id, zone_label = zonal$meta$zone_label,
      stage = stage, year = year, month = m, variable = v,
      mean = mn, spatial_sd = spatial_sd, spatial_variance = spatial_var, spatial_cv = spatial_cv,
      n_cells = counts, n_zone_cells = zonal$meta$n_zone_cells,
      coverage = ifelse(zonal$meta$n_zone_cells > 0, counts / zonal$meta$n_zone_cells, NA_real_))
  }
  do.call(rbind, out)
}

wb_zonal_annual <- function(monthly) {
  if (is.null(monthly) || !nrow(monthly)) return(NULL)
  flux <- c("cwd", "aet", "pet", "water", "ppt", "snowfall", "melt")
  key <- interaction(monthly$zone_id, monthly$year, monthly$stage, monthly$variable, drop = TRUE)
  parts <- split(monthly, key)
  do.call(rbind, lapply(parts, function(d) {
    if (nrow(d) != 12L || any(!is.finite(d$mean))) val <- NA_real_
    else val <- if (d$variable[1] %in% flux) sum(d$mean) else mean(d$mean)
    data.frame(zone_id = d$zone_id[1], zone_label = d$zone_label[1], stage = d$stage[1],
      year = d$year[1], variable = d$variable[1], mean = val,
      n_cells = min(d$n_cells), n_zone_cells = d$n_zone_cells[1],
      coverage = min(d$coverage), temporal_stat = if (d$variable[1] %in% flux) "sum" else "mean")
  }))
}

wb_zonal_climatology <- function(monthly) {
  if (is.null(monthly) || !nrow(monthly)) return(NULL)
  d <- monthly[monthly$stage == "analysis", , drop = FALSE]
  if (!nrow(d)) return(NULL)
  key <- interaction(d$zone_id, d$month, d$variable, drop = TRUE)
  parts <- split(d, key)
  do.call(rbind, lapply(parts, function(x) {
    ok <- is.finite(x$mean)
    sd_valid <- function(v) if (sum(is.finite(v)) > 1L) stats::sd(v[is.finite(v)]) else NA_real_
    avg_valid <- function(v) if (any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_
    data.frame(zone_id = x$zone_id[1], zone_label = x$zone_label[1],
      month = x$month[1], variable = x$variable[1],
      mean = avg_valid(x$mean), temporal_sd = sd_valid(x$mean),
      temporal_variance = if (sum(ok) > 1L) stats::var(x$mean[ok]) else NA_real_,
      mean_spatial_sd = avg_valid(x$spatial_sd),
      mean_spatial_variance = avg_valid(x$spatial_variance),
      mean_spatial_cv = avg_valid(x$spatial_cv),
      n_years = sum(ok), n_zone_cells = x$n_zone_cells[1],
      mean_coverage = avg_valid(x$coverage),
      year_start = if (any(ok)) min(x$year[ok]) else NA_integer_,
      year_end = if (any(ok)) max(x$year[ok]) else NA_integer_)
  }))
}

wb_zonal_export <- function(cfg, site, stages, zonal) {
  dest <- file.path(cfg$output_root, "zonal")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  meta <- zonal$meta
  utils::write.csv(meta, file.path(dest, "zonal_metadata.csv"), row.names = FALSE)
  files <- unlist(lapply(stages, function(st) {
    d <- file.path(site$out_dir, paste0(st, "_years"))
    file.path(d, paste0("zonal_", if (st == "baseline") cfg$base_years else cfg$years_annual, ".rds"))
  }))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing zonal year checkpoints: ", paste(missing, collapse = ", "))
  monthly <- do.call(rbind, lapply(files, readRDS))
  if (length(cfg$zonal_monthly))
    utils::write.csv(monthly, file.path(dest, "zonal_monthly.csv"), row.names = FALSE)
  if (isTRUE(cfg$zonal_annual))
    utils::write.csv(wb_zonal_annual(monthly), file.path(dest, "zonal_annual.csv"), row.names = FALSE)
  if (length(cfg$zonal_monthly))
    utils::write.csv(wb_zonal_climatology(monthly), file.path(dest, "zonal_monthly_climatology.csv"), row.names = FALSE)
  writeLines(c("In-memory model summaries, not summaries of saved raster files.",
    "Monthly mean is the spatial mean of each modeled month's cell values.",
    "Annual fluxes are sums of monthly spatial means; states/temperature/radiation are means of monthly spatial means.",
    "Soil and snowpack are end-of-month states; rs and rnet are daily radiation rates.",
    "n_zone_cells is the number of modeled cells assigned to each zone.",
    "n_cells is the number of valid cells; coverage is their ratio.",
    "Monthly spatial_sd and spatial_variance are sample dispersion across valid modeled cells within each zone/month (n-1 denominator).",
    "Monthly spatial_cv is spatial_sd / abs(mean); undefined where zone mean is zero.",
    "zonal_monthly_climatology.csv uses analysis years only: mean is the average monthly zone mean across years.",
    "temporal_sd/temporal_variance in climatology describe interannual variability; mean_spatial_sd and mean_spatial_variance describe spatial variability.",
    "Annual zonal spatial variance is not inferred from monthly variances, because cross-month cell covariance is needed.",
    "Polygons use cell-center assignment; overlapping polygons use last-feature assignment.",
    "Zones are evaluated only on modeled cells. Points use the containing cell.",
    "An incomplete twelve-month series produces NA for its annual value."), file.path(dest, "README.txt"))
  cat("In-memory zonal CSV results:", dest, "\n")
  invisible(dest)
}
