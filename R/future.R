# Private future-climate interface. Same engine file as public wbR.
# This file is intentionally absent from the public package candidate.
#' Run a LOCA2 future climate ensemble with the shared wbR v6.3 engine
#' @param config A wb_config(advanced = TRUE) list, or a compact wb_config() list.
#' @param climate_dir Directory containing monthly LOCA2 NetCDF files.
#' @param models Model names matching the LOCA2 NetCDF filenames.
#' @param scenario Scenario suffix, e.g. ssp245.
#' @param modes anchored, direct, or both.
#' @param future_periods List of c(start,end) year ranges.
#' @param annual_datatype INT2S or FLT4S for final annual stacks.
#' @param pr_floor_mm Ratio stabilization for anchored precipitation.
#' @param vpd_floor_kpa Ratio stabilization for anchored VPD.
#' @param engine_path Optional alternate v6.3 engine source for development.
#' @export
wbR_future <- function(config, climate_dir, models, scenario = "ssp245",
                       modes = c("anchored", "direct"),
                       future_periods = list(c(2031L, 2060L), c(2045L, 2074L)),
                       annual_datatype = "INT2S",
                       pr_floor_mm = 1, vpd_floor_kpa = 0.01, engine_path = NULL) {
  if (!is.list(config)) stop("config must be a list from wb_config().")
  path <- .wbr_engine_path(engine_path)
  defaults <- .wbr_config_defaults(path)
  unknown <- setdiff(names(config), names(defaults))
  if (length(unknown)) stop("Unknown wbR configuration setting(s): ", paste(unknown, collapse = ", "))
  defaults[names(config)] <- config
  # Source the shared scientific engine into an isolated environment. The future
  # extension rebinds only wb_climate_year() in that environment, not the package namespace.
  env <- new.env(parent = asNamespace("wbR"))
  sys.source(path, envir = env)
  future_path <- system.file("engine", "future_loca2_functions.R", package = "wbR")
  if (!nzchar(future_path)) stop("Private LOCA2 extension not installed.")
  sys.source(future_path, envir = env)
  defaults$code_md5 <- unname(tools::md5sum(path))
  fut <- list(climate_dir = climate_dir, models = models, scenario = scenario, modes = modes,
              future_periods = future_periods, annual_datatype = annual_datatype,
              pr_floor_mm = pr_floor_mm, vpd_floor_kpa = vpd_floor_kpa,
              code_md5 = unname(tools::md5sum(future_path)))
  env$fut_main(defaults, fut)
}
