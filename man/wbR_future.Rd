\name{wbR_future}
\alias{wbR_future}
\title{Run LOCA2 future climate ensembles}
\usage{wbR_future(config, climate_dir, models, scenario = "ssp245", modes = c("anchored", "direct"), future_periods = list(c(2031L, 2060L), c(2045L, 2074L)), annual_datatype = "INT2S", pr_floor_mm = 1, vpd_floor_kpa = 0.01, engine_path = NULL)}
\description{Private LOCA2 forcing extension using the same water-balance engine as the public wbR package.}
