\name{wb_run}
\alias{wb_run}
\title{Run the v5.5 water-balance engine}
\usage{wb_run(config, engine_path = NULL)}
\arguments{\item{config}{A named configuration list returned by wb_config.}\item{engine_path}{Optional editable engine source path.}}
\value{The runtime configuration, invisibly.}
\description{Runs the original v5.5 in-memory engine, with optional monthly zonal CSV summaries.}
