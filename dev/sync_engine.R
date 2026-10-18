# Run from the package root before rebuilding an implementation release.
source_file <- "dev/water_balance_v5.5_functions.R"
installed_snapshot <- "inst/engine/water_balance_v5.5_functions.R"
stopifnot(file.exists(source_file))
file.copy(source_file, installed_snapshot, overwrite = TRUE)
message("Updated installed engine snapshot; rebuild or reinstall wbR.")
