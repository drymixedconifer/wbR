source_file <- "dev/water_balance_v6_functions.R"
installed_snapshot <- "inst/engine/water_balance_v6_functions.R"
stopifnot(file.exists(source_file))
file.copy(source_file, installed_snapshot, overwrite = TRUE)
message("Updated installed engine snapshot; rebuild or reinstall wbR.")
