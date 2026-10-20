# Run from the package root after changing inst/config/defaults.R.
source_text <- readLines("inst/config/defaults.R", warn = FALSE)
output <- c("# Generated from inst/config/defaults.R. Keep synchronized with that source.",
            "# Installed-package configuration does not depend on system.file() lookups.",
            ".wbr_config_defaults <- function(.wb_functions_path) {",
            paste0("  ", source_text), "  cfg", "}")
writeLines(output, "R/config_defaults.R")
