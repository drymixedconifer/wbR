# Configuration lookup fix

`wb_config()` now evaluates the package defaults from the R namespace rather than locating `inst/config/defaults.R` through `system.file()`. The original `inst/config/defaults.R` is retained as the editable development reference. The generated `R/config_defaults.R` mirrors its expressions.

If the installed package still shows the old error after reinstalling, restart R and check `find.package("wbR")` and `packageVersion("wbR")`, as an older installation may be getting loaded.

The VPD hPa → kPa fix is preserved. No scientific equations changed.
