# Zonal balance output fix

Fix `subscript out of bounds` on the first baseline or analysis year when `balance` is selected in `zonal_monthly`.

`wb_year()` returns a numeric matrix with a fixed number of columns. Assigning `dx[, "balance_m01"] <- ...` cannot introduce new columns. The zonal path now creates a separate 12-column matrix calculated per modeled cell as monthly AET minus monthly CWD, then sends that matrix into the unchanged zonal statistics function (including spatial q10/q90).

No model physics or unit changes. No change to the annual raster `balance` path. Runtime tests in R are still required.
