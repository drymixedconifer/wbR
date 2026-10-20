# Monthly-years validation fix

`monthly_output = FALSE` or `write_rasters = FALSE` bypasses monthly raster year-list validation; zonal outputs are unaffected. The original monthly-year validation remains in effect whenever monthly rasters are enabled.
