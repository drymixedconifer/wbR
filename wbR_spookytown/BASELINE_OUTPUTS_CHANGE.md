# Baseline output policy

Default: baseline years are simulated to initialize the analysis, but baseline annual rasters and baseline means are not written. Analysis annual rasters/stacks remain. `cfg$write_baseline <- TRUE` opts into baseline annual rasters (seasonal snow requires a valid prior autumn state). `cfg$summarize <- TRUE` opts into period mean rasters. No scientific equations modified. R checks not executed.
