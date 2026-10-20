# Configuration cleanup — development build

- `wb_config()` now exposes frequently changed inputs, time periods and output decisions only.
- `wb_config(advanced = TRUE)` retains the full controls for scientific sensitivity and tuning.
- `wb_run()` expands either configuration to the complete engine configuration without changing existing equations.
- `analysis_initialization` is rejected, and analysis always receives the final baseline state. The periods must be consecutive.
- Baseline rasters and summary means are opt-in as before.
- Unknown cfg fields are rejected to catch misspellings.
- Existing full configuration lists still work unless they include removed `analysis_initialization`.

**Not numerical regression-tested**: R is unavailable in this build environment.
