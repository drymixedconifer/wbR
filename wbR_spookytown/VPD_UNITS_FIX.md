# PRISM VPD units fix

- Source: PRISM `vpdmax`, `vpdmin` in hPa.
- Internal Penman–Monteith input: monthly mean VPD in kPa.
- Conversion: `(vpdmax + vpdmin) / 20`.
- Applies to normals, spin-up, and historical time series through the shared `wb_vpd_stack()` function.
- Raw PRISM cache remains unchanged.
- Both installed engine and development engine source are updated.

Validation still needs an actual R smoke test; R is not available in the build environment.
