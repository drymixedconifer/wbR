# Blockwise latitude fix (development)

Replaced the DEM + cell-index raster stack with a disk-backed blockwise WGS84 transform. Only valid DEM cells are transformed. Added early modeled-cell latitude validation and explicit all-NA spin-up diagnostic. No water-balance equations changed.

**Not R-runtime tested.** Rebuild of prepared-environment cache is required to exercise updated latitude calculations; stale cached environments should not be reused. Check CRS, PRISM coverage, AWC if the all-NA check fires.
