# Zonal min/max update

- Replaces exact within-zone q10/q90 with `spatial_min` and `spatial_max` in zonal monthly output.
- Climatology adds `temporal_min` and `temporal_max` across annual zone means (per calendar month), and `mean_spatial_min`/`mean_spatial_max` averaged across annual spatial extrema.
- Retains spatial/temporal mean, SD, variance, CV and coverage.
- Existing output CSVs are rewritten by a new successful run; old CSV schemas are not migrated.
- No hydrological calculations were changed. R runtime validation pending.
