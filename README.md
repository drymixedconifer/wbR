# wbR

### Spatial water and energy balance modeling in R

**Where water falls, where it goes, and where it falls short.**

`wbR` estimates the climatic water balance of landscapes through time. It combines monthly climate data with elevation, terrain radiation, and soil water-storage capacity to model the movement of water through snowpack and soil, actual evapotranspiration (AET), and climatic water deficit (CWD). The result is a spatially explicit account of how water supply and atmospheric demand interact across complex terrain.

Developed for landscape ecology, forest science, and watershed research, `wbR` supports analyses ranging from individual watersheds to regions of the western United States. Its outputs can be used to explore hydroclimatic change, ecological drought, forest responses, and geographic differences in water availability.

**Authors:** James Johnston (University of Oregon) and Andreas Wion (New Mexico Forest and Watershed Restoration Institute)  
**Code:** [github.com/drymixedconifer/wbR](https://github.com/drymixedconifer/wbR) · **License:** MIT  
**Current package:** v1.1.0.9000 (development). The v6 scientific implementation is under active testing; results should be checked before scientific interpretation.

---

## What wbR calculates

At each modeled grid cell, wbR follows precipitation through rain–snow partitioning, snow accumulation and melt, soil-water recharge and depletion, evaporative demand, realized water use, and excess water. Soil moisture and snowpack carry forward between months and years.

| Product | Interpretation |
|---|---|
| **Climatic water deficit (CWD)** | Atmospheric evaporative demand unmet by available water |
| **Actual evapotranspiration (AET)** | Water returned to the atmosphere through evaporation and transpiration, represented jointly |
| **Potential evapotranspiration (PET)** | Reference atmospheric demand calculated from water–energy conditions |
| **Soil water** | Modeled soil-water storage through time |
| **Snow** | Snow accumulation, storage, and melt, including seasonal summaries |
| **Runoff** | Water exceeding the modeled available storage; *not routed streamflow* |

Additional output variables include precipitation, snowfall, snowmelt, water reaching the soil, temperature, incoming and net radiation, and vapor pressure deficit. Outputs can be written as annual GeoTIFFs, monthly GeoTIFFs, or monthly and annual **zonal tables** for watersheds, ecological classes, and other user-defined geography.

## Why water balance?

Temperature and precipitation alone cannot describe ecological water availability. The same precipitation regime can produce different outcomes on a cool north-facing slope, a sun-exposed ridge, or a deep-soiled valley. `wbR` integrates **water supply, energy demand, terrain, and soil storage** to quantify this heterogeneity rather than treating a climate grid cell as a uniform landscape.

The model is deliberately **vegetation-independent**: it represents climatic and physical controls on water balance rather than predicting separate vegetation-specific transpiration and soil evaporation. This provides a consistent environmental baseline for comparing landscapes and investigating vegetation response.

## Install

A working R installation (R ≥ 4.2) and C++ compilation toolchain are required. The scientific engine uses `terra`, `prism`, `Rcpp`, and `RcppParallel`.

```r
install.packages(c("remotes", "terra", "prism", "Rcpp", "RcppParallel"))
remotes::install_github("drymixedconifer/wbR", upgrade = "never")
library(wbR)

cfg <- wb_config()              # compact, editable configuration
names(cfg)
wb_config(advanced = TRUE)      # scientific and technical controls
```

> The GitHub repository must contain the intended release before installation. The development archive and currently published GitHub version may differ.

## Quick start

Point wbR at an elevation raster, soil available-water-capacity raster, AOI, and a PRISM cache directory. Choose the analysis period and the outputs you want.

```r
library(wbR)

cfg <- wb_config()

# Geography and environmental inputs
cfg$dem_path     <- "data/dem.tif"
cfg$awc_path     <- "data/awc_mm.tif"
cfg$aoi_path     <- "data/watershed.gpkg"
cfg$prism_root   <- "data/PRISM"
cfg$output_root  <- "outputs"

# Consecutive initialization and analysis periods
cfg$base_years   <- 1988:1989
cfg$years_annual <- 1990:2025

# Model resolution and PET formulation
cfg$dem_aggregate_fact <- 1L
cfg$pet_method <- "penman_monteith"
cfg$wind_path  <- "data/wind_u2_monthly.tif"  # or use wind_constant (advanced)

# Annual output rasters
cfg$outputs_to_write <- c("cwd", "aet")
cfg$write_rasters <- TRUE
cfg$monthly_output <- FALSE

wb_run(cfg)
```

`base_years` initialize the state preceding the analysis period; their outputs are **not** written by default. Analysis starts from the final soil and snow state of the preceding baseline year—there is no user-selectable reset to equilibrium.

### Choose the outputs you need

For additional annual products, change the selection and run again:

```r
cfg$outputs_to_write <- c("soil", "snow", "runoff")
wb_run(cfg)
```

Supported annual selections in the current engine are:

```r
c("cwd", "aet", "pet", "ppt", "soil", "pack", "snow",
  "snow_peak", "melt", "runoff", "snowfall", "water")
```

Annual **`snow`** and **`snow_peak`** are seasonal snow metrics, not equivalent to calendar-year precipitation sums. For detailed monthly water and energy variables, use `monthly_output`; see the advanced configuration and engine documentation for available variables and units.

Because very large regional runs can exceed RAM limits when assembling many outputs, wbR supports separate output batches. The development version has a **persistent baseline-end initialization cache** intended to reuse spin-up and baseline state when the scientific inputs have not changed; this pathway still needs end-to-end validation.

### Zonal-only analysis

Compute monthly and annual summaries by watersheds, vegetation types, or other geographic units **without writing annual raster products**:

```r
cfg <- wb_config()

cfg$dem_path     <- "data/dem.tif"
cfg$awc_path     <- "data/awc_mm.tif"
cfg$aoi_path     <- "data/watershed.gpkg"
cfg$prism_root   <- "data/PRISM"
cfg$output_root  <- "outputs_zonal"
cfg$base_years   <- 1988:1989
cfg$years_annual <- 1990:2025
cfg$wind_path    <- "data/wind_u2_monthly.tif"

cfg$write_rasters <- FALSE
cfg$zonal_path <- "data/subwatersheds.gpkg"
cfg$zonal_id_field <- "huc12"  # replace with a field in your vector layer

wb_run(cfg)
```

Zonal inputs can be categorical rasters or supported vector geometries. Zonal tables are collected during model simulation, so monthly GeoTIFFs are not required. `outputs_to_write` and `monthly_years` do not select the variables or years included in zonal CSVs; use `zonal_monthly` and `years_annual` respectively.

## Climate and environmental inputs

| Input | Description | Notes |
|---|---|---|
| **DEM** | Grid, elevation, slope, aspect, latitude, radiation | Verify horizontal CRS and vertical units |
| **AOI** | Study area boundary | Mask and geographic extent |
| **AWC** | Soil available-water storage | Typically supplied in millimeters; confirm data meaning and units |
| **PRISM** | Monthly precipitation and temperature; VPD for Penman–Monteith; atmospheric normals | Historical data may require substantial disk space |
| **Wind** | Twelve-month 2 m wind-speed climatology | Optional constant-speed fallback through `wind_constant` |
| **Zones** (optional) | Watersheds, classifications, polygons, or points | Used for zonal time series |

PRISM VPD is supplied in **hPa** and converted internally to **kPa** for the Penman–Monteith calculations. Do **not** manually modify the original PRISM data to perform that conversion.

### PRISM storage modes

```r
cfg$prism_mode <- "full"    # existing full-grid workflow (default)
# cfg$prism_mode <- "subset" # retain buffered AOI subsets locally
# cfg$prism_mode <- "local"  # existing cache only; no downloads
```

Subset mode is designed to save **persistent disk space** by downloading, extracting, cropping, and retaining buffered rasters. The current official ZIP-based acquisition path may still transfer complete source archives over the network; subset mode does **not** promise reduced download bandwidth. Buffered coverage is needed for interpolation and local lapse-rate fitting.

## Model approach

The model operates monthly on a grid derived from the DEM. Its core processes include:

1. **Spatial climate preparation:** Interpolate PRISM fields and adjust temperature and precipitation for elevation using locally fitted lapse relationships.
2. **Radiation and PET:** Compute terrain-dependent shortwave radiation and longwave terms. Penman–Monteith is the v6 default, using PRISM vapor pressure deficit and 2 m wind speed; Priestley–Taylor and Hamon are alternative methods.
3. **Snow balance:** Partition precipitation between rain and snow, maintain snowpack, and estimate melt using temperature and available radiation.
4. **Soil-water accounting:** Track available storage, AET, unmet demand (CWD), and excess water (runoff), carrying the modeled state between months and years.
5. **Output and aggregation:** Produce annual raster summaries, optional monthly spatial products, and in-simulation zonal tables.

This is a **monthly climatic water-balance model**, not a hydrodynamic streamflow-routing model. Its runoff output represents modeled excess water prior to channel routing, transmission losses, or groundwater processes not explicitly resolved by the model.

## Reuse, reproducibility, and performance

- **Prepared-environment cache:** Reuses terrain and soil preparation when inputs and geometry match.
- **Initialization cache (development):** Saves the state at the end of `base_years` for reuse across output configurations with matching scientific inputs.
- **Resume:** Retains checkpoints to recover interrupted runs when configurations remain compatible.
- **Scientific fingerprinting:** Prevents inappropriate reuse when model code or relevant inputs change.
- **Full-grid computation:** The current engine does **not** automatically tile the spatial calculation; memory scales with the number of modeled cells and requested intermediate/output arrays.

For regional work, initialization and output assembly can require substantial memory and temporary disk space. Start with a small AOI, confirm climate coverage and units, and increase spatial extent only after inspecting performance. Use separate output batches when the complete output set exceeds available memory.

## Current development status and limitations

`wbR` is research software under active development. The v6 integration includes Penman–Monteith PET, PRISM VPD handling, zonal-only execution, and initialization caching. Source-level fixes have been made during integration, but **this development package has not yet passed comprehensive cross-platform checks or end-to-end scientific regression tests**. Do not assume numerical equivalence between v5.5 and v6, or that all optional workflows have been validated.

Important interpretation limits:

- Monthly climate and static atmospheric normals cannot reproduce all event-scale weather, interception, soil, or snow processes.
- Water-balance estimates are vegetation-independent and should not be interpreted as direct measurements of canopy transpiration or vegetation-mediated hydrological change.
- Annual seasonal-snow products span a cold season; interpretation depends on correctly carried prior-autumn state. Missing values must not be treated as zero.
- PRISM coverage, elevation units, AWC definitions, wind assumptions, and domain resolution materially affect results.
- Modeled runoff is **not** equivalent to observed stream discharge.

For reproducible work, record the wbR package version or commit, input data versions, spatial resolution, PET method, years, and complete effective configuration (`wb_config(advanced = TRUE)` is the starting point for exploring available parameters).

## Citation and scientific foundations

If using wbR in research, cite the software repository **and the specific version or commit**. The package includes `CITATION.cff` for citation metadata.

The scientific framing follows the development of climatic water deficit and actual evapotranspiration as ecological measures of water availability, including:

- Stephenson, N. L. (1998). Actual evapotranspiration and deficit: biologically meaningful correlates of vegetation distribution across spatial scales. *Journal of Biogeography*, **25**, 855–870.
- Lutz, J. A., van Wagtendonk, J. W., & Franklin, J. F. (2010). Climatic water deficit, tree species ranges, and climate change in Yosemite National Park. *Journal of Biogeography*, **37**, 936–950. https://doi.org/10.1111/j.1365-2699.2009.02268.x
- McCune, B., & Keon, D. (2002). Equations for potential annual direct incident radiation and heat load. *Journal of Vegetation Science*, **13**, 603–606.

`wbR` implements its **own documented numerical formulations and parameterizations**; citing foundational literature does not imply that every historical equation is reproduced verbatim.

## Contributing and reporting issues

Issues, reproducible test cases, and constructive feedback are welcome through [GitHub Issues](https://github.com/drymixedconifer/wbR/issues). Include the wbR version or commit, operating system, R version, relevant configuration, a minimal reproducible example where possible, and the complete error message or traceback. Do not publish credentials or private filesystem information in issue reports.

---

*wbR · Water balance across landscapes, through time.*

### Zonal spatial quantiles

Zonal monthly CSV output includes `spatial_q10` and `spatial_q90`: the 10th and 90th percentiles of **valid modeled cells** in each zone for a given year, month, and variable (R `quantile()` type 7). Existing means, spatial SD, variance, CV, coverage, and temporal dispersion are unchanged. The monthly climatology includes `mean_spatial_q10` and `mean_spatial_q90`, averages of the corresponding yearly spatial percentiles—not pooled spatial/temporal quantiles. Annual spatial percentiles cannot be derived from monthly marginal quantiles; no such estimate is added to annual tables.

### Derived AET minus CWD

Use `balance` in `outputs_to_write` for an annual raster of **AET − CWD** (mm), or in `zonal_monthly` for monthly cell-level spatial statistics, including q10 and q90. The zonal quantiles are calculated **after subtracting CWD from AET in each modeled cell**; they are not the difference of the separate AET and CWD quantiles. Annual results respect `months_out`. This is a derived diagnostic, not a new water-balance equation. It is not currently a standalone `monthly_output` GeoTIFF option.

```r
cfg$outputs_to_write <- c("cwd", "aet", "balance")
cfg$zonal_monthly <- c("cwd", "aet", "balance")
```
