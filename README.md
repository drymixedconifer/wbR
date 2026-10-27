# wbR

### Spatial water and energy balance modeling in R

**Where water falls, where it goes, and where it falls short.**

[Source code](https://github.com/drymixedconifer/wbR) · [Issues](https://github.com/drymixedconifer/wbR/issues) · MIT License

**wbR** is an R package for reconstructing climatic water balance across landscapes and through time. It combines monthly climate data, terrain, radiation, and soil water-storage capacity to estimate **actual evapotranspiration (AET)**, **climatic water deficit (CWD)**, potential evapotranspiration (PET), snowpack dynamics, soil water, and excess water. The same physically informed, vegetation-independent framework can be applied across elevations, ecological communities, watersheds, and regions.

Where temperature and precipitation describe climate inputs, water balance describes their consequences: the water delivered to a landscape, the energy available to evaporate it, what can be stored, and where atmospheric demand remains unmet.

**Authors:** James Johnston (University of Oregon) and Andreas Wion (New Mexico Forest and Watershed Restoration Institute).  
**Software version:** `1.1.0.9000` (development). **Scientific implementation:** v6.1.  
**Status:** Under active testing. The v6.1 snow-cover formulation changes scientific results relative to v6.0; validate outputs before use in publications or management decisions.

---

## Contents

1. [Why water balance?](#why-water-balance)
2. [What wbR estimates](#what-wbr-estimates)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Input data](#input-data)
6. [Model processes and assumptions](#model-processes-and-assumptions)
7. [What's new in v6.1: snow-cover demand](#whats-new-in-v61-snow-cover-demand)
8. [Outputs and their interpretation](#outputs-and-their-interpretation)
9. [Zonal analysis](#zonal-analysis)
10. [PRISM acquisition and local storage](#prism-acquisition-and-local-storage)
11. [Initialization, caching, and reproducibility](#initialization-caching-and-reproducibility)
12. [Configuration reference](#configuration-reference)
13. [Worked analysis patterns](#worked-analysis-patterns)
14. [Performance and troubleshooting](#performance-and-troubleshooting)
15. [Limitations, validation, and scientific use](#limitations-validation-and-scientific-use)
16. [Citing wbR](#citing-wbr)

## Why water balance?

Two sites receiving the same precipitation may experience fundamentally different water availability. Elevation influences temperature and snow, aspect alters incoming radiation, and soils differ in the water they can store and release. An air-temperature trend does not uniquely determine drought stress; neither does precipitation alone. These drivers interact.

wbR represents this interaction at the landscape-cell level, month by month. Two complementary diagnostics are especially useful in ecology:

- **AET**, the amount of water actually evaporated or transpired under modeled water and energy constraints.
- **CWD**, the portion of modeled *effective evaporative demand* that cannot be met by available water.

AET and CWD are not interchangeable. A high-AET site may be warm and wet, while a low-AET site may be cool or severely water-limited. Interpreting them together helps separate water availability from atmospheric demand and facilitates comparisons across distinct hydroclimatic regimes.

The framework is deliberately **vegetation-independent**. It estimates climatic and physical water-balance conditions rather than simulating species-specific stomatal conductance, leaf area, root uptake, canopy interception, or forest demography. That separation is useful when relating a common environmental predictor to different ecological responses.

## What wbR estimates

| Variable | What it represents | Typical units |
|:--|:--|:--|
| `cwd` | Effective evaporative demand left unmet by available water | mm |
| `aet` | Actual evaporation + transpiration represented as a combined flux | mm |
| `pet` | Reference potential evapotranspiration / atmospheric demand | mm |
| `balance` | Derived diagnostic: **AET − CWD** | mm |
| `soil` | Soil water remaining at the end of each month | mm |
| `pack` | Snow water equivalent (SWE) remaining at month end | mm |
| `snow` | Cold-season snow summary across calendar years | mm-based seasonal metric |
| `snow_peak` | Seasonal snowpack maximum | mm SWE |
| `runoff` | Modeled water in excess of soil storage capacity | mm |
| `water` | Rain plus snowmelt delivered to the soil | mm |
| `ppt` / `snowfall` | Downscaled precipitation / snow portion | mm |
| `melt` | Modeled snowmelt | mm |
| `tmean` | Downscaled mean air temperature | °C |
| `vpd` | Vapor pressure deficit used for Penman–Monteith PET | kPa |
| `rs` / `rnet` | Incoming shortwave / net radiation at a snow surface | MJ m⁻² day⁻¹ |
| `scover` | Modeled fraction of the month with snow-cover suppression | 0–1 |

**Availability depends on output type.** The annual output selector is not the same as the monthly or zonal selector; see [Outputs](#outputs-and-their-interpretation). `snow` and `snow_peak` are seasonal metrics, not simply sums of January–December snowfall. `runoff` is *not* routed streamflow.

## Installation

Install R **4.2 or newer** and an R-compatible C++ toolchain. macOS users may need Xcode Command Line Tools. Install runtime dependencies and the source package:

```r
install.packages(c("remotes", "terra", "prism", "Rcpp", "RcppParallel"))
remotes::install_github("drymixedconifer/wbR", upgrade = "never")
library(wbR)

packageVersion("wbR")
wb_engine_info()                  # installed engine source and checksum
names(wb_config())                # compact user-facing settings
names(wb_config(advanced = TRUE)) # full settings
```

**Testing an unpushed ZIP:** unzip the package source, then install that local `wbR` directory using `remotes::install_local("/path/to/wbR", upgrade = "never", force = TRUE)`. The GitHub repository may be behind the development ZIP until changes are tested and pushed.

The compiled C++ model is built through `Rcpp`/`RcppParallel` when required. The first call can be slower because of compilation.

## Quick start

Provide a DEM, available water capacity (AWC) raster, AOI, and a location for PRISM climate data. Choose contiguous baseline and analysis periods, then select outputs.

```r
library(wbR)

cfg <- wb_config()

cfg$dem_path    <- "data/dem.tif"
cfg$awc_path    <- "data/awc_mm.tif"
cfg$aoi_path    <- "data/study_area.gpkg"
cfg$prism_root  <- "data/PRISM"
cfg$output_root <- "outputs"

# Initialization years must immediately precede analysis years
cfg$base_years   <- 1988:1989
cfg$years_annual <- 1990:2025

cfg$dem_aggregate_fact <- 1L
cfg$pet_method <- "penman_monteith"

# Optional wind climatology: 12 monthly layers of 2 m wind speed (m/s)
cfg$wind_path <- "data/wind_u2_monthly.tif"

cfg$outputs_to_write <- c("cwd", "aet")
cfg$write_rasters <- TRUE
cfg$monthly_output <- FALSE

wb_run(cfg)
```

`base_years` are **initialization years**. By default, they do not create baseline output rasters. The annual analysis inherits soil and snow state continuously from the end of the baseline; there is no option to reset the analysis state to spin-up equilibrium.

For a smaller test, shorten the periods while preserving contiguity and ensure all required PRISM files and normals exist.

## Input data

| Input | Required? | Guidance |
|:--|:--|:--|
| DEM (`dem_path`) | Yes | Provides model grid, elevation, slope, aspect, latitude and terrain radiation. Confirm CRS and vertical units. |
| AOI (`aoi_path` or `aoi_bbox`) | Yes | Study domain; can be a supported vector boundary or explicit geographic bounding box. |
| Soil AWC (`awc_path`) | Yes | Maximum available soil-water storage. Confirm whether source units are mm, cm, m or fraction with known depth. |
| PRISM (`prism_root`) | Yes | Monthly `ppt` and `tmean`, required normals, and for Penman–Monteith monthly `vpdmax`/`vpdmin`. |
| Wind (`wind_path`) | No | 12-layer monthly climatology of wind speed at 2 m in m/s; otherwise set `wind_constant` using advanced config. |
| Zones (`zonal_path`) | Optional | Categorical raster, polygons, or supported point geometry used for summaries. |
| Points (`points_path`) | Optional | Independent point sampling/validation workflow. |

The DEM defines spatial support. An elevation aggregation factor changes the cell size; it does **not** increase the native information content of the source climate grids. AWC should have physically meaningful positive storage and valid coverage where the simulation is intended to run.

### PRISM vapor pressure deficit: units matter

PRISM `vpdmin` and `vpdmax` are provided in **hPa**. wbR converts them to **kPa** internally (divide by 10) before Penman–Monteith calculations. **Do not modify the downloaded PRISM archives or divide them yourself** if they are being read by wbR.

The v6 Penman–Monteith formulation uses monthly VPD in conjunction with a wind climatology or constant wind. Shortwave transmittance and dew-point normals supply other radiation-related inputs. This is a monthly approximation to reference evapotranspiration, not a daily meteorological simulation.

## Model processes and assumptions

At each modeled cell, wbR advances a monthly water- and energy-balance state.

### 1. Climate downscaling

PRISM precipitation and temperature are interpolated to modeled cells and adjusted for elevation differences using fitted regional temperature and log-precipitation lapse relationships. Relevant fitting windows and coefficient limits are advanced parameters. Broad climate-grid resolution still constrains the representation of subgrid weather variability.

### 2. Solar radiation and evaporative demand

Monthly clear-sky shortwave radiation depends on latitude, elevation, slope and aspect. The terrain calculation integrates daily direct-beam geometry and includes diffuse and terrain-reflected radiation; it does **not** calculate horizon shading cast by distant terrain. PRISM transmittance normals adjust shortwave radiation, while mean temperature and dew-point normals contribute to the longwave term.

Available PET methods:

- `"penman_monteith"` (**default**): radiation and aerodynamic terms, with PRISM VPD and 2 m wind.
- `"priestley_taylor"`: radiation-based PET using a configurable coefficient.
- `"hamon"`: temperature/daylength formulation with heat-load adjustment.

These approaches can produce materially different estimates; report the method in manuscripts.

### 3. Rain, snow, and melt

Monthly precipitation is partitioned into rain and snow using mean-temperature thresholds, with mixed precipitation in the transition range. Snowpack persists across months and years. The default melt method combines a temperature term with net radiation at the snow surface; an alternative degree-day method is retained for comparison.

### 4. Snow-cover suppression of soil demand — **new in v6.1**

The 6.1 formulation estimates the fraction of a month affected by snow cover from available SWE and potential melt. Where modeled snow cover is present, potential **soil/plant** evaporative demand is reduced. The reference atmospheric PET output remains unchanged.

For monthly snow-cover fraction \(f_{snow}\), the effective demand used for soil-water accounting is:

\[
PET_{soil} = (1 - f_{snow})\,PET
\]

Monthly CWD is the unmet portion of **effective soil demand**, not necessarily of the full reported atmospheric PET. Consequently, in v6.1 **`AET + CWD` can be lower than reported `PET`** during snow-covered months. This is intentional, not an accounting error.

The operational model uses potential melt and available snowpack to approximate monthly snow-covered fraction. It is **not** based on satellite-observed snow-covered area and is not a subdaily energy-balance model.

Defaults:

```r
cfg <- wb_config(advanced = TRUE)
cfg$snow_cover_demand <- "exclude"  # v6.1 default
cfg$snow_cover_min_mm <- 10         # SWE threshold (mm)
```

The legacy setting is available for **controlled comparisons**:

```r
cfg$snow_cover_demand <- "soil"     # previous v6 treatment
```

Changing either snow parameter can change AET, CWD, soil storage, runoff, and subsequent-year states. **v6.1 results should not be assumed numerically identical to v6.0.** The snow-cover calculation is a new scientific assumption requiring empirical evaluation.

### 5. Monthly soil-water accounting

Rain and snowmelt provide incoming liquid water \(W\). If \(W\) meets or exceeds effective potential demand, the model satisfies AET, recharges soil water, and routes storage in excess of AWC to `runoff`. If demand exceeds incoming water, the model withdraws water from soil according to its storage-depletion equation and assigns unmet effective demand to CWD. Soil water and SWE carry forward continuously.

In this model, `runoff` means **excess from the modeled storage accounting**. It is not a prediction of channel discharge, groundwater recharge, stream connectivity, or water routed to neighboring cells.

### 6. What `balance` means

`balance = AET - CWD` is a convenient derived contrast (mm). It is **not** the conservation-of-mass water-balance residual, and it is not a dimensionless AET/CWD ratio. Positive and negative values should be interpreted alongside their constituent AET and CWD terms.

## Outputs and their interpretation

### Annual raster outputs

The current source validates these `outputs_to_write` names:

```r
c("cwd", "aet", "pet", "ppt", "soil", "pack", "snow",
  "snow_peak", "melt", "runoff", "snowfall", "water", "balance")
```

For example:

```r
cfg$outputs_to_write <- c("cwd", "aet", "balance")
cfg$write_rasters <- TRUE
cfg$monthly_output <- FALSE
wb_run(cfg)
```

`months_out` determines which months contribute to selected annual flux sums. Soil storage and snowpack are **state variables**, so their seasonal aggregation and interpretation differ from sums of fluxes. `snow` and `snow_peak` summarize a cold season crossing calendar years; interpret these alongside the model's seasonal conventions.

An all-`NA` band may be retained rather than misrepresented as zero. An `NA` output is **missing or undefined**, not zero snow or zero deficit. Inspect run warnings before interpreting any such band.

### Monthly raster outputs

Monthly raster selections use `monthly_output` (and advanced `monthly_years` for the requested years), independently of `outputs_to_write`:

```r
cfg <- wb_config(advanced = TRUE)
cfg$monthly_output <- c("cwd", "aet", "soil", "pack", "runoff", "scover")
cfg$monthly_years <- list(
  cwd = 2000:2005,
  aet = 2000:2005,
  soil = 2000:2005,
  pack = 2000:2005,
  runoff = 2000:2005,
  scover = 2000:2005
)
```

Other engine-supported monthly variables include `pet`, `water`, `ppt`, `snowfall`, `tmean`, `rs`, `rnet`, `melt`, and `vpd`. `scover` is a unitless 0–1 monthly fraction. `balance` is supported for annual rasters and monthly **zonal** tables, but not as an independent monthly GeoTIFF selector in this build.

Monthly raster files can be large. Disable them when you only need annual rasters or zonal CSVs:

```r
cfg$monthly_output <- FALSE
```

### Output switches are separate

| Setting | Controls |
|:--|:--|
| `write_rasters` | Whether annual spatial rasters are written |
| `outputs_to_write` | Which annual raster variables are written, if enabled |
| `monthly_output` | Which monthly GeoTIFF variables are written |
| `monthly_years` (advanced) | Years scheduled for those monthly raster products |
| `zonal_path` | Whether in-simulation zonal statistics are produced |
| `zonal_monthly` | Variables summarized in zonal tables |
| `summarize` | Optional period summary raster generation; off by default |
| `write_baseline` | Optional baseline output rasters; off by default |

`write_rasters = FALSE` does **not** prohibit zonal output. If no raster products are requested, specify a zonal or point output so the run has something to return.

## Zonal analysis

One of wbR's central features is the ability to summarize simulated water-balance values by **arbitrary user-defined geographic units**, without saving the corresponding monthly raster stacks.

Zones can represent catchments, ecoregions, vegetation classes, ownership units, elevation classes, or other valid raster/vector classifications. The statistical machinery does not depend on HUC conventions.

```r
cfg <- wb_config()

cfg$dem_path    <- "data/dem.tif"
cfg$awc_path    <- "data/awc_mm.tif"
cfg$aoi_path    <- "data/study_area.gpkg"
cfg$prism_root  <- "data/PRISM"
cfg$output_root <- "outputs_zonal"

cfg$base_years   <- 1988:1989
cfg$years_annual <- 1990:2025
cfg$wind_path <- "data/wind_u2_monthly.tif"

cfg$write_rasters <- FALSE
cfg$monthly_output <- FALSE
cfg$zonal_path <- "data/zones.gpkg"
cfg$zonal_id_field <- "zone_id"
cfg$zonal_monthly <- c("cwd", "aet", "balance", "soil", "pack", "scover")

wb_run(cfg)
```

The current implementation computes monthly zonal values directly from modeled cells. It uses a **persistent cell-to-zone mapping cache** to avoid repeating the same polygon-to-grid assignment across compatible runs. Polygon assignment is based on cell centers. **Overlapping polygons use last-feature assignment**, not fractional allocation to multiple overlapping zones; this matters for users analyzing overlapping treatments or fires.

### Zonal output files

Under `<output_root>/zonal/`:

| File | Meaning |
|:--|:--|
| `zonal_metadata.csv` | Zone IDs, labels, and modeled cell counts |
| `zonal_monthly.csv` | Per-zone, per-year, per-month, per-variable spatial statistics |
| `zonal_annual.csv` | Aggregated annual zonal time series |
| `zonal_monthly_climatology.csv` | Calendar-month climatologies from analysis years |
| `README.txt` | Definitions and aggregation rules |

The monthly table records `mean`, `spatial_min`, `spatial_max`, `spatial_sd`, `spatial_variance`, `spatial_cv`, `n_cells`, `n_zone_cells`, and `coverage`.

The climatology table distinguishes two kinds of variability:

- **Temporal envelope:** `temporal_min` and `temporal_max` are the extrema of annual *zone means* for each calendar month. `temporal_sd` and `temporal_variance` are interannual statistics.
- **Spatial envelope:** `mean_spatial_min` and `mean_spatial_max` average the per-year spatial extrema. They are **not** the absolute minimum and maximum across all cells and all years. `mean_spatial_sd` describes typical within-zone spatial dispersion.

Exact q10/q90 calculation was removed in favor of minimum and maximum; older README references to spatial quantiles are obsolete for this build.

Annual zonal fluxes are assembled from sums of monthly zone means; monthly state or rate variables use annual means. Annual spatial variances are **not** reconstructed from monthly variances because the covariance between months would be required.

### Plot a monthly climatology ribbon

This example uses the temporal range among annual zone means:

```r
library(ggplot2)
clim <- read.csv("outputs_zonal/zonal/zonal_monthly_climatology.csv")
clim <- subset(clim, variable == "cwd" & zone_id == unique(zone_id)[1])

ggplot(clim, aes(month, mean)) +
  geom_ribbon(aes(ymin = temporal_min, ymax = temporal_max), alpha = 0.18) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = 1:12) +
  labs(x = "Month", y = "CWD (mm)", title = "Monthly climatic water deficit") +
  theme_minimal()
```

For a **spatial** envelope instead, use `mean_spatial_min` and `mean_spatial_max` in `geom_ribbon()`. Spatial extrema are sensitive to individual extreme cells; neither ribbon is a confidence interval.

## PRISM acquisition and local storage

Climate data can be substantial in regional historical applications. wbR offers three modes:

```r
cfg$prism_mode <- "full"    # default: keep the usual complete local archive
cfg$prism_mode <- "subset"  # retain buffered AOI crops of acquired PRISM grids
cfg$prism_mode <- "local"   # read existing local inputs without downloads
```

For `"subset"`, wbR downloads the required PRISM material, extracts source rasters, crops them around a buffered AOI, and retains the smaller subset. **This saves retained disk space, but the ZIP download may still transfer the full source grid.** It is not yet a guaranteed remote range-request / bandwidth-saving download method.

```r
cfg$prism_mode <- "subset"
cfg$prism_buffer_km <- 25
```

The buffer should be large enough to support climate interpolation and lapse-rate estimation. For reproducibility, preserve the PRISM product names, resolution, dates, and the normals period used. A small example AOI can still need a wider climate neighborhood.

## Initialization, caching, and reproducibility

The model distinguishes three computational stages:

1. **Prepared environment:** DEM and AWC alignment, terrain derivatives, latitude, radiation and other reusable environmental inputs.
2. **Spin-up and baseline initialization:** establish and carry forward soil-water and snowpack states; save the final state at the end of `base_years`.
3. **Annual analysis:** begin at the baseline-end state and calculate requested years and products.

A development-stage **persistent initialization cache** is intended to avoid repeating stages 1–2 when only output variables, zones, or raster-writing choices change. Its scientific fingerprint should protect against incompatible reuse when climate, model code, parameters, grid or other scientifically relevant inputs differ.

```r
cfg$outputs_to_write <- c("cwd", "aet")
wb_run(cfg)

cfg$outputs_to_write <- c("soil", "snow", "runoff")
wb_run(cfg)   # should reuse compatible initialization; analysis still recalculates
```

Because this behavior is development-stage, **verify the cache-hit message and output integrity** before assuming initialization was reused. The second run may still recompute annual water balance; cache reuse is not the same as reusing a full saved historical simulation.

Zone mapping is cached separately from scientific initialization, so changing the output selection does not inherently require new spatial overlay. Different zone boundaries, IDs or grid alignment appropriately require remapping.

`resume = TRUE` is intended to recover supported intermediate computations after interruption. Cache validity, fingerprinting, and the integrity of checkpoint files should be checked when debugging changed code or inputs.

## Configuration reference

`wb_config()` exposes a compact set of commonly used settings, while `wb_config(advanced = TRUE)` returns the complete list of scientific and runtime parameters.

| Setting | Purpose |
|:--|:--|
| `dem_path`, `awc_path`, `aoi_path`, `aoi_bbox` | Spatial inputs and modeling domain |
| `prism_root`, `prism_mode`, `prism_buffer_km` | Climate storage and acquisition |
| `output_root` | Destination for saved output products |
| `base_years`, `years_annual`, `months_out` | Initialization period, analysis period, annual-month selection |
| `dem_aggregate_fact` | Aggregate DEM-derived resolution |
| `pet_method`, `wind_path` | PET selection and wind input |
| `outputs_to_write`, `write_rasters` | Annual raster output selection and switch |
| `monthly_output` | Monthly raster output selection |
| `zonal_path`, `zonal_id_field`, `zonal_monthly` | Zone geometry, labels and requested variables |
| `points_path` | Independent point workflow |
| `summarize`, `write_baseline` | Optional summary and baseline rasters |
| `wb_cores`, `resume` | Parallelism and restart behavior |

Advanced scientific controls include lapse fitting limits and windows, radiation parameters, wind constant and aerodynamic scaling, snow/rain thresholds, melt coefficients, soil and elevation units, spin-up tolerances, memory settings, output storage, `snow_cover_demand`, and `snow_cover_min_mm`.

```r
cfg <- wb_config(advanced = TRUE)
cfg$wind_path <- NULL
cfg$wind_constant <- 2             # m/s, when wind climatology is unavailable
cfg$snow_cover_demand <- "exclude"
cfg$snow_cover_min_mm <- 10
cfg$monthly_output <- FALSE
```

**Important:** `analysis_initialization` was removed. Analysis years **always** follow the final baseline state, and contiguous time periods are required. `base_years` are not a user-selectable alternative equilibrium reset.

## Worked analysis patterns

### Batch outputs on memory-limited machines

```r
cfg <- wb_config()
# ... supply required inputs, geography, periods and wind ...
cfg$monthly_output <- FALSE
cfg$outputs_to_write <- c("cwd", "aet")
wb_run(cfg)

cfg$outputs_to_write <- c("soil", "snow", "runoff")
wb_run(cfg)
```

Use this to avoid retaining many large output arrays simultaneously. Do not assume different output selections reduce the core per-cell simulation memory itself.

### Zonal time series without rasters

```r
cfg$write_rasters <- FALSE
cfg$monthly_output <- FALSE
cfg$zonal_path <- "data/ecoregions.gpkg"
cfg$zonal_id_field <- "ecoregion_id"
cfg$zonal_monthly <- c("cwd", "aet", "pet", "soil", "pack", "scover")
wb_run(cfg)
```

### Compare snow-cover formulations

```r
cfg_new <- wb_config(advanced = TRUE)
# ... set identical inputs, years, and separate output roots ...
cfg_new$snow_cover_demand <- "exclude"

cfg_old <- cfg_new
cfg_old$snow_cover_demand <- "soil"
cfg_old$output_root <- "outputs_legacy_snow_demand"

# Run both after supplying actual paths and required data.
# wb_run(cfg_new)
# wb_run(cfg_old)
```

A controlled comparison should verify that each configuration uses its own appropriate state cache and outputs. Differences in snow-covered seasons are expected.

## Performance and troubleshooting

wbR's full-grid engine is computationally efficient during C++ annual simulation, but environmental preparation, large raster matrices, and output assembly can become RAM or disk bottlenecks at regional scale. **This build does not perform automatic spatial tiling.** Resolution and domain area determine the number of modeled cells.

Practical steps:

- Start with a small AOI and validate input coverage, CRS and units before western-US-scale runs.
- At large extent, use coarser justified DEM aggregation, fewer monthly output variables, or separate annual-output batches.
- For zonal-only analysis, set `write_rasters = FALSE` and `monthly_output = FALSE`; monthly zonal values are still computed from model cells.
- Use `wb_config(advanced = TRUE)` when tuning `monthly_chunk`, memory limits, wind constants or scientific coefficients.
- Check the run report, cached files, and first-year values after package updates; a successful file write is not proof of scientifically valid values.

| Symptom | Inspect first |
|:--|:--|
| No valid DEM elevations / invalid latitude transform | DEM cell values, CRS, mask, transform, current engine revision; a projected grid may include invalid cells outside the study domain |
| All-`NA` model results | DEM/AWC coverage, valid modeled-cell set, climate-grid interpolation and temperature/precipitation coverage |
| VPD range rejected | Verify PRISM raw VPD is hPa, and internal conversion to kPa occurs before validation |
| Snow baseline/first-year missing | Prior-autumn state availability and correct seasonal initialization; `NA` is not zero |
| Run uses too much RAM | Modeled cell count, DEM aggregation, annual variables, monthly arrays, zonal selections |
| Long startup | Environmental preparation, terrain radiation, latitude processing, PRISM normals and C++ compilation |
| Zonal CSV missing | `zonal_path`, identifier field, overlap with modeled cells, and output-directory permissions |
| PRISM download failure | Archive completeness, file names, internet connection, API changes, remaining disk space |

When reporting errors, include `sessionInfo()`, package version, `wb_engine_info()`, the failing call, traceback, key configuration choices, modeled domain/resolution, and the run report. Redact private paths or credentials as necessary.

## Limitations, validation, and scientific use

This is a **monthly, vegetation-independent, spatially explicit reference water-balance model**. Its results are estimates conditional on the model structure and inputs, not direct observations of ecosystem water use or drought injury.

In particular:

- **Vegetation:** No explicit vegetation-specific canopy interception, transpiration physiology, dynamic rooting depth, forest mortality, or species competition.
- **Atmospheric forcing:** Monthly means and static atmospheric normals cannot reproduce event-scale humidity, rainfall timing, wind, or snow processes.
- **Snow:** Snow-cover fraction is inferred from SWE and melt potential, not directly measured from snow maps. The v6.1 demand-suppression formulation is a methodological choice that should be tested against independent observations.
- **Terrain:** Slope/aspect radiation is modeled; distant-terrain horizon shading is not.
- **Soils:** AWC represents available storage, not a mechanistic multilayer soil column with full vertical flow and capillary exchange.
- **Hydrology:** No overland/channel routing, explicit lateral water movement, or calibrated streamflow prediction.
- **Zonal geometry:** Cell-center assignment and overlapping-polygon conventions influence boundary results, especially at coarse resolution.
- **Scale:** Climate-input resolution and soil uncertainty do not disappear when the simulation is run at finer DEM resolution.

**Validation status:** The integrated v6.1 development archive has undergone source-level assembly and iterative user testing, but there is no claim here of comprehensive R package checking, independent field validation, or numerical equivalence across all PET/snow methods and output combinations. Record the exact package version, Git commit, configuration, data provenance and resolution when publishing results.

## Citing wbR

Please cite **wbR**, its version or specific Git commit, and the scientific foundations relevant to your analysis. Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

**Primary references**

- Stephenson, N. L. (1998). Actual evapotranspiration and deficit: biologically meaningful correlates of vegetation distribution across spatial scales. *Journal of Biogeography*, **25**, 855–870.
- Stephenson, N. L. (1990). Climatic control of vegetation distribution: the role of the water balance. *The American Naturalist*, **135**, 649–670.
- Lutz, J. A., van Wagtendonk, J. W., & Franklin, J. F. (2010). Climatic water deficit, tree species ranges, and climate change in Yosemite National Park. *Journal of Biogeography*, **37**, 936–950. https://doi.org/10.1111/j.1365-2699.2009.02268.x
- McCune, B., & Keon, D. (2002). Equations for potential annual direct incident radiation and heat load. *Journal of Vegetation Science*, **13**, 603–606.
- Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). *Crop evapotranspiration: Guidelines for computing crop water requirements*. FAO Irrigation and Drainage Paper 56.

These references provide scientific context. The particular wbR algorithms and coefficient choices should be described from the **actual version used**; they should not be assumed identical to every equation in the cited literature.

## Contributing

Bug reports, reproducible benchmarks, methodological questions, and contributions are welcome through [GitHub Issues](https://github.com/drymixedconifer/wbR/issues). If reporting an unexpected result, include a small representative AOI and input metadata when redistribution permits.

---

*wbR — landscape water balance, from environmental gradients to changing climates.*
