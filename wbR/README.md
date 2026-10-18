# wbR 1.0.0

**Spatial water and energy balance modeling in R**

wbR is a landscape-scale water-balance model for examining geographic and temporal variation in water availability and evaporative demand. It combines monthly climate inputs with topography, solar radiation, and soil water-storage capacity to estimate actual evapotranspiration (AET), climatic water deficit (CWD), potential evapotranspiration (PET), snow accumulation and melt, soil-water storage, and excess water (runoff). The model is designed for spatially explicit ecological and watershed analyses, from individual study areas to larger regions where computational resources permit.

**Authors:** James Johnston (University of Oregon) and Andreas Wion (New Mexico Forest and Watershed Restoration Institute)  
**Source:** [github.com/drymixedconifer/wbR](https://github.com/drymixedconifer/wbR) · **License:** MIT  
**Release:** wbR 1.0.0 — initial public source release. The core uses a full-grid model without automatic spatial tiling. Validate model outputs for your study area before scientific interpretation.

## Contents

1. [Installation](#installation)
2. [Inputs and spatial domain](#inputs-and-spatial-domain)
3. [Run a model](#run-a-model)
4. [How the model works](#how-the-model-works)
5. [Configuration reference](#configuration-reference)
6. [Raster and tabular outputs](#raster-and-tabular-outputs)
7. [Seasonal snow and soil products](#seasonal-snow-and-soil-products)
8. [Zonal time series](#zonal-time-series)
9. [Prepared environment and restarts](#prepared-environment-and-restarts)
10. [Computational requirements](#computational-requirements)
11. [Interpretation and limitations](#interpretation-and-limitations)
12. [Scientific references](#scientific-references)

## Installation

Install the package from a local checkout after installing its R dependencies. A working C++ toolchain is required to compile the embedded scientific engine on first use.

```r
install.packages(c("terra", "prism", "Rcpp", "RcppParallel"))
install.packages("/path/to/wbR", repos = NULL, type = "source")
library(wbR)

cfg <- wb_config()   # all installed defaults, as an editable named list
names(cfg)           # available settings
wb_engine_info()     # installed source and checksum
```

The installed package supplies `wb_config()`, `wb_run()`, and `wb_engine_info()`. Configuration does not execute the simulation; `wb_run(cfg)` does.

## Inputs and spatial domain

| Input | Purpose | Guidance |
|---|---|---|
| Digital elevation model (`dem_path`) | Model grid, elevation, slope, aspect and latitude | Confirm horizontal CRS and elevation units. Resolution follows the DEM and aggregation factor. |
| Area of interest (`aoi_path` or `aoi_bbox`) | Geographic extent / mask | Vector geography or configured longitude–latitude bounding box. |
| Available water capacity (`awc_path`) | Maximum modeled soil-water storage | Check units (`awc_units`) and valid-data mask. |
| PRISM monthly data (`prism_root`) | Precipitation and temperature; configured atmospheric normals | Local cache is reused; additional downloads may be needed. |
| Optional zones (`zonal_path`) | Spatial aggregation of monthly model results | One-layer categorical raster, polygon features, or points. |

**Missing-data handling:** Raster NoData codes such as `-9999` must be identified correctly in the inputs. Do not interpret them as real elevation or soil capacity. Set `dem_z_units` explicitly when the DEM's units are uncertain; automatic inference uses an elevation-range heuristic.

The model uses a projected metric working grid (`metric_crs`, normally `EPSG:5070`). `dem_aggregate_fact = 1L` retains the native input DEM resolution; `3L` aggregates 3 × 3 source cells. It does not directly specify a resolution in meters.

## Run a model

Use `wb_config()` to obtain the full set of installed defaults. Change only parameters needed for your domain and requested outputs. `wb_run(cfg)` executes the run.

The following is a complete **example of the principal settings**. Replace paths, years, and identifiers with those for your study area. All other parameters remain available through `wb_config()`.

```r
library(wbR)

cfg <- wb_config()

# Inputs
cfg$dem_path <- "/project/terrain/dem.tif"
cfg$aoi_path <- "/project/boundaries/study_area.gpkg"
cfg$awc_path <- "/project/soils/awc_mm.tif"
cfg$prism_root <- "/project/climate/PRISM"
cfg$output_root <- "/project/results"

# Grid and simulation years
cfg$dem_aggregate_fact <- 1L
cfg$dem_z_units <- "meters"  # or "feet"; "auto" is available
cfg$base_years <- 1988:1989
cfg$years_annual <- 1990:2025
cfg$months_out <- 1:12

# Persistent, reusable terrain preparation
cfg$prepared_environment <- TRUE
cfg$environment_cache <- "/project/cache/prepared_environment"

# Annual rasters and mean-of-period summary maps
cfg$write_rasters <- TRUE
cfg$outputs_to_write <- c("aet", "cwd", "pet", "ppt", "snow", "soil", "runoff")
cfg$summarize <- TRUE

# Monthly rasters are optional and independent of the zonal histories
cfg$monthly_output <- FALSE
cfg$monthly_years <- NULL

# Optional in-simulation watershed time series
cfg$zonal_path <- "/project/boundaries/HUC12.gpkg"
cfg$zonal_id_field <- "HUC12"
cfg$zonal_monthly <- c("aet", "cwd", "pet", "ppt", "soil", "pack", "runoff", "melt")
cfg$zonal_annual <- TRUE

# Retain compatible restart checkpoints after an interruption
cfg$resume <- TRUE

wb_run(cfg)
```

Set `cfg$zonal_path <- NULL` to disable zonal processing. For a categorical raster of zones, set `cfg$zonal_id_field <- NULL`. To enable monthly TIFFs instead, set `monthly_output` to a character vector of variables and supply corresponding `monthly_years` values. **If `monthly_output <- FALSE`, use `monthly_years <- NULL`.**

## How the model works

wbR follows the sequence below. Soil-water and snowpack states carry between months and years; the model is not a series of independent annual calculations.

1. **Validate settings and discover inputs.** Check paths, periods, raster products and available memory; locate cached PRISM data and download missing inputs if permitted.
2. **Prepare environmental geography.** Align the DEM, area of interest and available water capacity; derive slope, aspect, latitude, heat load and related terrain products. Reuse a matching prepared-environment cache when available.
3. **Adjust climate to terrain.** Use spatially varying elevation relationships fitted from PRISM climate and topography. Local estimation is governed by the lapse-window and coefficient-limit parameters. The full-grid core does not automatically substitute a single regional lapse rate.
4. **Calculate radiation and PET.** Monthly terrain-aware radiation uses solar geometry, slope, aspect, elevation and configured cloud-transmittance/humidity normals. Priestley–Taylor is the default PET formulation; Hamon is optional.
5. **Spin up initial states.** Repeated climatological months allow soil water and snowpack to approach a stable starting condition subject to configured convergence tolerances.
6. **Simulate baseline and analysis years.** Each calendar month partitions precipitation into rain and snow, calculates melt and water reaching the soil, estimates potential and actual evapotranspiration, updates soil storage, calculates deficit and excess water, and carries the remaining states into the next month.
7. **Collect outputs during execution.** Requested annual rasters and monthly rasters are written as configured. Monthly zonal means are calculated from model arrays, without requiring monthly TIFFs.
8. **Assemble period products.** The analysis-year raster stacks and optional mean-of-period rasters are saved; compatible checkpoints enable restarting interrupted runs.

The baseline years are **simulation/initialization years**; they are not necessarily the most appropriate scientific climatology. Analysis-period mean maps can be generated independently using `summarize = TRUE`.

## Configuration reference

| Setting | Description |
|---|---|
| `dem_path`, `aoi_path`, `awc_path`, `prism_root`, `output_root` | Required input and destination paths |
| `dem_aggregate_fact`, `dem_z_units`, `metric_crs` | Model grid and DEM interpretation |
| `base_years`, `years_annual`, `months_out`, `analysis_initialization` | Time periods and starting states |
| `lapse_window_km`, `lapse_t_limits`, `lapse_p_limits` | Spatially varying climate lapse-rate estimation |
| `pet_method`, `pt_alpha`, `cloud_transmittance`, `humidity_source` | PET and atmospheric controls |
| `snow_t_all`, `rain_t_all`, `melt_method`, `melt_mm_per_deg_day`, `melt_mm_per_mj`, `snow_albedo` | Precipitation phase and snowmelt |
| `spinup_max_iter`, `soil_tol_mm`, `pack_tol_mm`, `spinup_action` | Spin-up convergence |
| `outputs_to_write`, `summarize`, `write_rasters` | Annual raster products and reporting-period means |
| `monthly_output`, `monthly_years`, `monthly_chunk`, `monthly_datatype` | Monthly TIFFs and their memory/disk tradeoffs |
| `zonal_path`, `zonal_id_field`, `zonal_monthly`, `zonal_annual` | In-simulation zonal CSVs |
| `prepared_environment`, `environment_cache`, `resume` | Reuse of terrain data and restart behavior |
| `memory_limit_gb`, `terra_memfrac`, `wb_cores` | Computational resources |

To view **every parameter and its current default**, run `str(wb_config())`. The canonical full list is in `inst/config/defaults.R`.

## Raster and tabular outputs

| Annual output selector | Interpretation | Units |
|---|---|---|
| `cwd` | Sum of unmet potential evapotranspiration (climatic water deficit) | mm |
| `aet` | Sum of actual evapotranspiration | mm |
| `pet` | Sum of potential evapotranspiration | mm |
| `ppt` | Sum of modeled precipitation | mm |
| `runoff` | Sum of modeled excess water above soil storage capacity; **not routed streamflow** | mm |
| `melt`, `snowfall`, `water` | Sum of corresponding monthly water fluxes | mm |
| `snow` / `pack` | Mean September–May end-of-month snow-water equivalent | mm SWE |
| `snow_peak` | Maximum monthly snow-water equivalent in the September–May season | mm SWE |
| `soil` | Mean March–October end-of-month soil-water storage | mm |

`outputs_to_write` selects annual raster variables. The model also supports optional monthly TIFF variables such as `tmean`, `rs`, and `rnet`. Consult `names(wb_config())` and the installed defaults for supported controls and the `monthly_output` variable list.

**Period summaries.** When `cfg$summarize <- TRUE`, `results/summaries/` receives a mean raster for each selected output for the baseline and analysis periods. For flux variables, this is the **mean annual total**; for `snow` and `soil`, it is the **mean of their seasonal annual statistics**. The analysis-year multilayer stacks remain available separately. Missing incomplete snow seasons are excluded from means; valid zero-snow seasons are included.

## Seasonal snow and soil products

Snowpack and soil moisture are **stored states**, not precipitation or other monthly fluxes. A December snapshot is not a representative annual measure for many ecological applications.

- **Snow year:** September of the preceding calendar year through May of the labeled year. For example, snow year 2020 includes September–December 2019 and January–May 2020. The annual `snow`/`pack` output averages the nine end-of-month SWE states; `snow_peak` takes their maximum. These are not seasonal snowfall totals.
- **Soil-water year:** `soil` is the mean of March–October end-of-month soil-water storage in the labeled calendar year.
- **Incomplete first season:** When the preceding autumn is unavailable, the first snow-year product is missing (`NA`), not zero. This can occur in the first baseline year. Subsequent complete seasons are retained for climatological averages.

For finer seasonal analysis, request `pack` and `soil` in `zonal_monthly` and compute alternate date windows from the monthly history.

## Zonal time series

Zonal summaries are generated **during the simulation** from in-memory monthly results. They do not require saved monthly rasters and are suited to watersheds, ecological units and categorical maps.

```r
cfg$monthly_output <- FALSE
cfg$monthly_years <- NULL
cfg$zonal_path <- "/project/regions/ecological_types.tif"  # factor/category raster
cfg$zonal_id_field <- NULL
cfg$zonal_monthly <- c("cwd", "aet", "pet", "ppt", "soil", "pack", "runoff")
cfg$zonal_annual <- TRUE
```

The `zonal/` directory contains four primary tables:

| File | Description |
|---|---|
| `zonal_monthly.csv` | Monthly zonal means, spatial standard deviation, variance, coefficient of variation, and coverage. |
| `zonal_annual.csv` | Annual zonal summaries from the modeled monthly history. |
| `zonal_monthly_climatology.csv` | Mean seasonal cycle by zone, calendar month, and variable across the analysis years. Includes interannual variability (`temporal_sd`, `temporal_variance`) and mean within-zone dispersion (`mean_spatial_sd`, `mean_spatial_variance`), along with coverage and contributing-year counts. |
| `zonal_metadata.csv` | Zone identifiers, labels, and contributing-cell information. |

Monthly records identify zone, period (`baseline` or `analysis`), year, month, variable, spatial mean and valid-cell coverage. Monthly climatologies use **analysis years**, not the initialization baseline. Spatial variability measures differences among cells inside a zone; temporal variability measures differences among years in the zone's monthly mean. The coefficient of variation can be undefined where the mean is zero. Annual spatial variance must not be inferred by simply adding monthly variance estimates. For polygon zones, assignment is based on cell centers; disconnected cells with the same factor-raster category are treated as one zone. These summaries are spatial means, **not basin-integrated water volumes**. A zonal mean runoff depth is not a routed streamflow prediction.

## Prepared environment and restarts

With `prepared_environment = TRUE`, wbR stores prepared static geography for reuse between compatible simulations. The cache includes prepared DEM/terrain and soil-alignment products and is independent of the annual simulation output directory. A matching input geography and matching relevant preparation settings permit reuse; modifications to the input DEM, AOI, soil raster or grid parameters may invalidate the cache. **The monthly clear-sky radiation matrix is not yet persistently cached** and may be recalculated.

Restart checkpoints serve a different purpose: they record progress and carried soil/snow state during a particular simulation. A changed scientific-engine revision or incompatible configuration can invalidate those checkpoints even when the prepared-environment cache is still reusable.

## Computational requirements

The current core uses an **in-memory, full-grid simulation**. Its practical spatial limit depends on the available RAM, modeled cell count, raster alignment/reprojection requirements and concurrent temporary arrays. The observed peak resource demand is often during **initialization and environmental preparation**, when multiple large rasters are read or processed. The compiled monthly/yearly calculations are typically much faster after initialization. Writing and assembling many output rasters can impose another time, memory and disk-space cost.

Prepared-environment caching reduces repeated preprocessing but does **not** remove the RAM required to load prepared cell data for a run. wbR checks estimated memory demands and may refuse jobs too large for the configured budget. To reduce load, decrease the geographic extent, increase `dem_aggregate_fact`, reduce requested raster products or use a machine with more memory. The simple core does not promise automatic spatial tiling.

## Interpretation and limitations

wbR is a **vegetation-independent climatic water-balance model**. AET and CWD are modeled water-balance quantities, not direct observations of transpiration or plant water stress. The snow formulation uses monthly inputs and parameterized precipitation phase and melt. Terrain radiation incorporates slope and aspect but does not presently account for horizon shading by surrounding topography. Atmospheric transmittance and humidity may use monthly climatological normals rather than year-specific observations. The model's runoff variable represents water exceeding modeled storage capacity and is **not routed through a stream network**. Validate modeled results against measurements and independent data before interpreting local-scale predictions.

## Scientific references

- Stephenson, N. L. (1990). Climatic control of vegetation distribution: The role of the water balance. *The American Naturalist*, 135, 649–670. https://doi.org/10.1086/285067
- Stephenson, N. L. (1998). Actual evapotranspiration and deficit: Biologically meaningful correlates of vegetation distribution across spatial scales. *Journal of Biogeography*. https://doi.org/10.1046/j.1365-2699.1998.00233.x
- Lutz, J. A., van Wagtendonk, J. W., & Franklin, J. F. (2010). Climatic water deficit, tree species ranges, and climate change in Yosemite National Park. *Journal of Biogeography*, 37, 936–950. https://doi.org/10.1111/j.1365-2699.2009.02268.x
- McCune, B., & Keon, D. (2002). Equations for potential annual direct incident radiation and heat load. *Journal of Vegetation Science*, 13, 603–606. https://doi.org/10.1111/j.1654-1103.2002.tb02087.x
- Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large-scale parameters. *Monthly Weather Review*, 100, 81–92. https://doi.org/10.1175/1520-0493(1972)100%3C0081:OTAOSH%3E2.3.CO;2
- Hamon, W. R. (1961). Estimating potential evapotranspiration. *Journal of the Hydraulics Division*, 87, 107–120.
- Kustas, W. P., Rango, A., & Uijlenhoet, R. (1994). A simple energy budget algorithm for the snowmelt runoff model. *Water Resources Research*, 30, 1515–1527. https://doi.org/10.1029/94WR00152
- Hock, R. (2003). Temperature index melt modelling in mountain areas. *Journal of Hydrology*, 282, 104–115. https://doi.org/10.1016/S0022-1694(03)00257-9
- Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). *Crop evapotranspiration: Guidelines for computing crop water requirements*. FAO Irrigation and Drainage Paper 56. https://www.fao.org/4/X0490E/X0490E00.htm
- Daly, C., et al. (2008). Physiographically sensitive mapping of climatological temperature and precipitation across the conterminous United States. *International Journal of Climatology*, 28, 2031–2064. https://doi.org/10.1002/joc.1688

---

For questions or issue reports, use the [wbR repository issue tracker](https://github.com/drymixedconifer/wbR/issues). Cite the software using the repository's `CITATION.cff` file.
