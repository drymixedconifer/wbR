# wbR engine v6.4: merge notes

**Status: candidate. Built by static merge; not run in R.**

## What v6.4 is

The Oct 26 package engine (`50ef335`, "snow evaporation update") plus the v6.3 additions that the Nov 19 swap (`5df9625`) brought in. Everything the swap removed is kept: zonal series, `balance`, monthly `runoff`, the seasonal snow outputs (`snow`, `snow_peak`) and their first-year initialization, all-NA band handling, the initialization cache, the prepared-environment cache, PRISM subset mode, and the `summarize` / `write_baseline` outputs. The snow-cover logic in `50ef335` and v6.3 is identical (only comments differ), so nothing was lost there.

## Taken from v6.3 (57 of 94 differing hunks)

- Comment and doc rewording (about 40 hunks), code-neutral.
- AWC: constant storage (`AWC_CONSTANT_MM`), depth raster (`AWC_DEPTH_PATH`, units, resampling), `mm_per_m` units, soil provenance fields.
- PRISM cache residency check (`wb_check_prism_resident`, `PRISM_CHECK_RESIDENT`).
- `temp_root` setting.
- Provenance fields in `site_metadata.rds` (melt method, snow settings).
- VPD constant (`wb_vpd_hpa_to_kpa <- 0.1`): `(hi + lo) / 2 * 0.1` equals the old `(hi + lo) / 20`.
- `analysis_initialization` as an option (`"baseline_end"` default, `"equilibrium"`).

## Deliberate decisions (change if you disagree)

1. **`analysis_initialization = "equilibrium"` is refused together with the `snow` / `snow_peak` outputs.** Those need the preceding autumn's snow history, which only `"baseline_end"` carries into the first analysis year. Equilibrium also forces a fresh spin-up (it skips the cached baseline-end state).
2. **Cache directory kept as `water_balance_v6`.** v6.3 renamed it `water_balance_v6_2`, which would orphan restart caches.
3. **Both caches now key on the new soil options.** The prepared-environment key and the initialization-state key ignored the depth raster, constant storage and resampling, so changing them could have reused a stale soil grid. Environment key version bumped to 2, so existing prepared environments rebuild once.
4. **Not taken:** v6.3's extracted `wb_write_monthly()` (it removed the zonal collection from the monthly loop) and the v6.3 limit of annual outputs to `cwd`/`aet`. Hunks 12, 30, 58 (snow-cover validation) were skipped as duplicates of code already present.

## Check before relying on it

Run in R, in this order. Use a small AOI and a few years.

1. `wb_config(advanced = TRUE)`, then `stop_after_setup = TRUE`. Confirms it parses and sources, including the C++ block.
2. **Regression:** run `50ef335` and 6.4 on the same inputs with the AWC raster path (not the new options). Expect identical `cwd` / `aet` annual rasters, identical zonal CSVs, identical `snow` / `snow_peak`. Any difference is a merge error.
3. New AWC options: `AWC_CONSTANT_MM`, then a depth raster. After changing a soil option, confirm the environment cache is rebuilt.
4. `analysis_initialization = "equilibrium"` with `cwd` / `aet` only, then confirm the snow guard stops it with `snow`.
5. Future extension: it was written for v6.2, so a short LOCA2 test is still needed.

## Land it on a branch

```
git checkout -b wbr-6.4
# unzip the candidate over the repo (or: git apply wbR-6.4.patch)
git add -A && git commit -m "wbR engine 6.4: package engine + v6.3 additions"
git push -u origin wbr-6.4
```
