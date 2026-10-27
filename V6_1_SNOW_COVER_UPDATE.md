# wbR v6.1 snow-cover update (development)

- Default `snow_cover_demand = "exclude"`, `snow_cover_min_mm = 10`. Snow-covered fraction of each month reduces soil/plant atmospheric demand; reported PET remains total atmospheric demand. Thus CWD + AET can be below PET during snow cover.
- Set `snow_cover_demand = "soil"` to recover the previous v6 soil-demand treatment.
- Monthly `scover` (fraction 0–1) is selectable in monthly and zonal outputs.
- Preserved `runoff` at bit 8192, placed `scover` at bit 16384; expanded C++ mask to 32767.
- Model physics changed: historical results will not match v6 unless compatibility mode is selected. Initialization scientific fingerprint includes new settings.
- No Git actions performed. R runtime and numerical tests remain necessary.
