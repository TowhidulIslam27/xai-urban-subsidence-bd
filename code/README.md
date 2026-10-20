# Code for "Do the drivers of urban subsidence transfer between cities?" (manuscript_xai_v1)

One script per manuscript section. Shared helpers are in `common/`. Paths are relative to
`XAI_Paper` (see `common/paths.py`); data products go to `analysis/`, figures to
`manuscript/Figures` as PNG at 600 dpi.

Requirements: Python 3 with numpy, scipy, pandas, matplotlib, Pillow, tifffile (no GDAL needed).
Run everything in order with `python run_all.py` (about 3.5 minutes), or one script at a time;
each script takes city codes as arguments (`python s32b_geology.py dha`).

| Section | Script | What it does | Output |
|---|---|---|---|
| 3.1 / 4.1 | `s31_target_lattice.py` | fits the 100 m lon/lat lattice of the InSAR cells, extended to the city outlines in `analysis/City_Maps`; prediction area = cells within 1 km of a PS cell plus every cell inside the city (northern Chattogram: no PS, other S1 frame; predicted, never trained on) | `analysis/factors/<c>_grid.json`, `parts/<c>_lattice.npz` |
| 3.2 Terrain and hydrology | `s32a_terrain_hydrology.py` | FABDEM if `analysis/gee/exports/<c>_fabdem.tif` exists, else the provisional SRTM ground proxy; slope, aspect, curvature, TRI, TPI, relief, ERR; priority-flood + D8 → TWI, HAND (≥ 0), depression depth, streams | `parts/<c>_terrain.npz`, `<c>_hydrology.npz` |
| 3.2 Geology | `s32b_geology.py` | GSB 1:1M units, Holocene share, signed distance to the Holocene contact (±2 km) | `parts/<c>_geology.npz` |
| 3.2 Land cover, water | `s32c_landcover_fill.py` | 2024 city land cover, coverage-aware and check only: shares (NaN if < 50 % of the window mapped), censored distances; distance to permanent water from the GSW mask (GEE) | `parts/<c>_landcover.npz` |
| 3.2 Urban load | `s32d_urban_osm.py` | OSM road density, distance to main road, railway; facilities; city-specific: embankment, MRT/expressway (DHA), port (CTG); distances capped at 2 km | `parts/<c>_osm.npz` |
| 3.2 (check) | `s32e_ctg_kit_crosscheck.py` | Chattogram GLO-30 kit, soil and soil-moisture layers on the lattice | `parts/ctg_ctgkit.npz` |
| 3.2 (check) | `s32m1_wells_dha.py` | Dhaka borehole/well points (DW001-DW160): distance, 1 km count, 500 m kernel density; check only (survey network, not pumping) | `parts/dha_wells.npz`, columns in the Dhaka tables, `screening/wells_dha_check.csv` |
| 3.2 / 4.1 | `s32f_assemble_table.py` | joins all parts, target (vU …), paper building factors (_ps), shared geology classes, 1 km CV blocks, vU classes; records the DEM and water source | `analysis/factors/<c>_factors.csv`, `<c>_maps.npz`, `<c>_sources.json` |
| 3.2 GEE factors | `s32g_make_gee_script.py` → `s32g_gee_export.js` | writes the Earth Engine export script for the lattice; run the .js in the Code Editor | CSVs/tifs in Drive → `analysis/gee/exports` |
| 3.2 GEE factors | `s32h_merge_gee.py` | merges the GEE exports by lattice row/col | `analysis/factors/<c>_factors_full.csv` |
| 3.2, Table 3 | `s32i_factor_catalogue.py` | catalogue of every factor: source, resolution, status per city, set M0–M3, transfer, flags | `analysis/screening/factor_catalogue.csv` |
| 3.2 / 5.1, Table 4 | `s32j_screening.py` | on PS cells: resolution > 1 km, no spread, acts as a coordinate (trend R² > 0.7), Spearman \|ρ\| > 0.8, VIF > 10; transfer set | `analysis/screening/screen_<c>.csv`, `kept_factors.csv` |
| 5.1, Table 4 | `s32l_table4_tex.py` | LaTeX body of Table 4 | `analysis/screening/table4_screening.tex` |
| 4.2 / 5.2, Tables 5-6, Fig. 4 | `s52_models_rq1.py` | RQ1: nested sets M0-M3a on the PS cells, GBM (LightGBM if installed, else scikit-learn HistGradientBoosting), RF and OLS, plus a coordinates-only reference; spatial block CV with 2 km blocks, 5 folds x 10 repeats; block-size sensitivity (1/2/4 km) and residual variogram range | `analysis/models/` (cv_*, cv_summary, step_gain, block_sensitivity, range_check, oof_*, table5_cv.tex, table6_gain.tex), `Figures/Fig04_cv_skill.png` |
| 2, Fig. 1 | `s21_fig01_study_area.py` | location, surface geology (GSB 1:1M, nine grouped units), vU at the PS cells with the GNSS stations, standard error of vU | `Figures/Fig01_study_area.png` |
| 3.2, Fig. S | `s32k_figS_dem_check.py` | SRTM ground proxy vs GLO-30 derivatives (Chattogram) | `Figures/FigS_dem_check.png` |
| 5.1, Fig. 3 | `s51_fig03_factor_maps.py` | six key factors, Dhaka beside Chattogram | `Figures/Fig03_factor_maps.png` |
| 5.1, Fig. S | `s51_figS_factor_atlas.py` | all candidate factors in the Fig. 3 layout (six per sheet, clipped to the city boundaries), then the land-cover check factors and the PS-only building factors on their own sheets | `Figures/FigS_factor_atlas_*.png` |
| 5.1, Fig. S | `s51_figS_factor_correlation.py` | Spearman matrix of the candidate factors per city | `Figures/FigS_factor_correlation.png` |
| 5.1, Fig. S | `s51_figS_factor_association.py` | Spearman ρ of each factor with vU per city | `Figures/FigS_factor_association.png` |

`common/`: `paths.py` (folders), `lattice.py` (lattice, pixel-to-cell aggregation, part files),
`geo.py` (UTM, shapefile/GeoPackage readers, rasterisation, SRTM), `hydro.py` (flow routing and terrain
derivatives), `factor_groups.py` (screening priority and groups), `figlib.py` and `figstyle.py`
(figure style and maps), `gee_template.js` (GEE export template).

After the GEE export: put every CSV and the three tifs per city (`<c>_glo30.tif`, `<c>_fabdem.tif`,
`<c>_gsw_water.tif`) in `analysis/gee/exports` and run `python run_all.py`. Terrain switches to FABDEM
and `dist_water` appears by themselves; `<c>_sources.json` records what was used. To force another DEM:
`set DEM_DHA=path` (add `set DEM_BARE=1` if it is bare earth).

Map figures (v2): each factor is drawn over the prediction area (the cells the models use) on a dark
grey canvas (land #3b3b3b, water #1d2024, OSM roads, thana lines), ~1.8 km of canvas around the larger city.
Water comes from the GSW export when present, otherwise from the gaps between thana polygons (sea and
rivers in Chattogram) plus the open water of the 2024 city maps; lattice cells that are mostly water
are left to the basemap. Light grey = land without source data. Colour scales come from the
prediction areas. Code: `basemap()`, `draw_basemap()` and `draw_map()` in `common/figlib.py`.

Map figures (v3, 19 Oct 2026): every map figure (Fig. 1, Fig. 3, the atlas) uses the Fig. 3 layout:
rows of Dhaka | Chattogram pairs, two pairs per row, one colour bar or legend under each pair, maps
clipped to the city corporation outlines (`analysis/City_Maps`, white line, ward lines inside) with
`map_sheet(..., clip=True)` / `draw_map(..., clip=True)`. Colour scales come from the cells inside the
two cities. A factor defined for one city only is marked "not defined" in the other panel.

Changes in v2 (18 Oct 2026; v1 code and outputs are in `analysis/_v1_snapshot`):
- Model candidates must cover every cell. The 2024 city land-cover factors are check only: coverage
  separates core from fringe (median vU −0.4 vs +1.3 mm/yr in Dhaka), and a 2024 map shows surviving
  ponds, not fill. Shares need ≥ 50 % of the window mapped; distances are censored at the map edge.
- Distance to river (LC + Natural Earth) is replaced by distance to permanent water (GSW, GEE).
- All distances are capped at 2 km; the Holocene-contact distance at ±2 km.
- Embankment, MRT-6, expressway and port are city-specific (no transfer). In Chattogram the OSM
  embankments are mostly the Outer Ring Road; railway-named features are left out.
- Screening drops layers coarser than 1 km and factors that act as coordinates (trend R² > 0.7;
  distance to faults 0.97, to the coast 1.00); PS-only building factors are check only.
- GEE script: urban epoch without the 2030 fill value; water loss 1990–2022 from one Landsat method;
  GSW water lost inside the InSAR window; height change of existing buildings only; bilinear DEMs;
  Open Buildings as its own task (A_ob); Landsat scene counts; fuller cloud masks; GSW water mask.
