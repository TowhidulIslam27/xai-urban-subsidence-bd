# Do the drivers of urban subsidence transfer between cities?

Code for the manuscript `manuscript_xai_v1`. We wrote one script per step, with shared helpers in `common/`. All paths are relative to `XAI_Paper` (see `common/paths.py`). We save data products to `analysis/` and figures to `manuscript/Figures` as PNG at 600 dpi.

**What we used:** Python 3 with numpy, scipy, pandas, matplotlib, Pillow and tifffile. We did not need GDAL. We run the whole pipeline in order with `python run_all.py`, which takes about 3.5 minutes on our machine. Each script can also be run alone, and each one takes city codes as arguments (for example, `python geology.py dha`).

## Scripts

| Section | Script | What we did | Output |
|---|---|---|---|
| 3.1 / 4.1 | `make_grid.py` | We fitted a 100 m lon/lat lattice to the InSAR cells and extended it to the city outlines in `analysis/City_Maps`. Our prediction area is every cell within 1 km of a PS cell, plus every cell inside the city. Northern Chattogram has no PS (it falls in another S1 frame), so we only predicted there and never trained on it. | `analysis/factors/<c>_grid.json`, `parts/<c>_lattice.npz` |
| 3.2 Terrain and hydrology | `terrain_hydrology.py` | We used FABDEM where `analysis/gee/exports/<c>_fabdem.tif` existed and the provisional SRTM ground proxy where it did not. From the DEM we derived slope, aspect, curvature, TRI, TPI, relief and ERR. We then ran priority-flood and D8 routing to get TWI, HAND (≥ 0), depression depth and streams. | `parts/<c>_terrain.npz`, `<c>_hydrology.npz` |
| 3.2 Geology | `geology.py` | We took the GSB 1:1M units and computed the Holocene share and the signed distance to the Holocene contact (±2 km). | `parts/<c>_geology.npz` |
| 3.2 Land cover, water | `landcover_water.py` | We used the 2024 city land cover as a check only, and we made it coverage-aware: shares are NaN where less than 50 % of the window is mapped, and distances are censored. We also computed distance to permanent water from the GSW mask (GEE). | `parts/<c>_landcover.npz` |
| 3.2 Urban load | `osm_factors.py` | We built OSM road density, distance to main road and railway, and facilities. We added city-specific layers for the embankment, MRT/expressway (DHA) and port (CTG). All distances are capped at 2 km. | `parts/<c>_osm.npz` |
| 3.2 (check) | `ctg_kit_check.py` | We put the Chattogram GLO-30 kit, soil and soil-moisture layers on the lattice. | `parts/ctg_ctgkit.npz` |
| 3.2 (check) | `wells_dhaka.py` | We added the Dhaka borehole/well points (DW001-DW160): distance, 1 km count and 500 m kernel density. We used these as a check only, since they come from a survey network and do not measure pumping. | `parts/dha_wells.npz`, columns in the Dhaka tables, `screening/wells_dha_check.csv` |
| 3.2 / 4.1 | `build_table.py` | We joined all the parts and added the target (vU …), the paper building factors (_ps), shared geology classes, 1 km CV blocks and vU classes. We also recorded which DEM and water source went in. | `analysis/factors/<c>_factors.csv`, `<c>_maps.npz`, `<c>_sources.json` |
| 3.2 GEE factors | `make_gee_script.py` → `gee_export.js` | We wrote the Earth Engine export script for the lattice and ran the `.js` in the Code Editor. | CSVs/tifs in Drive → `analysis/gee/exports` |
| 3.2 GEE factors | `merge_gee.py` | We merged the GEE exports by lattice row/col. | `analysis/factors/<c>_factors_full.csv` |
| 3.2, Table 3 | `factor_catalogue.py` | We catalogued every factor with its source, resolution, status per city, set (M0–M3), transfer and flags. | `analysis/screening/factor_catalogue.csv` |
| 3.2 / 5.1, Table 4 | `screening.py` | On the PS cells we screened out factors with resolution coarser than 1 km, no spread, coordinate-like behaviour (trend R² > 0.7), Spearman \|ρ\| > 0.8 or VIF > 10. What survived became our transfer set. | `analysis/screening/screen_<c>.csv`, `kept_factors.csv` |
| 5.1, Table 4 | `table4.py` | We generated the LaTeX body of Table 4. | `analysis/screening/table4_screening.tex` |
| 4.2 / 5.2, Tables 5-6, Fig. 4 | `cv_models.py` | For RQ1 we fitted nested sets M0-M3a on the PS cells with GBM (LightGBM if installed, otherwise scikit-learn HistGradientBoosting), RF and OLS, plus a coordinates-only reference. We used spatial block CV with 2 km blocks (5 folds × 10 repeats), and also checked block-size sensitivity (1/2/4 km) and the residual variogram range. | `analysis/models/` (`cv_`, `cv_summary`, `step_gain`, `block_sensitivity`, `range_check`, `oof_`, `table5_cv.tex`, `table6_gain.tex`), `Figures/Fig04_cv_skill.png` |
| 2, Fig. 1 | `fig1_study_area.py` | We mapped the location, surface geology (GSB 1:1M, nine grouped units), vU at the PS cells with the GNSS stations, and the standard error of vU. | `Figures/Fig01_study_area.png` |
| 3.2, Fig. S | `fig_dem_check.py` | We compared the SRTM ground proxy with GLO-30 derivatives for Chattogram. | `Figures/FigS_dem_check.png` |
| 5.1, Fig. 3 | `fig3_factor_maps.py` | We drew six key factors, with Dhaka beside Chattogram. | `Figures/Fig03_factor_maps.png` |
| 5.1, Fig. S | `fig_factor_atlas.py` | We plotted all candidate factors in the Fig. 3 layout (six per sheet, clipped to the city boundaries). The land-cover check factors and the PS-only building factors each got their own sheets. | `Figures/FigS_factor_atlas_*.png` |
| 5.1, Fig. S | `fig_factor_correlation.py` | We computed the Spearman matrix of the candidate factors for each city. | `Figures/FigS_factor_correlation.png` |
| 5.1, Fig. S | `fig_factor_association.py` | We computed the Spearman ρ of each factor with vU for each city. | `Figures/FigS_factor_association.png` |

**Shared helpers in `common/`:** `paths.py` (folders), `lattice.py` (lattice, pixel-to-cell aggregation, part files), `geo.py` (UTM, shapefile/GeoPackage readers, rasterisation, SRTM), `hydro.py` (flow routing and terrain derivatives), `factor_groups.py` (screening priority and groups), `figlib.py` and `figstyle.py` (figure style and maps), and `gee_template.js` (GEE export template).

## After the GEE export

Once the export finished, we copied every CSV and the three tifs per city (`<c>_glo30.tif`, `<c>_fabdem.tif`, `<c>_gsw_water.tif`) into `analysis/gee/exports` and re-ran `python run_all.py`. The terrain then switched to FABDEM and `dist_water` appeared on its own, and `<c>_sources.json` recorded what we had used. When we wanted to try another DEM, we set `DEM_DHA=path` (and `DEM_BARE=1` if it was bare earth).

## Map figures

**v2.** We draw each factor over the prediction area (the cells the models use) on a dark grey canvas (land `#3b3b3b`, water `#1d2024`, OSM roads, thana lines), with about 1.8 km of canvas around the larger city. Water comes from the GSW export when it is present. Otherwise we take it from the gaps between thana polygons (sea and rivers in Chattogram) plus the open water of the 2024 city maps, and we leave lattice cells that are mostly water to the basemap. Light grey means land without source data. Colour scales come from the prediction areas. The code is `basemap()`, `draw_basemap()` and `draw_map()` in `common/figlib.py`.

**v3 (19 Oct 2026).** We moved every map figure (Fig. 1, Fig. 3 and the atlas) to the Fig. 3 layout: rows of Dhaka | Chattogram pairs, two pairs per row, with one colour bar or legend under each pair. Maps are clipped to the city corporation outlines (`analysis/City_Maps`, white line, ward lines inside) using `map_sheet(..., clip=True)` and `draw_map(..., clip=True)`. Colour scales now come from the cells inside the two cities. Where a factor exists for only one city, we mark the other panel "not defined".

## Changes in v2 (18 Oct 2026)

The v1 code and outputs are kept in `analysis/_v1_snapshot`.

- We required model candidates to cover every cell. The 2024 city land-cover factors became check only: coverage separates core from fringe (median vU −0.4 vs +1.3 mm/yr in Dhaka), and a 2024 map shows surviving ponds, not fill. We only compute shares where at least 50 % of the window is mapped, and we censor distances at the map edge.
- We replaced distance to river (LC + Natural Earth) with distance to permanent water (GSW, GEE).
- We capped all distances at 2 km, and the Holocene-contact distance at ±2 km.
- We made embankment, MRT-6, expressway and port city-specific, so they do not transfer. In Chattogram the OSM embankments are mostly the Outer Ring Road, so we left out railway-named features.
- In screening, we dropped layers coarser than 1 km and factors that act as coordinates (trend R² > 0.7; distance to faults 0.97, to the coast 1.00). The PS-only building factors are check only.
- In the GEE script we used the urban epoch without the 2030 fill value, took water loss 1990–2022 from one Landsat method, counted GSW water lost inside the InSAR window, and used height change of existing buildings only. We also switched to bilinear DEMs, ran Open Buildings as its own task (A_ob), added Landsat scene counts and fuller cloud masks, and added a GSW water mask.
