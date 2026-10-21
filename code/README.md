# Code for "Do the drivers of urban subsidence transfer between cities?" (manuscript_xai_v1)

This folder has the code I used for the paper. I split the workflow into one script per step and
put the shared helpers in `common/`. All paths are relative to the `XAI_Paper` folder (set in
`common/paths.py`). Data products are written to `analysis/`, and figures are saved to
`manuscript/Figures` as 600 dpi PNG files.

I ran everything with Python 3 and numpy, scipy, pandas, matplotlib, Pillow and tifffile. GDAL is
not needed. `python run_all.py` runs the whole workflow in order, which took about 3.5 minutes on my
machine. The scripts can also be run one at a time. Each one takes city codes as arguments, for
example `python geology.py dha`.

| Section | Script | What I did | Output |
|---|---|---|---|
| 3.1 / 4.1 | `make_grid.py` | I fitted the 100 m lon/lat lattice of the InSAR cells and extended it to the city outlines in `analysis/City_Maps`. The prediction area is every cell within 1 km of a PS cell, plus every cell inside the city. Northern Chattogram has no PS cells because it falls in another Sentinel-1 frame, so the models predict there but are never trained there. | `analysis/factors/<c>_grid.json`, `parts/<c>_lattice.npz` |
| 3.2 Terrain and hydrology | `terrain_hydrology.py` | I used FABDEM when `analysis/gee/exports/<c>_fabdem.tif` is present and fell back to a provisional SRTM ground proxy when it is not. From the DEM I computed slope, aspect, curvature, TRI, TPI, relief and ERR. A priority-flood fill with D8 routing gave TWI, HAND (≥ 0), depression depth and streams. | `parts/<c>_terrain.npz`, `<c>_hydrology.npz` |
| 3.2 Geology | `geology.py` | I took the units from the GSB 1:1M map and computed the Holocene share and the signed distance to the Holocene contact (±2 km). | `parts/<c>_geology.npz` |
| 3.2 Land cover, water | `landcover_water.py` | I kept the 2024 city land cover as a check only. Shares are set to NaN where less than half of the window was mapped, and distances are censored at the map edge. Distance to permanent water comes from the GSW mask I exported from GEE. | `parts/<c>_landcover.npz` |
| 3.2 Urban load | `osm_factors.py` | From OSM I derived road density, distance to the nearest main road and railway, and facilities. A few factors exist in one city only: embankments, MRT and expressway in Dhaka, and the port in Chattogram. All distances are capped at 2 km. | `parts/<c>_osm.npz` |
| 3.2 (check) | `ctg_kit_check.py` | I put the Chattogram GLO-30 kit, soil and soil-moisture layers on the lattice. | `parts/ctg_ctgkit.npz` |
| 3.2 (check) | `wells_dhaka.py` | I used the Dhaka borehole and well points (DW001–DW160) to compute distance to the nearest well, the well count within 1 km and a 500 m kernel density. These are a check only, because the points come from a survey network and say nothing about pumping. | `parts/dha_wells.npz`, extra columns in the Dhaka tables, `screening/wells_dha_check.csv` |
| 3.2 / 4.1 | `build_table.py` | I joined all the parts into one table per city. It holds the target (vU and related fields), the building factors from the paper (`_ps`), the geology classes shared by both cities, the 1 km CV blocks and the vU classes. The script also records which DEM and which water source were used. | `analysis/factors/<c>_factors.csv`, `<c>_maps.npz`, `<c>_sources.json` |
| 3.2 GEE factors | `make_gee_script.py` → `gee_export.js` | I generated the Earth Engine export script from the lattice and ran the `.js` file in the GEE Code Editor. | CSVs and tifs in Google Drive, copied to `analysis/gee/exports` |
| 3.2 GEE factors | `merge_gee.py` | I merged the GEE exports into the factor table by lattice row and column. | `analysis/factors/<c>_factors_full.csv` |
| 3.2, Table 3 | `factor_catalogue.py` | I listed every factor with its source, resolution, status in each city, model set (M0–M3), whether it transfers, and any flags. | `analysis/screening/factor_catalogue.csv` |
| 3.2 / 5.1, Table 4 | `screening.py` | I screened the factors on the PS cells. I dropped layers coarser than 1 km, layers with no spread, layers that behave like a coordinate (trend R² > 0.7), and pairs with Spearman \|ρ\| > 0.8 or VIF > 10. What was left became the transfer set. | `analysis/screening/screen_<c>.csv`, `kept_factors.csv` |
| 5.1, Table 4 | `table4.py` | I wrote out the LaTeX body of Table 4. | `analysis/screening/table4_screening.tex` |
| 4.2 / 5.2, Tables 5–6, Fig. 4 | `cv_models.py` | This is the RQ1 step. I fitted the nested sets M0–M3a on the PS cells with GBM, RF and OLS, and added a coordinates-only model as a reference. GBM uses LightGBM if it is installed and scikit-learn's HistGradientBoosting otherwise. I used spatial block CV with 2 km blocks and 5 folds repeated 10 times, then tested 1, 2 and 4 km blocks and checked the range of the residual variogram. | `analysis/models/` (cv_*, cv_summary, step_gain, block_sensitivity, range_check, oof_*, table5_cv.tex, table6_gain.tex), `Figures/Fig04_cv_skill.png` |
| 2, Fig. 1 | `fig1_study_area.py` | I mapped the study area: location, surface geology (GSB 1:1M, grouped into nine units), vU at the PS cells with the GNSS stations, and the standard error of vU. | `Figures/Fig01_study_area.png` |
| 3.2, Fig. S | `fig_dem_check.py` | I compared derivatives from the SRTM ground proxy and GLO-30 in Chattogram. | `Figures/FigS_dem_check.png` |
| 5.1, Fig. 3 | `fig3_factor_maps.py` | I mapped six key factors, with Dhaka next to Chattogram. | `Figures/Fig03_factor_maps.png` |
| 5.1, Fig. S | `fig_factor_atlas.py` | I mapped all candidate factors in the Fig. 3 layout, six per sheet and clipped to the city boundaries. The land-cover check factors and the PS-only building factors have their own sheets. | `Figures/FigS_factor_atlas_*.png` |
| 5.1, Fig. S | `fig_factor_correlation.py` | I plotted the Spearman matrix of the candidate factors for each city. | `Figures/FigS_factor_correlation.png` |
| 5.1, Fig. S | `fig_factor_association.py` | I plotted the Spearman ρ between each factor and vU for each city. | `Figures/FigS_factor_association.png` |

The shared code in `common/` is split as follows: `paths.py` holds the folder locations, `lattice.py`
the lattice, pixel-to-cell aggregation and part files, `geo.py` the UTM conversion, shapefile and
GeoPackage readers, rasterisation and SRTM, and `hydro.py` the flow routing and terrain derivatives.
`factor_groups.py` sets the factor groups and screening priority, `figlib.py` and `figstyle.py` handle
figure style and map drawing, and `gee_template.js` is the GEE export template.

Once the GEE export had finished, I put every exported CSV in `analysis/gee/exports` together with
three tifs per city (`<c>_glo30.tif`, `<c>_fabdem.tif`, `<c>_gsw_water.tif`) and ran
`python run_all.py` again. With those files in place, the terrain step picks up FABDEM and
`dist_water` is added without any other change. `<c>_sources.json` records which inputs were used.
A different DEM can be forced with `set DEM_DHA=path`, adding `set DEM_BARE=1` if it is bare earth.

## Map figures

In v2 I drew each factor over the prediction area, meaning the cells the models use. The background
is a dark grey canvas (land #3b3b3b, water #1d2024) with OSM roads and thana boundaries on top, and it
extends about 1.8 km beyond the larger city. Water comes from the GSW export when I have it.
Otherwise I took it from the gaps between thana polygons (the sea and rivers around Chattogram) and
from the open water in the 2024 city maps. Lattice cells that are mostly water are left blank so the
basemap shows through, and light grey marks land with no source data. Colour scales were set from
the prediction areas. The code is in `basemap()`, `draw_basemap()` and `draw_map()` in
`common/figlib.py`.

In v3 (19 Oct 2026) I moved all the map figures (Fig. 1, Fig. 3 and the atlas) to the Fig. 3 layout.
Each row holds two Dhaka | Chattogram pairs, with one colour bar or legend under each pair. I clipped
the maps to the city corporation outlines in `analysis/City_Maps`, drawn as a white line with the
ward boundaries inside, using `map_sheet(..., clip=True)` and `draw_map(..., clip=True)`. Colour
scales now come from the cells inside the two cities. Where a factor exists for only one city, the
other panel says "not defined".

## What changed in v2 (18 Oct 2026)

The v1 code and outputs are kept in `analysis/_v1_snapshot`.

- I decided that a model candidate has to cover every cell, so the 2024 city land-cover factors
  became check only. Their coverage splits the city core from the fringe: in Dhaka the median vU is
  −0.4 mm/yr where the map has data and +1.3 mm/yr where it does not. A 2024 map also shows the
  ponds that are still there, not the ones that were filled. I kept shares only where at least half
  of the window was mapped and censored distances at the map edge.
- I replaced distance to river (from the land cover and Natural Earth) with distance to permanent
  water from GSW, exported through GEE.
- I capped every distance at 2 km, and the Holocene-contact distance at ±2 km.
- I treated embankment, MRT-6, expressway and port as city-specific, so they are not used for
  transfer. Most of the OSM embankments in Chattogram turned out to be the Outer Ring Road, and I
  left out features named as railways.
- In the screening I dropped layers coarser than 1 km and factors that behave like coordinates
  (trend R² > 0.7). Distance to faults (0.97) and distance to the coast (1.00) went for this reason.
  The PS-only building factors became check only.
- I changed the GEE script in several places. The urban epoch no longer uses the 2030 fill value,
  water loss from 1990 to 2022 now comes from a single Landsat method, and I added the GSW water lost
  inside the InSAR window. Building height change now covers existing buildings only, the DEMs are
  resampled bilinearly, and Open Buildings runs as its own task (`A_ob`). I also added Landsat scene
  counts, more complete cloud masks and an export of the GSW water mask.
