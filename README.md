# Do the drivers of urban subsidence transfer between cities?

**Explainable machine learning on Sentinel-1 PS-InSAR in Dhaka and Chattogram, Bangladesh**

Towhidul Islam · SAR.Sense Geointelligence Lab, Chittagong · work in progress (manuscript in preparation)

> The figures on this page are previews of a paper that has not been published yet.
> They are reduced in resolution and watermarked. Please don't reuse or cite them.

Dhaka and Chattogram are both sinking, but for different reasons. This project asks what
controls where the ground goes down in each city, and whether a model of those controls
trained in one city still works in the other. The target is the GNSS-validated vertical
velocity (vU, 2017–2023) from my Sentinel-1 PS-InSAR processing, on 100 m cells; the same
deformation data is on the [Bangladesh Ground Motion Observatory](https://github.com/TowhidulIslam27/bGMO).

![Study area](figures/Fig01_study_area.jpg)
*Fig. 1. (a) Location; (b) surface geology; (c) vertical velocity vU at the PS cells with the GNSS
stations used for validation; (d) its standard error.*

## Research questions

1. **Skill.** How much of the subsidence pattern do nested groups of factors explain, under
   spatial cross-validation: terrain, hydrology and soil (M0), + geology (M1), + urban load and
   land cover (M2), + fill history and groundwater (M3)?
2. **Drivers.** Which factors matter, and how (SHAP: global, grouped, dependence and local
   explanations for the known hotspots)?
3. **Transfer.** Does a model trained in Dhaka predict Chattogram, and the other way round, using
   only factors that mean the same thing in both cities?

## What is done

- **A common 100 m lattice** for both cities, covering the whole city corporation area
  (Dhaka 68 443 cells, 17 622 with vU; Chattogram 32 752 cells, 5 142 with vU).
- **64 candidate factors per city** from FABDEM, JRC Global Surface Water, Landsat, Dynamic World,
  Sentinel-2, Google Open Buildings, GHSL, OpenStreetMap, OpenLandMap, HYSOGs and the Geological
  Survey of Bangladesh map, built with Python and Google Earth Engine.
- **Transparent screening** before any model sees vU: resolution, coverage, near-constant values,
  a test for factors that only act as map coordinates, Spearman |ρ| > 0.8 and VIF > 10.
  Dhaka keeps 52 factors, Chattogram 45, and 40 are shared for the transfer test.
- **RQ1 models**: gradient-boosted trees, random forest and linear regression on the nested
  sets, with 2 km spatial-block cross-validation (5 folds × 10 repeats). The block size comes
  from the residual semivariogram, and 1 km and 4 km blocks are run as a sensitivity check.

![Factor maps](figures/Fig03_factor_maps.jpg)
*Fig. 3. Six of the factors, Dhaka beside Chattogram.*

## First results (RQ1)

![Cross-validated skill](figures/Fig04_cv_skill.jpg)
*Fig. 4. Cross-validated R² of vU for the nested factor sets (mean ± s.d. over 10 repeats).
Dashed line: a model given only the cell coordinates.*

| Set (gradient-boosted trees) | Dhaka R² | Chattogram R² |
|---|---|---|
| M0 terrain, hydrology, soil | 0.11 | 0.38 |
| M1 + geology | 0.22 | 0.44 |
| M2 + urban load, land cover | 0.46 | 0.64 |
| M3a + fill history, groundwater proxies | 0.50 | 0.64 |
| Coordinates only (reference) | 0.61 | 0.74 |

- The two cities are controlled differently. In Chattogram, terrain alone explains a large part
  of the pattern (low ground subsides, the hills are stable). In Dhaka it explains almost
  nothing, and urban load, land cover and fill history carry the signal.
- Every step adds skill in every repeat, and the order holds for 1, 2 and 4 km blocks. Tree
  models beat linear regression by 0.3 in R², so the effects are non-linear.
- A model that only knows where a cell is still does better than all the factors together.
  Part of the subsidence follows a broad spatial pattern that the current factors miss. In
  Dhaka the obvious candidate is the cone of depression in the Dupi Tila aquifer, so the next
  model (M3b) adds groundwater heads digitised from published piezometric maps. I report this
  baseline on purpose: it is the honest test of what the factors explain.

## Next

Groundwater heads for Dhaka (M3b), SHAP explanations, the Dhaka ↔ Chattogram transfer test,
susceptibility maps with an area of applicability, and a 50 m / 200 m scale check.

## Code

`code/` holds the full pipeline, one script per manuscript section (see
[code/README.md](code/README.md)). It runs with NumPy, SciPy, pandas, scikit-learn and
Matplotlib, without GDAL; the Earth Engine scripts export the satellite factors onto the same
lattice. The input data and the derived cell tables are not included while the papers are in
preparation.

| Step | Scripts |
|---|---|
| Lattice and study area | `s31_*`, `s21_fig01_study_area.py` |
| Factors (terrain, geology, land cover, OSM, soil, GEE) | `s32a`–`s32h`, `s32m1` |
| Catalogue and screening | `s32i`, `s32j`, `s32l` |
| Figures | `s51_*` |
| RQ1 models and cross-validation | `s52_models_rq1.py` |

## Status and use

This is ongoing research shared to show the work and the methods. Please contact me before
reusing the code or figures; a licence will be added when the paper is published.

© 2026 Towhidul Islam · [portfolio](https://towhidulislam27.github.io/) · m.towhid92@gmail.com
