"""Section 3.2  Google Earth Engine factors - write the export script for the GEE Code Editor.

Fills common/gee_template.js with the lattice of each city (analysis/factors/<c>_grid.json) and
writes s32g_gee_export.js (all exports) and s32g_gee_export_part3.js (only the exports still
missing after two runs: B_gsw, B_ls, B_dw, C_s2, C_lst). Part 2 (A_ob1/2, C_dw) is done. Run that script in code.earthengine.google.com and start the tasks:
  A_urban   Open Buildings 2.5D Temporal, GHSL built surface/volume/height, WorldPop, VIIRS
  B_fill    JRC GSW 1.4 fill year, Landsat 5 water 1990/2005, Dynamic World water/built change
  C_surface Sentinel-2 dry-season indices, Dynamic World 2017-2023, Landsat 8/9 LST
  D_soil    OpenLandMap clay/sand/bulk density/organic carbon/texture, MERIT HAND, CHIRPS, SMAP L4
  E_dem     Copernicus GLO-30 and FABDEM at 30 m UTM 46N (GeoTIFF)
Every layer is aggregated onto the lattice and exported with lattice row/col, so the merge is exact.
Put the files from the Drive folder XAI_factors into analysis/gee/exports.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json
from paths import FACTORS, CODE, CITIES

g = {c: {k: v for k, v in json.load(open(FACTORS / f'{c}_grid.json')).items() if k != 'note'} for c in CITIES}
tpl = open(CODE / 'common' / 'gee_template.js').read().replace('__GRIDS__', json.dumps(g, indent=1))
ALL = ['A_ob', 'A_urban', 'B_gsw', 'B_ls', 'B_dw', 'C_s2', 'C_lst', 'C_dw', 'D_soil', 'E_dem']
# full script, and part 3 = the exports still missing after two runs (19 Oct 2026):
# run 1 gave A_urban, D_soil and the GeoTIFFs; run 2 gave A_ob1/A_ob2 and C_dw;
# B_fill and C_s2 failed and are now split into B_gsw, B_ls, B_dw, C_s2, C_lst
RUNS = {'s32g_gee_export.js': ALL, 's32g_gee_export_part3.js': ['B_gsw', 'B_ls', 'B_dw', 'C_s2', 'C_lst']}
for fn, on in RUNS.items():
    flags = json.dumps({k: k in on for k in ALL})
    open(CODE / fn, 'w').write(tpl.replace('__EXPORT__', flags))
    print('written', CODE / fn, '| exports:', ', '.join(on))
