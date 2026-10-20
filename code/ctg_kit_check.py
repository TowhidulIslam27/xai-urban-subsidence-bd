"""Section 3.2  Chattogram GLO-30 kit layers, sampled onto the lattice as cross-checks.

Source    Dataset/outputs/Final_Derivatives (30 m, EPSG:32646) and Dataset/SM_*.tif (soil moisture).
Use       glo_* terrain for the DEM check (Fig. S dem_check); hsg, clay, bsr, sm_* as CTG-only
          checks; dist_fault and dist_coast are recorded but left out of the models
          (regional gradients that act as coordinates; no coast in Dhaka).
Output    analysis/factors/parts/ctg_ctgkit.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import numpy as np
from geo import Grid, sample_regular
from lattice import Ctx, Agg, save_part, read_tif
from paths import DATASET

LAYERS = [('F06_Elevation', 'glo_elev'), ('F07_Slope_degrees_new', 'glo_slope'), ('F09_Profile_Curvature', 'glo_curv_prof'),
          ('F10_Plan_Curvature', 'glo_curv_plan'), ('F11_TWI', 'glo_twi'), ('F12_TRI', 'glo_tri'),
          ('F13_ElevationReliefRatio', 'glo_err'), ('F14_Dist_NearestFault', 'dist_fault'),
          ('F16_Clay_weighted_mean', 'clay'), ('F21_Distance_to_Rivers', 'glo_dist_river'),
          ('F23_DrainageDensity', 'glo_drain_dens'), ('F24_Distance_to_Coastline', 'dist_coast'), ('F30_BSR', 'bsr')]
SM = [('SM_annual_mean_2017_2024', 'sm_annual'), ('SM_annual_std_2017_2024', 'sm_std'),
      ('SM_wetseason_mean_MayOct_2017_2024', 'sm_wet'), ('SM_dryseason_mean_NovApr_2017_2024', 'sm_dry'),
      ('SM_rootzone_mean_2017_2024', 'sm_root')]


def run(c='ctg'):
    x = Ctx(c)
    fd = DATASET / 'outputs' / 'Final_Derivatives'
    a, tx, ty, sx, sy = read_tif(fd / 'F06_Elevation.tif')
    AF = Agg(x.L, Grid(tx, ty, sx, a.shape[1], a.shape[0]))
    F = {nm: AF.mean(read_tif(fd / f'{fn}.tif')[0]) for fn, nm in LAYERS}
    F['hsg'] = AF.majority(np.nan_to_num(read_tif(fd / 'F15_HSG_HydroSoilGroup.tif')[0], nan=-1), 16)
    for fn, nm in SM:
        s, l0, la0, dx, dy = read_tif(DATASET / f'{fn}.tif')
        F[nm] = sample_regular(np.nan_to_num(s, nan=np.nanmean(s)), l0, la0, dx, dy, x.LON, x.LAT)
    save_part(c, 'ctgkit', F); x.log('CTG kit', list(F))


if __name__ == '__main__':
    run('ctg')
