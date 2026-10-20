"""Section 3.2  Land cover and open water: city land-cover maps (check only) and distance to water.

Source    2024 land-cover maps of the city corporations (analysis/inputs/3_City_LC_2024), 10 m.
          They cover only part of each prediction area (lc_cover = QA), and the uncovered fringe
          has a very different vU from the covered core. Every factor from these maps is therefore
          coverage-aware and kept for checks on the covered area only; none is a model candidate
          (see common/factor_groups.py, LC_CHECK). A 2024 map also shows the ponds that survived,
          not the ones that were filled, so pond and development shares are land cover, not fill
          history; fill history comes from the GEE layers (GSW, Landsat, Dynamic World, GHSL).
Classes   open water (river, canal, lake, jheel, creek, mudflat), ponds/ditches/aquaculture,
          area under development, built, green/open, road, other.
Factors   lc_cover (QA); shares within 300 m (500 m for open water) of the mapped area, NaN where
          less than half the disk is mapped; distance to pond, to land under development and to
          mapped open water (dist_lc_water), censored where unmapped ground is closer than the
          nearest feature, capped at 2 km.
          dist_water: distance to permanent surface water (JRC GSW occurrence >= 50 %), full
          coverage, from the GEE export analysis/gee/exports/<c>_gsw_water.tif (skipped until the
          export is there). This replaces the old LC + Natural Earth distance to river.
Output    analysis/factors/parts/<c>_landcover.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import numpy as np
from geo import Grid, read_shp_polygons, read_dbf, rasterize_polys_fast
from lattice import Ctx, Agg, save_part, dist_m, dist_censored, coverage_share, polys_utm, read_tif, MIN_VALID
from paths import LC, GEE_EXPORTS, CITIES

WATER_BIG = {'River', 'Canal', 'Lake', 'Jheel', 'Creek', 'Mudflat'}
WATER_SMALL = {'Pond', 'Ditch', 'Aquaculture'}
DEV = {'Area Under Development'}
BUILT = {'Urban Built-up Area', 'Slum', 'Industrial Area', 'Commercial Area', 'Administrative area',
         'Educational Institutions', 'Terminal', 'Port', 'Air Port', 'Stadium', 'Brick-field',
         'Brick Field', 'Fort', 'Cattle market'}
GREEN = {'Green Space/ Open Space', 'Green Space/ Open space', 'Crop Land', 'Forest',
         'Rural Settlment With Homestead Vegetation', 'Graveyard'}


def lc_code(cl):
    if cl in WATER_BIG: return 1
    if cl in WATER_SMALL: return 2
    if cl in DEV: return 3
    if cl in BUILT: return 4
    if cl in GREEN: return 5
    if cl == 'Road': return 6
    return 7


def gsw_distance(x):
    """Distance to permanent water from the GEE GSW mask (1 = water), or None if not exported yet."""
    fn = GEE_EXPORTS / f'{x.c}_gsw_water.tif'
    if not fn.exists():
        x.log('no', fn.name, '- dist_water skipped until the GEE export is in place')
        return None
    a, tx, ty, sx, sy = read_tif(fn)
    water = np.nan_to_num(a, nan=0) > 0.5
    A = Agg(x.L, Grid(tx, ty, sx, a.shape[1], a.shape[0]))
    return A.mean(dist_m(water, sx))


def run(c):
    x = Ctx(c); g10, A = x.g10, x.A10
    stem = str(LC / f"{x.cfg['lc']}_City_LC_2024")
    polys = read_shp_polygons(stem + '.shp'); recs = read_dbf(stem + '.dbf')
    v0 = next(p for p in polys if len(p))[0][0]
    if abs(v0[0]) < 180:                                     # lon/lat shapefile -> UTM
        polys = [polys_utm(p) for p in polys]
    codes = [lc_code(r['Class']) for r in recs]
    order = np.argsort([-sum(len(r) for r in p) for p in polys])
    lc = rasterize_polys_fast(g10, [polys[i] for i in order], [codes[i] for i in order], fill=0, dtype=np.int8)
    cover = lc > 0
    F = {'lc_cover': A.mean(cover.astype(float))}
    for nm, code, r in [('built', 4, 3), ('green', 5, 3), ('pond', 2, 3), ('dev', 3, 3), ('water_big', 1, 5)]:
        F[f'{nm}_share_{r * 100}'] = coverage_share(A, lc == code, cover, r)
    for nm, code in [('dist_pond', 2), ('dist_dev', 3), ('dist_lc_water', 1)]:
        F[nm] = A.mean(dist_censored(lc == code, cover, 10.), MIN_VALID)
    dw = gsw_distance(x)
    if dw is not None:
        F['dist_water'] = dw
    save_part(c, 'landcover', F)
    ps = x.has_ps
    x.log('land cover', list(F), '| PS cells with a value:',
          {k: round(float(np.isfinite(v[ps]).mean()), 3) for k, v in F.items() if k != 'lc_cover'})


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
