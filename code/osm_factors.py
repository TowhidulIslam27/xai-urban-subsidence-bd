"""Section 3.2  Conditioning factors - Urban load and access from OpenStreetMap (M2).

Source    analysis/inputs/<c>_infra_lines.sqlite (OSM extract of the three-city paper), 10 m.
Factors   distance to main road, road density within 500 m (km/km2), distance to railway,
          density of schools/health/power facilities (per km2, 500 m). All distances are capped at
          2 km (DIST_CAP in common/lattice.py).
          City-specific (not used in the transfer test): distance to MRT-6 and the elevated
          expressway (Dhaka), to the port (CTG), and distance to embankments. The OSM embankment
          layer is not the same object in the two cities: in Chattogram 18 of 23 features are the
          City Outer Ring Road (a coastal embankment-road built around the InSAR window; dates to
          confirm) plus a stretch of coastal road;
          in Dhaka it is a few beribadh/embankment roads (the flood embankment system is largely
          missing). Features named as railways are left out.
Output    analysis/factors/parts/<c>_osm.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import sqlite3, numpy as np
import hydro
from geo import ll2utm, rasterize_lines, rasterize_polys_fast, read_gpkg_table
from lattice import Ctx, save_part, dist_m, lines_utm, polys_utm
from paths import INPUTS, CITIES


def run(c):
    x = Ctx(c); g10, A = x.g10, x.A10
    db = str(INPUTS / f'{c}_infra_lines.sqlite')
    tabs = {t for t, in sqlite3.connect(db).execute("select name from sqlite_master where type='table'")}
    main = rasterize_lines(g10, lines_utm(db, 'roads', lambda r: r[0] == 'main', ('road_class',)))
    allr = rasterize_lines(g10, lines_utm(db, 'roads'))
    rail = rasterize_lines(g10, lines_utm(db, 'railways', lambda r: r[0] in ('rail', 'narrow_gauge'), ('railway',)))
    F = {'dist_mainroad': A.mean(dist_m(main, 10.)),
         # burnt 10 m pixels x 10 m = road length; per 100 m cell (0.01 km2), mean over 500 m
         'road_dens_500': np.clip(hydro.focal_mean(A.mean(allr.astype(float)) * 100 * 10 / 1000 / 0.01, 5), 0, None),
         'dist_rail': A.mean(dist_m(rail, 10.))}
    if 'embankments' in tabs:
        emb = lines_utm(db, 'embankments', lambda r: 'railway' not in (r[0] or '').lower(), ('name',))
        F['dist_embank'] = A.mean(dist_m(rasterize_lines(g10, emb), 10.))
    if 'mrt6' in tabs:
        F['dist_mrt'] = A.mean(dist_m(rasterize_lines(g10, lines_utm(db, 'mrt6')), 10.))
    if 'expressway' in tabs:
        F['dist_expressway'] = A.mean(dist_m(rasterize_lines(g10, lines_utm(db, 'expressway')), 10.))
    if 'port' in tabs:
        rows = read_gpkg_table(db, 'port')
        pm = np.zeros((g10.ny, g10.nx), bool)
        pl = [polys_utm(r[0][1]) for r in rows if r[0][0] == 'pg']
        if pl:
            pm |= rasterize_polys_fast(g10, pl, [1] * len(pl), fill=0, dtype=np.int8) > 0
        pm |= rasterize_lines(g10, [q for r in rows for q in polys_utm(r[0][1])])
        F['dist_port'] = A.mean(dist_m(pm, 10.))
    fac = read_gpkg_table(db, 'facilities', ('group',))
    # points as they are; polygons at the mean of their outer-ring vertices
    pts = np.array([r[0][1][0][:, :2].mean(0) if r[0][0] == 'pg' else r[0][1][0][0, :2]
                    for r in fac if r[0][0] in ('pt', 'pg')])
    if len(pts):
        e, n = ll2utm(pts[:, 0], pts[:, 1]); i, j = g10.ij(e, n)
        k = (i >= 0) & (i < g10.ny) & (j >= 0) & (j < g10.nx)
        fm = np.zeros((g10.ny, g10.nx), float); np.add.at(fm, (i[k], j[k]), 1)
        F['facility_dens_500'] = np.clip(hydro.focal_mean(A.mean(fm) * A.cnt.reshape(A.shape), 5) * 100, 0, None)
    save_part(c, 'osm', F); x.log('osm', list(F))


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
