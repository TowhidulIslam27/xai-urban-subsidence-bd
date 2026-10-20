"""Section 3.2  Groundwater proxy (M3, Dhaka only) - borehole / well locations.

Source    analysis/City_Maps/Borelogs_wellsLocation_Dha.gpkg (layer 'wells', points digitised by the
          authors; EPSG:32646 or 4326). The layer has an id only: no type, depth, aquifer or year.
Factors   dist_well      distance from the cell centre to the nearest point (m, capped at 2 km)
          well_dens_1k   points within 1 km of the cell centre, per km2
          well_kde_500   Gaussian kernel density, sigma 500 m (points per km2)
          All three are Dhaka only (TRANSFER = False: Chattogram has no equivalent layer).
          Cells farther than 2 km from every point lie outside the surveyed network; their values
          are kept (dist = 2 km, density 0) and the flag well_extent marks the network area.
Checks    Clark-Evans ratio of the point pattern inside the network (R < 1 clustered, as a
          production wellfield would be; R > 1 regular, as a borehole or monitoring grid would be),
          Spearman rho with vU and the trend R2 used by the screening, printed and written to
          analysis/screening/wells_dha_check.csv.
Output    analysis/factors/parts/dha_wells.npz; the columns are also added to
          analysis/factors/dha_factors.csv, dha_factors_full.csv and dha_maps.npz (by row/col),
          so the screening can run without rebuilding the whole table.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json, sqlite3
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from matplotlib.path import Path as MPath

from geo import gpkg_geom, ll2utm
from lattice import Ctx, save_part, DIST_CAP
from paths import CITY_MAPS, FACTORS, SCREEN

SRC = CITY_MAPS / 'Borelogs_wellsLocation_Dha.gpkg'
KEYS = ['dist_well', 'well_dens_1k', 'well_kde_500']


def read_points(fn):
    db = sqlite3.connect(str(fn))
    t, = db.execute("select table_name from gpkg_contents where data_type = 'features'").fetchone()
    col, srs = db.execute('select column_name, srs_id from gpkg_geometry_columns where table_name = ?',
                          (t,)).fetchone()
    xy = []
    for g, in db.execute(f'select "{col}" from "{t}"'):
        if g:
            _, parts = gpkg_geom(g)
            xy.append(np.asarray(parts[0], float).ravel()[:2])
    xy = np.array(xy)
    if srs == 4326:
        return np.c_[ll2utm(xy[:, 0], xy[:, 1])]
    assert srs == 32646, f'unexpected CRS {srs}'
    return xy


def clark_evans(P, area):
    d = cKDTree(P).query(P, k=2)[0][:, 1]
    return d.mean() / (0.5 / np.sqrt(len(P) / area)), d


def run():
    x = Ctx('dha')
    P = read_points(SRC)
    P = np.unique(np.round(P, 1), axis=0)
    tree = cKDTree(P)
    XY = np.c_[x.E.ravel(), x.N.ravel()]
    d, _ = tree.query(XY)
    dens = np.array([len(i) for i in tree.query_ball_point(XY, 1000.)]) / (np.pi * 1.0 ** 2)
    # Gaussian kernel, sigma 500 m, truncated at 3 sigma; per km2
    s = 500.
    kde = np.zeros(len(XY))
    for i, nb in enumerate(tree.query_ball_point(XY, 3 * s)):
        if nb:
            r2 = ((P[nb] - XY[i]) ** 2).sum(1)
            kde[i] = np.exp(-r2 / (2 * s * s)).sum() / (2 * np.pi * (s / 1000) ** 2)
    sh = x.E.shape
    F = {'dist_well': np.minimum(d, DIST_CAP).reshape(sh),
         'well_dens_1k': dens.reshape(sh),
         'well_kde_500': kde.reshape(sh)}
    ext = (d <= DIST_CAP).reshape(sh)
    save_part('dha', 'wells', F, n_points=int(len(P)), source=SRC.name)

    # ---- add the columns to the existing tables and maps
    add = pd.DataFrame(dict(row=x.R.ravel(), col=x.C.ravel(), well_extent=ext.ravel().astype(int),
                            **{k: F[k].ravel() for k in KEYS}))
    for fn in ['dha_factors.csv', 'dha_factors_full.csv']:
        f = FACTORS / fn
        if not f.exists():
            continue
        df = pd.read_csv(f, low_memory=False).drop(columns=KEYS + ['well_extent'], errors='ignore')
        n0 = len(df)
        df = df.merge(add, on=['row', 'col'], how='left')
        assert len(df) == n0 and df[KEYS].notna().all().all()
        df.to_csv(f, index=False)
        x.log('updated', fn)
    fm = FACTORS / 'dha_maps.npz'
    if fm.exists():
        z = np.load(fm, allow_pickle=False)
        maps = {k: z[k] for k in z.files}
        meta = json.loads(str(maps['meta']))
        assert maps['aoi'].shape == sh and meta['r1'] == x.L.r1 and meta['c0'] == x.L.c0
        maps.update({k: F[k].astype(np.float32) for k in KEYS}, well_extent=ext)
        np.savez_compressed(fm, **maps)
        x.log('updated dha_maps.npz')

    # ---- checks
    aoi, ps = x.aoi.ravel(), x.has_ps.ravel()
    inside = (d <= 1000.)
    # Clark-Evans inside the network: area = AOI cells within 1 km of a point
    area = float((aoi & inside).sum()) * 100 * 100
    R, nn = clark_evans(P, area)
    cells = pd.read_csv(FACTORS / 'dha_factors_full.csv', usecols=['row', 'col', 'has_ps', 'vU', 'E', 'N', 'in_city'])
    cells = cells.merge(add, on=['row', 'col'])
    p = cells[cells.has_ps.astype(bool)]
    from screening import trend_r2
    rows = []
    for k in KEYS:
        rho = spearmanr(p[k], p.vU, nan_policy='omit')[0]
        q = p[p.well_extent == 1]
        rho_in = spearmanr(q[k], q.vU, nan_policy='omit')[0]
        rows.append(dict(factor=k, rho_vU=rho, rho_vU_in_network=rho_in, trend_r2=trend_r2(p, k),
                         mode_share=p[k].value_counts(normalize=True).iloc[0]))
    chk = pd.DataFrame(rows)
    chk['n_points'] = len(P); chk['nn_median_m'] = np.median(nn); chk['clark_evans_R'] = R
    chk['ps_cells_in_network'] = float(p.well_extent.mean())
    chk['city_cells_in_network'] = float(cells[cells.in_city.astype(bool)].well_extent.mean())
    chk.to_csv(SCREEN / 'wells_dha_check.csv', index=False)
    print(chk.round(3).to_string(index=False))
    # vU near vs far
    for lo, hi in [(0, 250), (250, 500), (500, 1000), (1000, 1999), (1999, 1e9)]:
        m = (p.dist_well >= lo) & (p.dist_well < hi)
        print(f'dist {lo:>5.0f}-{hi:<6.0f} m  PS cells {m.sum():6d}  median vU {p.vU[m].median():6.2f}')


if __name__ == '__main__':
    run()
