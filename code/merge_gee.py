"""Section 3.2  Merge the Google Earth Engine exports into the cell tables.

usage: python merge_gee.py [folder]   (default analysis/gee/exports)
Joins on lattice row/col (exact), derives fill age and change bands, and writes
analysis/factors/<c>_factors_full.csv; the full-rectangle arrays of every band are added to
<c>_maps.npz for the maps. Missing themes are reported, not fatal.
The DEM and water-mask tifs in the same folder are read by s32a and s32c, not here."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import glob, json
import numpy as np, pandas as pd

from paths import FACTORS, GEE_EXPORTS
D = str(FACTORS) + os.sep
src = sys.argv[1] if len(sys.argv) > 1 else str(GEE_EXPORTS)
# A_ob matches A_ob1/A_ob2 (and per-band A_ob_b* files); B_fill / C_surface are older single exports
THEMES = ['A_ob', 'A_urban', 'B_gsw', 'B_ls', 'B_dw', 'B_fill', 'C_s2', 'C_lst', 'C_dw', 'C_surface', 'D_soil']


def read_export(f, meta):
    """One GEE CSV: check the row/col bands against the pixel centre in .geo, return row, col + bands."""
    g = pd.read_csv(f).drop(columns=['system:index'], errors='ignore')
    if '.geo' in g:
        xy = np.array([json.loads(s)['coordinates'] for s in g['.geo']], float)
        col = np.rint((xy[:, 0] - meta['a']) / meta['b']).astype(int)
        row = np.rint((xy[:, 1] - meta['c']) / meta['d']).astype(int)
        bad = float(((row != g.row) | (col != g.col)).mean())
        if bad > 0:
            print('   ', os.path.basename(f), f'row/col bands differ from the pixel centres in {bad:.1%} of rows;'
                  ' using the centres')
        g['row'], g['col'] = row, col
        g = g.drop(columns=['.geo'])
    return g.drop_duplicates(['row', 'col'])
for c in ['dha', 'ctg']:
    df = pd.read_csv(D + f'{c}_factors.csv', low_memory=False)
    z = np.load(D + f'{c}_maps.npz', allow_pickle=False)
    maps = {k: z[k] for k in z.files}
    meta = json.loads(str(z['meta']))
    for th in THEMES:
        fs = sorted(glob.glob(os.path.join(src, f'{c}_{th}*.csv')))
        if not fs:
            print(c, th, 'MISSING'); continue
        # one or several files per theme (split bands or split areas): combine on row/col,
        # filling gaps from later files and adding their new bands
        g = None
        for f in fs:
            h = read_export(f, meta).set_index(['row', 'col'])
            g = h if g is None else g.combine_first(h)
        g = g.reset_index()
        bands = [k for k in g.columns if k not in ('row', 'col')]
        n0 = len(df)
        df = df.merge(g, on=['row', 'col'], how='left', suffixes=('', '_' + th))
        assert len(df) == n0
        print(c, th, len(bands), 'bands | cells with a value (first band, %s):' % bands[0],
              round(float(df[bands[0]].notna().mean()), 3))
        # the exports cover the whole lattice rectangle: keep full 2-D arrays for the maps
        for k in bands:
            v = pd.to_numeric(g[k], errors='coerce').to_numpy(float)
            a = np.full((meta['nrows'], meta['ncols']), np.nan)
            i = (meta['r1'] - g.row).to_numpy(int); j = (g.col - meta['c0']).to_numpy(int)
            ok = (i >= 0) & (i < meta['nrows']) & (j >= 0) & (j < meta['ncols'])
            a[i[ok], j[ok]] = v[ok]
            maps[k] = a.astype(np.float32)
    if 'gsw_fill_year' in df:
        df['fill_age_2023'] = 2023 - df.gsw_fill_year
    if 'ob_bfrac_2023' in df and 'ob_bfrac_2016' in df:
        df['ob_bfrac_change'] = df.ob_bfrac_2023 - df.ob_bfrac_2016
    if 'dw_built_2023' in df and 'dw_built_2018' in df:
        df['dw_built_change'] = df.dw_built_2023 - df.dw_built_2018
    if 'water_share_1990' in df and 'water_share_2022' in df:   # same Landsat MNDWI method in both years
        df['water_loss_1990_2022'] = df.water_share_1990 - df.water_share_2022
    for k in ('l5_nobs_1990', 'l5_nobs_2005', 'l8_nobs_2022'):
        if k in df:
            q = df.loc[df.has_ps.astype(bool), k]
            print(c, k, 'clear dry-season scenes per pixel at PS cells: median', q.median(), '| < 5:', round(float((q < 5).mean()), 3))
    for new, (x, y, sign) in {'fill_age_2023': ('gsw_fill_year', None, 0), 'ob_bfrac_change': ('ob_bfrac_2023', 'ob_bfrac_2016', -1),
                              'dw_built_change': ('dw_built_2023', 'dw_built_2018', -1),
                              'water_loss_1990_2022': ('water_share_1990', 'water_share_2022', -1)}.items():
        if x in maps and (y is None or y in maps):
            maps[new] = (2023 - maps[x]) if y is None else (maps[x] + sign * maps[y])
    np.savez_compressed(D + f'{c}_maps.npz', **maps)
    df.to_csv(D + f'{c}_factors_full.csv', index=False)
    print(c, 'written', df.shape)
