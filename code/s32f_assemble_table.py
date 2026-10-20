"""Section 3.2 / 4.1  Assemble the cell table for modelling.

Joins all factor parts (terrain, hydrology, geology, landcover, osm, ctgkit) with the target
(vU, vU_se, vU_h2, ...), the building factors of the three-city paper (PS cells only, suffix _ps),
shared geology classes, 1 km blocks for spatial CV and the vU classes (-2 / -5 mm/yr).
If analysis/gee/exports holds the GEE CSVs, run s32h_merge_gee.py afterwards.
Output    analysis/factors/<c>_factors.csv   one row per prediction cell (has_ps = vU available)
          analysis/factors/<c>_maps.npz      2-D factor arrays for the figures
          analysis/factors/<c>_sources.json  which DEM / water source the run used
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json, numpy as np, pandas as pd
from lattice import Ctx, load_part
from paths import FACTORS, INSAR, CITIES

PARTS = ['terrain', 'hydrology', 'geology', 'landcover', 'osm', 'hsg', 'ctgkit', 'wells']
GEO_SHARED = {'H2O': 'channel_1990', 'asc': 'holocene_fine', 'asl': 'holocene_fine', 'ava': 'holocene_fine',
              'ppc': 'holocene_marsh', 'csd': 'holocene_sand', 'rm': 'older_consolidated',
              'QTdi': 'older_consolidated', 'QTdt': 'older_consolidated', 'QTdd': 'older_consolidated',
              'Tt': 'older_consolidated', 'Tbb': 'older_consolidated', 'lake': 'water'}


def run(c):
    x = Ctx(c); L = x.L
    F, ucode, src = {}, {}, {}
    for p in PARTS:
        part, extra = load_part(c, p)
        if part is None:
            if p not in ('ctgkit', 'wells'):
                x.log('missing part', p, '- run its script first')
            continue
        F.update(part); ucode.update(extra.get('units', {}))
        if 'dem' in extra:
            src['terrain'] = extra['dem']
    src['dist_water'] = 'JRC GSW occurrence >= 50 % (GEE)' if 'dist_water' in F else 'pending GEE export'
    json.dump(src, open(FACTORS / f'{c}_sources.json', 'w'), indent=1)
    df = pd.DataFrame(dict(row=x.R.ravel(), col=x.C.ravel(), lon=x.LON.ravel(), lat=x.LAT.ravel(),
                           E=x.E.ravel(), N=x.N.ravel(), has_ps=x.has_ps.ravel(),
                           near_ps=x.near_ps.ravel(), in_city=x.in_city.ravel()))
    for k, v in F.items():
        df[k] = np.asarray(v, float).ravel()
    df = df[x.aoi.ravel()].copy()
    df['cell_id'] = (df.row * 100000 + df.col).astype(np.int64)        # same id as the three-city paper
    df['block_1km'] = (df.row // 10).astype(str) + '_' + (df.col // 10).astype(str)
    df['geo_glg'] = df.geo_unit.map({float(v): k for k, v in ucode.items()})
    df['geo_shared'] = df.geo_glg.map(GEO_SHARED)
    tgt = x.cells[['cell_id', 'n_asc', 'n_dsc', 'vU', 'vU_se', 'vU_h2', 'vU_h2_se', 'vE_dr', 'sigU', 'ampU',
                   'gi_z', 'gi_class', 'thana']]
    df = df.merge(tgt, on='cell_id', how='left')
    pap = pd.read_csv(INSAR / 'factors_from_paper' / f'construction_cells_{c}.csv')
    pap = pap[['cell_id', 'glg', 'density', 'n_bld', 'height', 'growth', 'new_n', 'old_n', 'h23_ob', 'eff_late']]
    pap = pap.rename(columns={'glg': 'glg_paper', 'density': 'bld_density_ps', 'n_bld': 'n_bld_ps',
                              'height': 'bld_height_ps', 'growth': 'bld_growth_ps', 'h23_ob': 'h23_ob_ps',
                              'eff_late': 'eff_late_ps'})
    pap['new_share_ps'] = pap.new_n / (pap.new_n + pap.old_n).replace(0, np.nan)
    df = df.merge(pap.drop(columns=['new_n', 'old_n']), on='cell_id', how='left')
    df['class_vU'] = pd.cut(df.vU, [-np.inf, -5, -2, np.inf], labels=['high', 'moderate', 'stable'])
    df['city'] = c
    df.to_csv(FACTORS / f'{c}_factors{os.environ.get("TAG", "")}.csv', index=False)
    lat_, _ = load_part(c, 'lattice')
    np.savez_compressed(FACTORS / f'{c}_maps.npz', aoi=x.aoi, has_ps=x.has_ps, near_ps=x.near_ps, in_city=x.in_city,
                        **{k: np.asarray(v, np.float32) for k, v in F.items()},
                        vU=lat_['vU'].astype(np.float32), meta=json.dumps(L.meta()), units=json.dumps(ucode))
    chk = df[df.has_ps]
    x.log('table', df.shape, '| PS rows', len(chk), '| geology agrees with paper',
          round(float((chk.geo_glg == chk.glg_paper).mean()), 3))


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
