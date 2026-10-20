"""Section 3.2  Conditioning factors - Geology (factor set M1).

Geological units of the 1:1M map of the Geological Survey of Bangladesh, as the polygons used in
the three-city paper (InSAR_data/factors_from_paper/geo_<c>.json), rasterised at 30 m.
Factors   geological unit (cell majority), Holocene share, signed distance to the contact between
          Holocene and older deposits (+ inside the Holocene). The distance is capped at +/-2 km, so it
          describes the contact zone; uncapped it reached 10 km and behaved like a coordinate
          (R2 0.8 against a quadratic trend in E, N). Pixels outside the map polygons are NaN.
          Shared classes for the transfer test are added in s32f_assemble_table.py.
Output    analysis/factors/parts/<c>_geology.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json, numpy as np
from geo import rasterize_polys_fast
from lattice import Ctx, save_part, dist_m, polys_utm, MIN_VALID
from paths import INSAR, CITIES

GEO_YOUNG = {'H2O', 'asc', 'asl', 'ava', 'ppc', 'csd', 'asd', 'dsl', 'dsd', 'lake'}


def run(c):
    x = Ctx(c); A = x.A30
    gj = json.load(open(INSAR / 'factors_from_paper' / f'geo_{c}.json'))
    units = sorted({g['glg'] for g in gj})
    ucode = {u: i + 1 for i, u in enumerate(units)}
    feats = sorted(gj, key=lambda g: -sum(len(r) for r in g['rings']))     # big first, small on top
    geo30 = rasterize_polys_fast(x.g30, [polys_utm(g['rings']) for g in feats], [ucode[g['glg']] for g in feats], fill=0)
    young = np.isin(geo30, [ucode[u] for u in units if u in GEO_YOUNG])
    mapped = geo30 > 0
    old = mapped & ~young
    dc = np.where(young, dist_m(old, 30.), -dist_m(young, 30.))
    F = {'geo_unit': A.majority(geo30, len(units) + 1),
         'd_contact': A.mean(np.where(mapped, dc, np.nan), MIN_VALID),
         'young_share': A.mean(np.where(mapped, young, np.nan), MIN_VALID)}
    save_part(c, 'geology', F, units=ucode); x.log('geology', units)


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
