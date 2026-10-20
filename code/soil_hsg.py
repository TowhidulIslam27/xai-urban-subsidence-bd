"""Section 3.2  Conditioning factors - Hydrologic soil group (M0), both cities.

Source    HYSOGs250m (Ross et al. 2018), global 250 m (1/480 deg) hydrologic soil groups, from
          Dataset/HYSOGs250m.tif. The two city windows were cut out beforehand (the file is
          570 MB) into analysis/inputs/hysogs/<c>_hysogs.npy + .json (upper-left corner, pixel size).
Classes   1 A, 2 B, 3 C, 4 D; dual classes 11 A/D, 12 B/D, 13 C/D, 14 D/D (high water table: the
          first letter applies when drained, D when not); 255 = no data.
Factors   HSG pixels are larger than the 100 m cells, so each cell takes its nearest HSG pixel and
          the factors are shares within 300 m (as for the land-cover shares):
            hsg_dual_share_300   share of dual (high water table) soils
            hsg_poor_share_300   share of poorly drained soils (D, C/D, D/D)
          hsg_class (nearest pixel) is kept for checks. No-data pixels are left out of the shares.
Output    analysis/factors/parts/<c>_hsg.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json, numpy as np
import hydro
from lattice import Ctx, save_part
from paths import INPUTS, CITIES


def run(c):
    x = Ctx(c)
    src = INPUTS / 'hysogs' / f'{c}_hysogs.npy'
    if not src.exists():
        x.log('no', src.name, '- cut the HYSOGs window first'); return
    a = np.load(src); m = json.load(open(INPUTS / 'hysogs' / f'{c}_hysogs.json'))
    j = np.floor((x.LON - m['lon0']) / m['dx']).astype(int)
    i = np.floor((m['lat0'] - x.LAT) / m['dy']).astype(int)
    ok = (i >= 0) & (i < a.shape[0]) & (j >= 0) & (j < a.shape[1])
    cls = np.full(x.LON.shape, 255, np.uint8); cls[ok] = a[i[ok], j[ok]]
    valid = cls != 255
    dual = np.isin(cls, [11, 12, 13, 14]); poor = np.isin(cls, [4, 13, 14])
    F = {'hsg_class': np.where(valid, cls, np.nan)}
    for nm, msk in [('hsg_dual_share_300', dual), ('hsg_poor_share_300', poor)]:
        num = hydro.focal_mean(np.where(valid, msk, 0).astype(float), 3)
        den = hydro.focal_mean(valid.astype(float), 3)
        with np.errstate(invalid='ignore', divide='ignore'):
            F[nm] = np.where(den >= 0.5, num / den, np.nan)
    save_part(c, 'hsg', F)
    ps = x.has_ps
    v, n = np.unique(cls[ps], return_counts=True)
    x.log('HSG classes at PS cells', dict(zip(v.tolist(), n.tolist())), '| shares valid at PS cells',
          round(float(np.isfinite(F['hsg_dual_share_300'][ps]).mean()), 3))


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
