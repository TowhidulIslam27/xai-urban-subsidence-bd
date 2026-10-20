"""Section 3.1  City boundaries by dissolving the 2024 city land-cover maps.

The land-cover maps of the two city corporations (analysis/inputs/3_City_LC_2024) tile the city
area with ~12 000 (Dhaka) and ~17 000 (Chattogram) polygons. All polygons are dissolved into one
outline per city: rasterised at 10 m, small slivers between polygons closed (30 m), interior holes
filled, and parts smaller than 1 km2 dropped. The outline is traced at 10 m and written in lon/lat.
The boundary is used to clip the factor maps for display; the models still use the whole
prediction area (cells within 1 km of a PS cell), which reaches beyond the city.
Output    analysis/factors/<c>_city_boundary.geojson
          analysis/factors/city_boundary_preview.png  (quick look, not a manuscript figure)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json
import numpy as np
from scipy import ndimage as ndi
from skimage import measure
from geo import read_shp_polygons, rasterize_polys_fast, utm2ll
from lattice import Ctx, polys_utm
from paths import LC, FACTORS, CITIES

CLOSE_M, MIN_KM2 = 30., 1.0


def dissolve(c):
    x = Ctx(c); g = x.g10
    polys = read_shp_polygons(str(LC / f"{x.cfg['lc']}_City_LC_2024") + '.shp')
    polys = [p for p in polys if len(p)]
    if abs(polys[0][0][0][0]) < 180:                         # lon/lat shapefile -> UTM
        polys = [polys_utm(p) for p in polys]
    m = rasterize_polys_fast(g, polys, [1] * len(polys), dtype=np.int8) > 0
    r = int(CLOSE_M / g.res)
    k = np.hypot(*np.ogrid[-r:r + 1, -r:r + 1]) <= r
    m = ndi.binary_closing(np.pad(m, r + 1), k)[r + 1:-r - 1, r + 1:-r - 1]
    m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    area = ndi.sum(m, lab, range(1, n + 1)) * g.res ** 2 / 1e6
    keep = np.isin(lab, 1 + np.flatnonzero(area >= MIN_KM2))
    rings = []
    for cnt in measure.find_contours(np.pad(keep, 1).astype(float), 0.5):
        if len(cnt) < 20:
            continue
        i, j = cnt[:, 0] - 1, cnt[:, 1] - 1                  # back to unpadded pixel index
        E = g.x0 + (j + .5) * g.res; N = g.y0 - (i + .5) * g.res
        lon, lat = utm2ll(E, N)
        rings.append(np.round(np.c_[lon, lat], 6).tolist())
    gj = {'type': 'FeatureCollection', 'features': [{'type': 'Feature', 'properties': {
        'city': c, 'source': 'dissolved 2024 city land-cover map', 'area_km2': round(float(keep.sum()) * g.res ** 2 / 1e6, 1)},
        'geometry': {'type': 'MultiPolygon', 'coordinates': [[r] for r in rings]}}]}
    json.dump(gj, open(FACTORS / f'{c}_city_boundary.geojson', 'w'))
    x.log('city boundary:', len(rings), 'ring(s),', gj['features'][0]['properties']['area_km2'], 'km2')
    return rings


def preview(R):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(7, 4.5))
    for ax, c in zip(axs, CITIES):
        for r in R[c]:
            r = np.asarray(r); ax.fill(r[:, 0], r[:, 1], fc='#3b3b3b', ec='k', lw=0.6)
        ax.set_aspect(1 / np.cos(np.radians(np.mean([np.mean(np.asarray(r)[:, 1]) for r in R[c]]))))
        ax.set_title({'dha': 'Dhaka', 'ctg': 'Chattogram'}[c], fontsize=9); ax.tick_params(labelsize=6)
    fig.tight_layout(); fig.savefig(FACTORS / 'city_boundary_preview.png', dpi=200)


if __name__ == '__main__':
    R = {c: dissolve(c) for c in (sys.argv[1:] or CITIES)}
    if len(R) == 2:
        preview(R)
