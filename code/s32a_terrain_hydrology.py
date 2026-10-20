"""Section 3.2  Conditioning factors - Terrain and hydrology (factor sets M0).

DEM choice (automatic, logged and stored in analysis/factors/<c>_sources.json):
  1. FABDEM bare earth from the GEE export (analysis/gee/exports/<c>_fabdem.tif), if present.
     Ground = FABDEM (3x3 mean); surface = GLO-30 (<c>_glo30.tif) or SRTM, so
     dsm_minus_ground is roughly the 2011-2015 building and tree height.
  2. Otherwise SRTM 1 arc-second (Feb 2000) with a ground proxy: 20th percentile in a 270 m
     window, then a 3x3 mean. PROVISIONAL: on the Chattogram hills this proxy cuts the ridges
     (slope at PS cells is ~0.2-0.3 x the GLO-30 slope; on slopes >= 10 deg the ground is ~7 m
     lower, partly tree canopy in GLO-30), so terrain from SRTM is not used for the final models.
  Override with DEM_<C>=path (bare earth if DEM_BARE=1, e.g. another DEM).

Terrain   ground elevation, DSM minus ground, slope, aspect (sin, cos), profile and plan
          curvature (Zevenbergen-Thorne), TRI (Riley 1999), TPI 300 m and 1 km, relief 1 km,
          elevation-relief ratio (450 m).
Hydrology priority-flood fill with epsilon gradient, D8 receivers, flow accumulation;
          TWI = ln(a / tan b); streams = contributing area >= 0.5 km2; HAND; depression depth;
          distance to stream; drainage density (1 km disk).
Output    analysis/factors/parts/<c>_terrain.npz, <c>_hydrology.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import map_coordinates
import hydro
from geo import utm2ll, read_srtm, sample_regular
from lattice import Ctx, save_part, dist_m, read_tif
from paths import DEMS, GEE_EXPORTS, CITIES


def utm_tif_on_grid(fn, x):
    """A 30 m UTM 46N GeoTIFF (GEE export) bilinearly sampled on the work grid."""
    X, Y = x.g30.centers()
    a, tx, ty, sx, sy = read_tif(str(fn))
    a = np.nan_to_num(a, nan=np.nanmedian(a))
    return map_coordinates(a, [(ty - Y) / sy - .5, (X - tx) / sx - .5], order=1, mode='nearest').astype(np.float64)


def srtm_on_grid(x):
    X, Y = x.g30.centers()
    lon30, lat30 = utm2ll(X, Y)
    a, l0, la0, dx, dy = read_srtm(str(DEMS / f"{x.cfg['dem']}_1arc_v3.tif"))
    return sample_regular(np.nan_to_num(a, nan=np.nanmedian(a)), l0, la0, dx, dy, lon30, lat30).astype(np.float64)


def ground_surface(x):
    """Returns surface model, ground, and a short record of the sources."""
    own = os.environ.get(f'DEM_{x.c.upper()}')
    fab, glo = GEE_EXPORTS / f'{x.c}_fabdem.tif', GEE_EXPORTS / f'{x.c}_glo30.tif'
    if own:
        z = utm_tif_on_grid(own, x)
        if os.environ.get('DEM_BARE') == '1':
            dsm, ground, src = srtm_on_grid(x), ndi.uniform_filter(z, 3, mode='nearest'), dict(ground=own, surface='SRTM')
        else:
            dsm = z
            ground = ndi.uniform_filter(ndi.percentile_filter(z, 20, size=9, mode='nearest'), 3, mode='nearest')
            src = dict(ground=f'{own}, 20th percentile proxy', surface=own)
    elif fab.exists():
        ground = ndi.uniform_filter(utm_tif_on_grid(fab, x), 3, mode='nearest')
        dsm = utm_tif_on_grid(glo, x) if glo.exists() else srtm_on_grid(x)
        src = dict(ground='FABDEM (bare earth)', surface='GLO-30' if glo.exists() else 'SRTM (Feb 2000)')
    else:
        dsm = srtm_on_grid(x)
        ground = ndi.uniform_filter(ndi.percentile_filter(dsm, 20, size=9, mode='nearest'), 3, mode='nearest')
        src = dict(ground='SRTM 20th-percentile proxy (PROVISIONAL)', surface='SRTM (Feb 2000)')
    x.log('terrain DEM:', src)
    return dsm, ground, src


def run(c):
    x = Ctx(c); A = x.A30
    dsm, ground, src = ground_surface(x)
    T = hydro.terrain(ground, 30.)
    F = {'elev': A.mean(ground), 'dsm_minus_ground': A.mean(dsm - ground)}
    for k in ('slope', 'asp_sin', 'asp_cos', 'curv_prof', 'curv_plan', 'tri'):
        F[k] = A.mean(T[k])
    F['tpi_300'] = A.mean(ground - hydro.focal_mean(ground, 10))
    F['tpi_1000'] = A.mean(ground - hydro.focal_mean(ground, 33))
    F['relief_1k'] = A.mean(ndi.maximum_filter(ground, 33) - ndi.minimum_filter(ground, 33))
    F['err_450'] = A.mean(hydro.err(ground, 15))
    save_part(c, 'terrain', F, dem=src); x.log('terrain', list(F))

    zf = hydro.fill_epsilon(ground)
    rec, _ = hydro.d8_receivers(zf, 30.)
    acc, order = hydro.accumulate(zf, rec)
    tanb = np.maximum(np.tan(np.radians(T['slope'])), 1e-3)
    stream = acc * 900 >= 0.5e6
    H = {'twi': A.mean(np.log(acc * 30. / tanb)),
         'fill_depth': A.mean(zf - ground),
         # HAND >= 0: a pit below its stream outlet is at the drainage level, not under it
         'hand': A.mean(np.maximum(hydro.hand(zf, ground, rec, stream, order), 0)),
         'dist_stream': A.mean(dist_m(stream, 30.)),
         'drain_dens_1k': A.mean(np.clip(hydro.focal_mean(stream.astype(float), 17), 0, None) * 1e6 / 900 * 30 / 1e3)}
    save_part(c, 'hydrology', H); x.log('hydrology', list(H), '| stream px', int(stream.sum()))


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
