"""Shared figure code for the factor figures: labels, colour maps, lattice maps (Dhaka beside
Chattogram, equal-size boxes at one map scale, north arrow and scale bar on every panel).
Maps sit on a dark grey canvas (land, water, OSM roads, thana lines) and show each factor over
the prediction area, the cells the models use (within 1 km of a PS cell); water cells are left
to the canvas and cells without source data are light grey. All figures are written as PNG at 600 dpi to
manuscript/Figures."""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import TwoSlopeNorm, Normalize
from scipy.stats import spearmanr
import figstyle as fs
from factor_groups import PRIORITY, GROUP, LC_CHECK, PS_ONLY

import os
from paths import FACTORS, FIGS, SCREEN, INSAR, INPUTS, LC, GEE_EXPORTS, CITY_MAPS, CITIES
D, FIG, O = str(FACTORS) + os.sep, str(FIGS), str(SCREEN) + os.sep
fs.BASE = str(INSAR / 'basedata')     # thana boundaries
fs.set_style()
M = T = None
MM = fs.MM

LABEL = dict(
    elev='Ground elevation (m)', dsm_minus_ground='DSM minus ground (m)', slope='Slope (°)',
    asp_sin='Aspect, eastness', asp_cos='Aspect, northness', curv_prof='Profile curvature (1/m)',
    curv_plan='Plan curvature (1/m)', tri='TRI (m)', tpi_300='TPI 300 m (m)', tpi_1000='TPI 1 km (m)',
    relief_1k='Relief 1 km (m)', err_450='Elevation-relief ratio', twi='TWI', hand='HAND (m)',
    fill_depth='Depression depth (m)', dist_stream='Distance to stream (km)', drain_dens_1k='Drainage density (km/km²)',
    dist_river='Distance to river (km)', water_big_share_500='Open water share, 500 m',
    d_contact='Distance to Holocene contact (km)', young_share='Holocene share',
    built_share_300='Built share, 300 m', green_share_300='Green share, 300 m',
    pond_share_300='Pond/ditch share, 300 m', dev_share_300='Land under development, 300 m',
    dist_pond='Distance to pond (km)', dist_dev='Distance to development (km)',
    road_dens_500='Road density (km/km²)', dist_mainroad='Distance to main road (km)',
    dist_rail='Distance to railway (km)', dist_embank='Distance to embankment (km)',
    facility_dens_500='Facilities (per km²)',
    bld_density_ps='Building density', bld_height_ps='Building height (m)', h23_ob_ps='Height 2023 (m)',
    new_share_ps='New-construction share', bld_growth_ps='Built-load growth',
    dist_water='Distance to permanent water (km)', dist_lc_water='Distance to mapped water, LC (km)',
    dist_mrt='Distance to MRT-6 (km)', dist_expressway='Distance to expressway (km)',
    dist_port='Distance to port (km)', lc_cover='Land-cover map coverage',
    hsg_dual_share_300='High-water-table soils, 300 m', hsg_poor_share_300='Poorly drained soils, 300 m',
    ob_bfrac_2023='Building footprint 2023', ob_height_2023='Building height 2023 (m)',
    ob_volume_2023='Building volume (m³/m²)', ob_new_frac='New footprint 2016-23',
    ob_dheight_existing='Height change, existing (m)', ob_bld_per_ha='Buildings per ha',
    ghs_bv_2020='GHSL built volume 2020', ghs_bh_2018='GHSL height 2018 (m)',
    ghs_bs_d2000_2020='GHSL built surface gain 2000-20', ghs_urban_epoch='First urban epoch',
    pop_density_2020='Population density (per km²)', viirs_rad='Night lights (nW/cm²/sr)',
    gsw_lost_share='Water lost (share)', gsw_lost_window='Water lost 2017-19 (share)',
    gsw_lost_after2010='Water lost after 2010 (share)', fill_age_2023='Fill age (yr)',
    water_share_1990='Water share 1990', water_loss_1990_2022='Water loss 1990-2022',
    dw_built_change='Built probability change', dw_built='DW built', dw_trees='DW trees',
    dw_grass='DW grass', dw_crops='DW crops', dw_bare='DW bare', dw_water='DW water',
    dw_flooded_vegetation='DW flooded vegetation', s2_ndvi='NDVI', s2_ndbi='NDBI', s2_mndwi='MNDWI',
    s2_bsi='BSI', s2_ndvi_sd='NDVI s.d.', merit_hand='MERIT HAND (m)', olm_clay='Clay 0-200 cm (%)',
    olm_sand='Sand 0-200 cm (%)', olm_bulkdens='Bulk density (kg/m³)', olm_orgc='Organic carbon (g/kg)')
KM = {k for k in LABEL if k.startswith('dist_') or k == 'd_contact'}
CMAP = {'Terrain': 'YlOrBr', 'Hydrology': 'Blues', 'Geology': 'BrBG', 'Fill': 'Purples', 'Urban load': 'Reds',
        'Urban load (PS only)': 'Reds', 'City-specific': 'Oranges', 'Land cover': 'Greens',
        'Land cover (check)': 'Greens', 'Soil': 'copper_r', 'Groundwater proxy': 'PuRd'}
# copper_r ends in black, which reads as a hole on the dark canvas: use its lighter three quarters
_cu = plt.get_cmap('copper_r')
CMAP['Soil'] = plt.matplotlib.colors.LinearSegmentedColormap.from_list('soil', _cu(np.linspace(0, 0.72, 256)))
DIVERGING = {'d_contact', 'tpi_300', 'tpi_1000', 'asp_sin', 'asp_cos', 'curv_prof', 'curv_plan'}


def init():
    global M, T
    M, T = load()


def save(fig, name):
    fs.save(fig, FIG, name, formats=('png',))


def load():
    M, T = {}, {}
    for c in CITIES:
        z = np.load(D + f'{c}_maps.npz', allow_pickle=True)
        M[c] = {k: z[k] for k in z.files}
        M[c]['meta'] = json.loads(str(z['meta']))
        f = D + f'{c}_factors_full.csv'
        T[c] = pd.read_csv(f if os.path.exists(f) else D + f'{c}_factors.csv', low_memory=False)
    return M, T


# ------------------------------------------------------------------ dark basemap
LAND, WATER, ROAD, MAINROAD, BOUND = '#3b3b3b', '#1d2024', '#535353', '#6b6b6b', '#8c8c8c'
NODATA = '#d2d2d2'                       # land inside the city rectangle without source data
SEA_LABEL = {'ctg': ('Bay of Bengal', 91.712, 22.30, 90)}
_BM = {}


CITY_FILES = {'dha': ('Dhaka_city_bound.gpkg', 'Dhaka_city_wards.gpkg'),
              'ctg': ('Chattogram_city_boundary.gpkg', 'Chattogram_city_wards.gpkg')}
_CITY = {}


def _gpkg_rings(fn):
    """Rings of every polygon in the first non-empty feature table of a GeoPackage (lon/lat)."""
    import sqlite3
    from geo import gpkg_geom
    db = sqlite3.connect(str(fn))
    for t, in db.execute("select table_name from gpkg_contents where data_type = 'features'"):
        col = db.execute("select column_name from gpkg_geometry_columns where table_name = ?", (t,)).fetchone()[0]
        rows = [r[0] for r in db.execute(f'select "{col}" from "{t}"') if r[0]]
        if rows:
            return [[p[:, :2] for p in gpkg_geom(b)[1]] for b in rows]
    return []


def city(c):
    """City corporation outline (list of rings) and ward polygons from analysis/City_Maps."""
    if c not in _CITY:
        fb, fw = CITY_FILES[c]
        _CITY[c] = dict(outline=[r for f in _gpkg_rings(CITY_MAPS / fb) for r in f],
                        wards=_gpkg_rings(CITY_MAPS / fw))
    return _CITY[c]


def city_path(c):
    from matplotlib.path import Path
    rings = city(c)['outline']
    verts = np.vstack([np.r_[r, r[:1]] for r in rings])
    codes = np.concatenate([[Path.MOVETO] + [Path.LINETO] * (len(r) - 1) + [Path.CLOSEPOLY] for r in rings])
    return Path(verts, codes)


def city_cells(c):
    """Lattice cells whose centre lies inside the city outline."""
    m = M[c]['meta']
    R, C = np.meshgrid(np.arange(m['r1'], m['r0'] - 1, -1), np.arange(m['c0'], m['c1'] + 1), indexing='ij')
    pts = np.c_[(m['a'] + m['b'] * C).ravel(), (m['c'] + m['d'] * R).ravel()]
    return city_path(c).contains_points(pts).reshape(R.shape)


def basemap(c, half, clip=False):
    """Land/water raster over the panel box (lon/lat, ~0.0003 deg) and OSM road lines.
    Water = permanent water from the GSW export when present (full coverage); otherwise open
    water of the 2024 city land-cover map plus, where thana polygons leave gaps (Chattogram:
    sea, Karnaphuli, Halda), the gaps. Dhaka's thanas include the rivers, so before the GSW
    export Dhaka shows only the mapped rivers, canals and lakes."""
    key = (c, tuple(half), clip)
    if key in _BM:
        return _BM[key]
    from geo import Grid, rasterize_polys_fast, read_shp_polygons, read_dbf, read_gpkg_table, ll2utm
    from lattice import CITY, read_tif
    (w, e, s_, n), lat0 = box(c, half, clip)
    res = 0.0003
    g = Grid(w, n, res, int((e - w) / res) + 1, int((n - s_) / res) + 1)
    th = fs.thanas(c)
    land = rasterize_polys_fast(g, [r for _, r in th], [1] * len(th), dtype=np.int8) > 0
    water = ~land
    fn = GEE_EXPORTS / f'{c}_gsw_water.tif'
    if fn.exists():
        a, tx, ty, sx, sy = read_tif(fn)
        X, Y = g.centers(); E, N = ll2utm(X, Y)
        i = ((ty - N) / sy).astype(int); j = ((E - tx) / sx).astype(int)
        ok = (i >= 0) & (i < a.shape[0]) & (j >= 0) & (j < a.shape[1])
        gw = np.zeros(water.shape, bool); gw[ok] = np.nan_to_num(a[i[ok], j[ok]]) > 0.5
        water |= gw
        canal = np.zeros(water.shape, bool)
        src = 'GSW'
    else:
        stem = str(LC / f"{CITY[c]['lc']}_City_LC_2024")
        recs = read_dbf(stem + '.dbf')
        allp = read_shp_polygons(stem + '.shp')
        if abs(next(p for p in allp if len(p))[0][0][0]) > 180:         # UTM shapefile -> lon/lat
            from geo import utm2ll
            allp = [[np.c_[utm2ll(r[:, 0], r[:, 1])] for r in p] for p in allp]
        # canals are drawn on the basemap but do not blank lattice cells (thin lines would speckle the data)
        cls = {'River': 2, 'Lake': 2, 'Jheel': 2, 'Creek': 2, 'Mudflat': 2, 'Canal': 1}
        sel = [(p, cls[r['Class']]) for p, r in zip(allp, recs) if r['Class'] in cls and len(p)]
        lcw = rasterize_polys_fast(g, [p for p, _ in sel], [v for _, v in sel], dtype=np.int8)
        canal = lcw == 1
        water |= lcw == 2
        src = 'thana gaps + city LC water'
    ext = [w, w + g.nx * res, n - g.ny * res, n]
    drawn = water | canal if src != 'GSW' else water
    roads = {'main': [], 'other': []}
    db = str(INPUTS / f'{c}_infra_lines.sqlite')
    for geom, rc in read_gpkg_table(db, 'roads', ('road_class',)):
        roads['main' if rc == 'main' else 'other'] += [p[:, :2] for p in geom[1]]
    # lattice cells that are mostly water: left to the basemap
    m = M[c]['meta']; X, Y = g.centers()
    col = np.rint((X - m['a']) / m['b']).astype(int); row = np.rint((Y - m['c']) / m['d']).astype(int)
    ii, jj = m['r1'] - row, col - m['c0']
    ok = (ii >= 0) & (ii < m['nrows']) & (jj >= 0) & (jj < m['ncols'])
    k = ii[ok] * m['ncols'] + jj[ok]
    nw = np.bincount(k, water[ok], minlength=m['nrows'] * m['ncols'])
    nt = np.bincount(k, minlength=m['nrows'] * m['ncols'])
    cell_water = (nw > 0.5 * np.maximum(nt, 1)).reshape(m['nrows'], m['ncols'])
    _BM[key] = dict(water=drawn, ext=ext, roads=roads, cell_water=cell_water, src=src)
    return _BM[key]


def draw_basemap(ax, c, half, clip=False):
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap
    b = basemap(c, half, clip)
    ax.set_facecolor(LAND)
    ax.imshow(np.where(b['water'], 1.0, np.nan), extent=b['ext'], cmap=ListedColormap([WATER]), vmin=0, vmax=2,
              interpolation='nearest', zorder=0.1, rasterized=True)
    ax.add_collection(LineCollection(b['roads']['other'], colors=ROAD, linewidths=0.12, zorder=0.2, rasterized=True))
    ax.add_collection(LineCollection(b['roads']['main'], colors=MAINROAD, linewidths=0.3, zorder=0.25, rasterized=True))
    if c in SEA_LABEL:
        t, x, y, rot = SEA_LABEL[c]
        ax.text(x, y, t, color='#5aa9e6', fontsize=5.6, style='italic', rotation=rot, ha='center', va='center',
                zorder=4, clip_on=True)
    return b


def full_array(c, k):
    """Factor k over the full lattice rectangle when the maps file has it (local and GEE factors);
    otherwise from the table (PS-only factors, prediction area only)."""
    if k in M[c]:
        return np.asarray(M[c][k], float)
    return lattice_array(c, T, k)


def lattice_array(c, T, k):
    """2-D array of factor k on the lattice from the table (covers PS-only columns too)."""
    m = M[c]['meta']; a = np.full((m['nrows'], m['ncols']), np.nan)
    if k not in T[c]:
        return a                                  # factor not built for this city: all grey
    t = T[c]; i = (m['r1'] - t.row).to_numpy(); j = (t.col - m['c0']).to_numpy()
    a[i, j] = t[k].to_numpy(float)
    return a


class G:                      # minimal grid object for figstyle.scalebar
    def __init__(self, lat): self.latm = np.radians(lat)
    def km_per_deg(self): return 111.195 * np.cos(self.latm), 111.195


def box(c, half_km, clip=False):
    """Panel limits: centred on the lattice, or on the city outline when clipping to the city."""
    if clip:
        p = np.vstack(city(c)['outline']); lon0, lat0 = (p[:, 0].min() + p[:, 0].max()) / 2, (p[:, 1].min() + p[:, 1].max()) / 2
    else:
        m = M[c]['meta']; lon0 = (m['west'] + m['east']) / 2; lat0 = (m['south'] + m['north']) / 2
    kx = 111.195 * np.cos(np.radians(lat0))
    return (lon0 - half_km[0] / kx, lon0 + half_km[0] / kx, lat0 - half_km[1] / 111.195, lat0 + half_km[1] / 111.195), lat0


def common_half_km(clip=False):
    w = h = 0
    if clip:                          # tight around the larger city outline, so the cities fill the boxes
        for c in CITIES:
            p = np.vstack(city(c)['outline']); lat0 = p[:, 1].mean()
            w = max(w, np.ptp(p[:, 0]) * 111.195 * np.cos(np.radians(lat0)) / 2)
            h = max(h, np.ptp(p[:, 1]) * 111.195 / 2)
        return w + 0.7, h + 0.7
    for c in CITIES:
        m = M[c]['meta']; lat0 = (m['south'] + m['north']) / 2
        w = max(w, (m['east'] - m['west']) * 111.195 * np.cos(np.radians(lat0)) / 2)
        h = max(h, (m['north'] - m['south']) * 111.195 / 2)
    return w + 1.8, h + 1.8           # margin so the basemap shows around the larger city too


def draw_map(ax, c, arr, cmap, norm, half, tag=True, nticks=4, clip=False):
    """Map of one factor on the dark grey canvas: over the prediction area, or (clip=True) clipped
    to the city corporation outline, with ward lines and a smooth outline instead of thana lines."""
    if clip:
        return draw_map_city(ax, c, arr, cmap, norm, half, tag, nticks)
    m = M[c]['meta']
    lims, lat0 = box(c, half)
    b = draw_basemap(ax, c, half)
    land = ~b['cell_water'] & M[c]['aoi']
    a = np.where(land, arr, np.nan)
    ext = [m['west'], m['east'], m['south'], m['north']]
    nodata = np.where(land & ~np.isfinite(arr), 1.0, np.nan)            # prediction area, no source data
    ax.imshow(nodata, extent=ext, cmap=plt.matplotlib.colors.ListedColormap([NODATA]), vmin=0, vmax=2,
              interpolation='nearest', zorder=0.9, rasterized=True)
    im = ax.imshow(a, extent=ext, cmap=cmap, norm=norm, interpolation='nearest', zorder=1, rasterized=True)
    for nm, rings in fs.thanas(c):
        for r in rings:
            ax.plot(r[:, 0], r[:, 1], color=BOUND, lw=0.22, zorder=3)
    ax.set_xlim(lims[0], lims[1]); ax.set_ylim(lims[2], lims[3])
    ax.set_aspect(1 / np.cos(np.radians(lat0)))
    fs.geo_ticks(ax, lims, n=nticks)
    ax.tick_params(labelsize=5)
    fs.north_arrow(ax, loc=(0.9, 0.84), size=0.07, fs=5.5, dark=True)
    fs.scalebar(ax, G(lat0), length_km=5, loc=(0.06, 0.04), fs=5.2, dark=True)
    if tag:
        fs.city_tag(ax, c, fs=6.2, dark=True)
    return im


def draw_map_city(ax, c, arr, cmap, norm, half, tag=True, nticks=4):
    from matplotlib.patches import PathPatch
    from matplotlib.colors import ListedColormap
    m = M[c]['meta']
    lims, lat0 = box(c, half, True)
    b = draw_basemap(ax, c, half, True)
    for nm, rings in fs.thanas(c):                      # thana lines on the canvas only, under the city
        for r in rings:
            ax.plot(r[:, 0], r[:, 1], color=BOUND, lw=0.18, alpha=0.6, zorder=0.3)
    cp = city_path(c)
    # inside the city: light grey (no source data) -> water -> factor, all clipped to the outline
    ax.add_patch(PathPatch(cp, fc=NODATA, ec='none', zorder=0.85))
    wat = ax.imshow(np.where(b['water'], 1.0, np.nan), extent=b['ext'], cmap=ListedColormap([WATER]), vmin=0,
                    vmax=2, interpolation='nearest', zorder=0.95, rasterized=True)
    ext = [m['west'], m['east'], m['south'], m['north']]
    im = ax.imshow(np.where(~b['cell_water'], arr, np.nan), extent=ext, cmap=cmap, norm=norm,
                   interpolation='nearest', zorder=1, rasterized=True)
    for art in (wat, im):
        art.set_clip_path(PathPatch(cp, transform=ax.transData))
    for poly in city(c)['wards']:
        for r in poly:
            ax.plot(r[:, 0], r[:, 1], color='0.35', lw=0.2, zorder=3)
    for r in city(c)['outline']:
        ax.plot(r[:, 0], r[:, 1], color='white', lw=0.6, zorder=3.5)
    ax.set_xlim(lims[0], lims[1]); ax.set_ylim(lims[2], lims[3])
    ax.set_aspect(1 / np.cos(np.radians(lat0)))
    fs.geo_ticks(ax, lims, n=nticks)
    ax.tick_params(labelsize=5)
    fs.north_arrow(ax, loc=(0.88, 0.86), size=0.06, fs=5.5, dark=True)
    fs.scalebar(ax, G(lat0), length_km=5, loc=(0.06, 0.035), fs=5.2, dark=True)
    if tag:
        fs.city_tag(ax, c, fs=6.2, dark=True)
    return im


def norm_for(k, vals):
    lo, hi = np.nanpercentile(vals, [2, 98])
    if k in DIVERGING:
        v = max(abs(lo), abs(hi)) or 1
        return TwoSlopeNorm(0, -v, v)
    return Normalize(lo, hi)


def map_sheet(keys, name, title=None, clip=False):
    """Rows of factor pairs: [Dhaka | Chattogram] x 2 per row, one colour bar per pair.
    clip=True: clip to the city corporation outlines (analysis/City_Maps), boxes fitted to the cities."""
    half = common_half_km(clip)
    ncol_pairs = 2
    nrow = int(np.ceil(len(keys) / ncol_pairs))
    pw = 39.5 * MM; ph = pw * half[1] / half[0]
    gx, gpair, top, cbh, rowgap = 5.5 * MM, 8 * MM, 6 * MM, 2.2 * MM, 15 * MM
    W = 2 * ncol_pairs * pw + ncol_pairs * gx + (ncol_pairs - 1) * gpair + 10 * MM
    H = top + nrow * (ph + rowgap) + 2 * MM
    fig = plt.figure(figsize=(W, H))
    letters = 'abcdefghijklmnopqrstuvwxyz'
    for n, k in enumerate(keys):
        r, q = divmod(n, ncol_pairs)
        x0 = 9 * MM + q * (2 * pw + gx + gpair)
        y0 = H - top - (r + 1) * ph - r * rowgap
        arrs = {c: full_array(c, k) * (1e-3 if k in KM else 1) for c in CITIES}
        # colour scale from the mapped area (prediction area, or the city when clipping)
        area = {c: city_cells(c) if clip else M[c]['aoi'] for c in CITIES}
        pooled = np.concatenate([arrs[c][area[c] & np.isfinite(arrs[c])] for c in CITIES])
        grp = GROUP.get(k, 'Terrain')
        cm_ = 'RdBu_r' if (k in DIVERGING and grp != 'Geology') else CMAP[grp]
        cmap = plt.get_cmap(cm_) if isinstance(cm_, str) else cm_
        nm = norm_for(k, pooled)
        for i, c in enumerate(CITIES):
            ax = fig.add_axes([(x0 + i * (pw + gx)) / W, y0 / H, pw / W, ph / H])
            im = draw_map(ax, c, arrs[c], cmap, nm, half, tag=True, clip=clip)
            if not np.isfinite(arrs[c][area[c]]).any():           # city-specific factor of the other city
                ax.text(0.5, 0.5, f'not defined for\n{fs.CITY_NAME[c]}' if hasattr(fs, 'CITY_NAME') else
                        'not defined here', transform=ax.transAxes, ha='center', va='center', fontsize=6,
                        color='0.25', zorder=6, bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
            if i == 0:
                ax.text(-0.02, 1.03, letters[n], transform=ax.transAxes, fontsize=8, fontweight='bold',
                        ha='right', va='bottom')
        cax = fig.add_axes([(x0 + 0.15 * pw) / W, (y0 - 8.2 * MM) / H, (2 * pw + gx - 0.3 * pw) / W, cbh / H])
        cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both', extendfrac=0.03)
        cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5.5, length=2, pad=1)
        cb.set_label(f'{LABEL.get(k, k)}  [{grp}]', fontsize=6.2, labelpad=1.5)
    save(fig, name)


def display_order(keys):
    """Keys grouped by theme (groups in order of first appearance in PRIORITY), PRIORITY order within."""
    groups = list(dict.fromkeys(GROUP[k] for k in PRIORITY))
    return sorted(keys, key=lambda k: (groups.index(GROUP[k]) if GROUP[k] in groups else 99, PRIORITY.index(k)
                                       if k in PRIORITY else 999))


def corr_fig():
    """Spearman correlation of the candidate factors at PS cells: Dhaka below the diagonal,
    Chattogram above it, one shared matrix so the labels stay legible with many factors.
    Grey labels: dropped in screening in both cities."""
    keep = pd.read_csv(O + 'kept_factors.csv', index_col=0)
    KEYS = display_order([k for k in PRIORITY if k in keep.index and all(k in T[c] for c in CITIES)])
    R = {c: T[c][T[c].has_ps.astype(bool)][KEYS].rank().corr().to_numpy() for c in CITIES}
    n = len(KEYS)
    Mx = np.where(np.tril(np.ones((n, n), bool), -1), R['dha'], np.where(np.triu(np.ones((n, n), bool), 1), R['ctg'], np.nan))
    side = min(170, 2.4 * n + 40)
    fig, ax = plt.subplots(figsize=(180 * MM, (side + 10) * MM))
    im = ax.imshow(Mx, cmap='RdBu_r', vmin=-1, vmax=1, interpolation='nearest')
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    lab = [LABEL.get(k, k).split(' (')[0] for k in KEYS]
    fsz = 4.6 if n <= 45 else 4.0
    ax.set_xticklabels(lab, rotation=90, fontsize=fsz); ax.set_yticklabels(lab, fontsize=fsz)
    for tl in (ax.get_xticklabels(), ax.get_yticklabels()):
        for t, k in zip(tl, KEYS):
            if all(keep.loc[k, c] != 'kept' for c in CITIES): t.set_color('0.6')
    g = [GROUP[k] for k in KEYS]
    for i in range(1, n):
        if g[i] != g[i - 1]:
            ax.axhline(i - .5, color='k', lw=0.4); ax.axvline(i - .5, color='k', lw=0.4)
    ax.plot([-.5, n - .5], [-.5, n - .5], color='0.3', lw=0.5)
    ax.tick_params(length=0, pad=1)
    nd, nc = int(T['dha'].has_ps.astype(bool).sum()), int(T['ctg'].has_ps.astype(bool).sum())
    ax.text(0.02, 0.02, f"Dhaka (PS cells, n = {nd:,})", transform=ax.transAxes, fontsize=6.5, ha='left', va='bottom')
    ax.text(0.98, 0.98, f"Chattogram (n = {nc:,})", transform=ax.transAxes, fontsize=6.5, ha='right', va='top')
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("Spearman's ρ", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    fig.subplots_adjust(left=0.2, right=0.95, bottom=0.17, top=0.99)
    save(fig, 'FigS_factor_correlation')


def assoc_fig():
    """Univariate Spearman rho with vU per factor and city: a first look, not a model."""
    keep = pd.read_csv(O + 'kept_factors.csv', index_col=0)
    rows = []
    for c in CITIES:
        s = pd.read_csv(O + f'screen_{c}.csv', index_col=0)
        rows.append(s.rho_vU.rename(c))
    R = pd.concat(rows, axis=1)
    R = R.loc[display_order([k for k in PRIORITY if k in R.index])]
    fig, ax = plt.subplots(figsize=(88 * MM, max(128, 3.3 * len(R) + 22) * MM))
    y = np.arange(len(R))[::-1]
    for c, off in zip(CITIES, (0.17, -0.17)):
        kept = (keep[c] == 'kept').reindex(R.index).fillna(False).to_numpy(bool)
        col = fs.CITY_COL[c]
        ax.hlines(y + off, 0, R[c], color=col, lw=0.8)
        ax.scatter(R[c][kept], (y + off)[kept], s=12, color=col, zorder=3, label=fs.CITY_NAME[c], ec='white', lw=0.4)
        ax.scatter(R[c][~kept], (y + off)[~kept], s=12, facecolor='white', ec=col, lw=0.7, zorder=3)
    ax.axvline(0, color='0.3', lw=0.5)
    ax.set_yticks(y); ax.set_yticklabels([LABEL.get(k, k).split(' (')[0] for k in R.index], fontsize=5.6)
    g = [GROUP[k] for k in R.index]
    for i in range(1, len(g)):
        if g[i] != g[i - 1]:
            ax.axhline(y[i] + .5, color='0.75', lw=0.4)
    for i, grp in enumerate(g):
        if i == 0 or g[i] != g[i - 1]:
            ax.text(0.99, y[i] + 0.35, grp, transform=ax.get_yaxis_transform(), ha='right', va='top',
                    fontsize=5.5, color='0.35', style='italic')
    ax.set_xlabel("Spearman's ρ with vertical rate $v_U$ (PS cells)")
    ax.set_xlim(-0.6, 0.6); ax.set_ylim(-0.8, len(R) - 0.2)
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 0.0), fontsize=6, handletextpad=0.3)
    ax.text(0.98, 0.015, 'open marker = dropped in screening', transform=ax.transAxes, fontsize=5.3, color='0.35', ha='right')
    fs.strip_spines(ax)
    ax.annotate('negative ρ: higher factor value goes with faster subsidence', xy=(0.5, 0), xycoords='axes fraction',
                xytext=(0, -30), textcoords='offset points', fontsize=5.2, color='0.35', ha='center', va='top')
    save(fig, 'FigS_factor_association')


def dem_fig():
    # PS cells: the cells the models are trained on (agreement over the whole area is higher)
    t = T['ctg']; t = t[t.has_ps.astype(bool)]
    try:
        own = json.load(open(D + 'ctg_sources.json'))['terrain']['ground']
    except Exception:
        own = 'this study'
    own = 'SRTM ground proxy' if own.startswith('SRTM') else own.split(' (')[0]
    pairs = [('elev', 'glo_elev', 'Elevation (m)'), ('slope', 'glo_slope', 'Slope (°)'),
             ('twi', 'glo_twi', 'TWI'), ('curv_prof', 'glo_curv_prof', 'Profile curvature')]
    fig, axs = plt.subplots(1, 4, figsize=(180 * MM, 46 * MM))
    for ax, (a, b, lab) in zip(axs, pairs):
        m = t[[a, b]].dropna(); x, yv = m[a], m[b]
        lo, hi = np.nanpercentile(np.r_[x, yv], [1, 99])
        hb = ax.hexbin(x, yv, gridsize=45, extent=(lo, hi, lo, hi), mincnt=1, cmap='Greys', bins='log', lw=0.1)
        ax.plot([lo, hi], [lo, hi], color=fs.CITY_COL['ctg'], lw=0.7)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
        r = spearmanr(x, yv)[0]
        ax.text(0.04, 0.95, f'{lab}\nρ = {r:.2f}', transform=ax.transAxes, va='top', fontsize=6)
        ax.set_xlabel(f'{own} (this study)', fontsize=5.8)
        if ax is axs[0]: ax.set_ylabel('Copernicus GLO-30 (CTG kit)', fontsize=5.8)
        ax.tick_params(labelsize=5.5)
    fig.subplots_adjust(wspace=0.35, left=0.07, right=0.99, bottom=0.22, top=0.97)
    save(fig, 'FigS_dem_check')


