"""Section 2, Fig. 1  Study areas, in the layout of Fig. 3 (Dhaka beside Chattogram, one map scale,
clipped to the city corporation boundaries on the dark grey canvas).

  (a) location of the two cities in Bangladesh (ADM1 divisions, Natural Earth coast and rivers)
  (b) surface geology, GSB 1:1M map as used in the three-city paper, grouped into nine units
  (c) vertical rate vU of the 100 m PS cells (2017-2023), with the GNSS stations
  (d) standard error of vU
Inputs    InSAR_data/basedata (BGD_ADM1, ne_coast, ne_rivers), factors_from_paper/geo_<c>.json,
          validation/gnss_station_fits.csv, analysis/factors (vU, vU_se per cell)
Output    manuscript/Figures/Fig01_study_area.png
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, TwoSlopeNorm, Normalize, BoundaryNorm
from matplotlib.patches import Polygon, PathPatch, Patch
from matplotlib.lines import Line2D
import figlib
import figstyle as fs
from figlib import MM, CITIES, LAND, WATER, BOUND
from geo import Grid, rasterize_polys_fast
from paths import INSAR

figlib.init()
M, T = figlib.M, figlib.T

# ------------------------------------------------------------------ geology units (grouped)
UNITS = [  # (label, colour, GSB codes)
    ('Alluvial silt and clay', '#9fd8c4', ['asc', 'asl']),
    ('Alluvial and deltaic sand, deltaic silt', '#f3e39c', ['asd', 'dsl', 'dsd']),
    ('Marsh clay and peat', '#3f9b8a', ['ppc']),
    ('Valley alluvium and colluvium', '#c5e6a6', ['ava']),
    ('Beach and dune sand', '#fbc98a', ['csd']),
    ('Madhupur Clay (Pleistocene)', '#c4683a', ['rm']),
    ('Dupi Tila and Dihing (Plio-Pleistocene)', '#e3a35a', ['QTdi', 'QTdt', 'QTdd']),
    ('Tipam and Boka Bil (Neogene)', '#8f5c35', ['Tt', 'Tbb']),
    ('Water in the 1990 map', '#6a9fc8', ['H2O', 'lake']),
]
CODE = {g: i + 1 for i, (_, _, gs) in enumerate(UNITS) for g in gs}
GEO_CMAP = ListedColormap([u[1] for u in UNITS])
GEO_NORM = BoundaryNorm(np.arange(0.5, len(UNITS) + 1), len(UNITS))


def geology_raster(c, half):
    """GSB polygons rasterised over the panel box (same ~0.0003 deg grid as the basemap)."""
    (w, e, s, n), _ = figlib.box(c, half, True)
    res = 0.0003
    g = Grid(w, n, res, int((e - w) / res) + 1, int((n - s) / res) + 1)
    feats = json.load(open(INSAR / 'factors_from_paper' / f'geo_{c}.json'))
    feats = sorted(feats, key=lambda f: -sum(len(r) for r in f['rings']))       # big first, small on top
    feats = [f for f in feats if f['glg'] in CODE]
    a = rasterize_polys_fast(g, [[np.asarray(r)[:, :2] for r in f['rings']] for f in feats],
                             [CODE[f['glg']] for f in feats], fill=0).astype(float)
    a[a == 0] = np.nan
    return a, [w, w + g.nx * res, n - g.ny * res, n]


def draw_geology(ax, c, half):
    figlib.draw_map(ax, c, np.full((M[c]['meta']['nrows'], M[c]['meta']['ncols']), np.nan),
                    'Greys', Normalize(0, 1), half, tag=True, clip=True)
    a, ext = geology_raster(c, half)
    xl, yl, asp = ax.get_xlim(), ax.get_ylim(), ax.get_aspect()
    im = ax.imshow(a, extent=ext, cmap=GEO_CMAP, norm=GEO_NORM, interpolation='nearest', zorder=0.9,
                   rasterized=True)
    im.set_clip_path(PathPatch(figlib.city_path(c), transform=ax.transData))
    ax.set_xlim(xl); ax.set_ylim(yl); ax.set_aspect(asp)   # imshow resets limits and aspect
    return sorted(np.unique(a[np.isfinite(a)]).astype(int))


# ------------------------------------------------------------------ location map
def coast_sea():
    """Sea polygon: the mainland coastline (Pakistan to Viet Nam) closed far to the south."""
    f = fs.load_geojson('ne_coast.geojson')
    a = np.asarray(f[1]['geometry']['coordinates'][0]); b = np.asarray(f[0]['geometry']['coordinates'][0])
    chain = np.vstack([a, b[1:]])
    sea = np.vstack([chain, [[chain[-1, 0], -20], [chain[0, 0], -20]]])
    islands = [np.asarray(l) for k, ft in enumerate(f) if k > 1 for l in ft['geometry']['coordinates']]
    return sea, islands, chain


def draw_location(ax, W_mm, H_mm):
    lat0, lon0, dlat = 23.65, 90.35, 6.6
    dlon = dlat * (W_mm / H_mm) / np.cos(np.radians(lat0))
    lims = (lon0 - dlon / 2, lon0 + dlon / 2, lat0 - dlat / 2, lat0 + dlat / 2)
    ax.set_facecolor('#2e2e2e')                                           # neighbouring countries
    sea, islands, chain = coast_sea()
    ax.add_patch(Polygon(sea, fc=WATER, ec='none', zorder=0.5))
    for r in islands:
        ax.add_patch(Polygon(r, fc='#2e2e2e', ec='none', zorder=0.6))
    for f in fs.load_geojson('BGD_ADM1.geojson'):
        for r in fs.rings_of(f['geometry']):
            ax.add_patch(Polygon(r, fc='#4f4f4f', ec='none', zorder=1))
    for f in fs.load_geojson('BGD_ADM1.geojson'):
        for r in fs.rings_of(f['geometry']):
            ax.plot(r[:, 0], r[:, 1], color=BOUND, lw=0.3, zorder=2)
    for f in fs.load_geojson('ne_rivers.geojson'):
        for r in fs.rings_of(f['geometry']):
            ax.plot(r[:, 0], r[:, 1], color='#3f6f96', lw=0.5, zorder=1.5)
    ax.plot(chain[:, 0], chain[:, 1], color='#5c5c5c', lw=0.4, zorder=1.6)
    halo = [pe.withStroke(linewidth=2.2, foreground='#1a1a1a')]
    for c, (dx, dy, ha) in {'dha': (-0.35, 0.25, 'right'), 'ctg': (0.45, -0.1, 'left')}.items():
        for r in figlib.city(c)['outline']:
            ax.add_patch(Polygon(r, fc=fs.CITY_COL[c], ec='white', lw=0.5, zorder=4))
        p = np.vstack(figlib.city(c)['outline']); x, y = p[:, 0].mean(), p[:, 1].mean()
        ax.annotate(fs.CITY_NAME[c], (x, y), (x + dx, y + dy), ha=ha, va='center', fontsize=6.5,
                    fontweight='bold', color='white', path_effects=halo, zorder=6,
                    arrowprops=dict(arrowstyle='-', color='white', lw=0.5))
    ax.text(89.6, 21.0, 'Bay of Bengal', color='#5aa9e6', fontsize=6, style='italic', ha='center', zorder=5)
    ax.text(90.15, 24.45, 'BANGLADESH', color='0.75', fontsize=5.5, ha='center', zorder=5, path_effects=halo)
    ax.text(87.0, 25.6, 'INDIA', color='0.6', fontsize=5.2, ha='center', zorder=5)
    ax.text(93.9, 22.6, 'MYANMAR', color='0.6', fontsize=5.2, ha='center', zorder=5, rotation=70)
    ax.set_xlim(lims[0], lims[1]); ax.set_ylim(lims[2], lims[3])
    ax.set_aspect(1 / np.cos(np.radians(lat0)))
    fs.geo_ticks(ax, lims, n=4, dec=0)
    ax.tick_params(labelsize=5)
    fs.north_arrow(ax, loc=(0.94, 0.86), size=0.06, fs=5.5, dark=True)
    fs.scalebar(ax, figlib.G(lat0), length_km=100, loc=(0.05, 0.035), fs=5.2, dark=True)


# ------------------------------------------------------------------ GNSS
def gnss(c):
    f = INSAR / 'validation' / 'gnss_station_fits.csv'
    if not f.exists():
        return pd.DataFrame(columns=['station', 'kind', 'lon', 'lat'])
    g = pd.read_csv(f)
    return g[g.city == fs.CITY_NAME[c]]


def draw_gnss(ax, c):
    halo = [pe.withStroke(linewidth=2, foreground='#1a1a1a')]
    for _, s in gnss(c).iterrows():
        mk = '^' if s.kind == 'continuous' else 's'
        ax.plot(s.lon, s.lat, mk, ms=4.2 if mk == '^' else 3.4, mfc='white', mec='black', mew=0.6, zorder=8)
        ax.text(s.lon + 0.006, s.lat + 0.004, s.station, fontsize=5, color='white', zorder=8,
                path_effects=halo, fontweight='bold')


# ------------------------------------------------------------------ layout (as in figlib.map_sheet)
half = figlib.common_half_km(True)
pw = 39.5 * MM; ph = pw * half[1] / half[0]
gx, gpair, top, cbh, rowgap = 5.5 * MM, 8 * MM, 6 * MM, 2.2 * MM, 20 * MM   # room for the geology legend
W = 4 * pw + 2 * gx + gpair + 10 * MM
H = top + 2 * (ph + rowgap) + 2 * MM
fig = plt.figure(figsize=(W, H))


def slot(n):
    r, q = divmod(n, 2)
    return 9 * MM + q * (2 * pw + gx + gpair), H - top - (r + 1) * ph - r * rowgap


def letter(ax, s):
    ax.text(-0.02, 1.03, s, transform=ax.transAxes, fontsize=8, fontweight='bold', ha='right', va='bottom')


def cbar_axes(x0, y0):
    return fig.add_axes([(x0 + 0.15 * pw) / W, (y0 - 8.2 * MM) / H, (2 * pw + gx - 0.3 * pw) / W, cbh / H])


# (a) location
x0, y0 = slot(0)
ax = fig.add_axes([x0 / W, y0 / H, (2 * pw + gx) / W, ph / H])
draw_location(ax, (2 * pw + gx) / MM, ph / MM)
letter(ax, 'a')

# (b) geology
x0, y0 = slot(1); used = set()
for i, c in enumerate(CITIES):
    ax = fig.add_axes([(x0 + i * (pw + gx)) / W, y0 / H, pw / W, ph / H])
    used |= set(draw_geology(ax, c, half))
    if i == 0:
        letter(ax, 'b')
lax = fig.add_axes([x0 / W, (y0 - 19 * MM) / H, (2 * pw + gx) / W, 10.5 * MM / H]); lax.axis('off')
hs = [Patch(fc=UNITS[k - 1][1], ec='0.3', lw=0.3, label=UNITS[k - 1][0]) for k in range(1, len(UNITS) + 1) if k in used]
lax.legend(handles=hs, loc='upper center', ncol=2, fontsize=5.4, frameon=False, handlelength=1.4,
           handleheight=0.9, columnspacing=1.2, labelspacing=0.35, borderaxespad=0)
lax.text(0.5, 1.02, 'Surface geology (GSB 1:1M)', transform=lax.transAxes, ha='center', va='bottom', fontsize=6.2)

# (c) vU and (d) its standard error at the PS cells
specs = [('vU', 'RdBu', 'Vertical rate $v_U$ (mm/yr), 2017–2023; negative = subsidence', True),
         ('vU_se', 'YlGnBu', 'Standard error of $v_U$ (mm/yr)', False)]
for n, (k, cm, lab, div) in enumerate(specs):
    x0, y0 = slot(2 + n)
    arrs = {c: figlib.lattice_array(c, T, k) for c in CITIES}
    pooled = np.concatenate([arrs[c][np.isfinite(arrs[c])] for c in CITIES])
    lo, hi = np.nanpercentile(pooled, [2, 98])
    nm = TwoSlopeNorm(0, -max(-lo, hi), max(-lo, hi)) if div else Normalize(0, hi)
    for i, c in enumerate(CITIES):
        ax = fig.add_axes([(x0 + i * (pw + gx)) / W, y0 / H, pw / W, ph / H])
        im = figlib.draw_map(ax, c, arrs[c], plt.get_cmap(cm), nm, half, tag=True, clip=True)
        if k == 'vU':
            draw_gnss(ax, c)
        if i == 0:
            letter(ax, 'cd'[n])
    cax = cbar_axes(x0, y0)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both' if div else 'max', extendfrac=0.03)
    cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5.5, length=2, pad=1)
    cb.set_label(lab, fontsize=6.2, labelpad=1.5)
    if k == 'vU':                                    # GNSS key beside the colour bar
        cax.legend(handles=[Line2D([], [], marker='^', ls='', mfc='white', mec='k', mew=0.6, ms=4.2, label='GNSS, continuous'),
                            Line2D([], [], marker='s', ls='', mfc='white', mec='k', mew=0.6, ms=3.4, label='GNSS, campaign')],
                   loc='upper center', bbox_to_anchor=(0.5, -2.9), ncol=2, fontsize=5.4, frameon=False,
                   handletextpad=0.3, columnspacing=1.2, borderaxespad=0)
figlib.save(fig, 'Fig01_study_area')
