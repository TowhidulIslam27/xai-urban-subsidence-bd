"""Shared style and map helpers for the paper figures."""
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, Polygon, Rectangle
from matplotlib.collections import LineCollection, PatchCollection
import matplotlib.patheffects as pe

try:
    import cmcrameri.cm as cmc
    CM_DIV = cmc.vik_r          # subsidence = warm (brown/red), uplift = blue
    CM_DIV_E = plt.get_cmap("PuOr_r")   # east-west: purple = west, orange = east
    CM_CYC = cmc.romaO          # cyclic, for seasonal phase
    CM_SEQ = cmc.batlow
    CM_SEQ_L = cmc.lajolla
except Exception:                # fallback without cmcrameri
    CM_DIV = plt.get_cmap("RdBu")
    CM_DIV_E = plt.get_cmap("PuOr")
    CM_CYC = plt.get_cmap("twilight_shifted")
    CM_SEQ = plt.get_cmap("viridis")
    CM_SEQ_L = plt.get_cmap("YlOrRd")

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "basedata")

CITY_ORDER = ["dha", "ctg", "khu"]
CITY_NAME = {"dha": "Dhaka", "ctg": "Chattogram", "khu": "Khulna"}
CITY_COL = {"dha": "#0072B2", "ctg": "#D55E00", "khu": "#009E73"}   # Okabe-Ito
MM = 1 / 25.4
W1, W15, W2 = 88 * MM, 140 * MM, 180 * MM   # single, 1.5, double column


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "FreeSans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.titleweight": "bold", "axes.titlepad": 3,
        "legend.frameon": False, "legend.handlelength": 1.4,
        "lines.linewidth": 0.9, "patch.linewidth": 0.5,
        "mathtext.default": "regular",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 600, "figure.dpi": 150,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "image.interpolation": "nearest",
    })


def save(fig, outdir, name, formats=("pdf", "png")):
    os.makedirs(outdir, exist_ok=True)
    for f in formats:
        fig.savefig(os.path.join(outdir, f"{name}.{f}"), dpi=600 if f != "pdf" else 300)
    plt.close(fig)
    print("  saved", name)


def panel_label(ax, s, x=0.0, y=1.0, dx=-2, dy=2, **kw):
    ax.annotate(s, (x, y), xycoords="axes fraction", xytext=(dx, dy),
                textcoords="offset points", ha="right" if dx < 0 else "left",
                va="bottom", fontsize=8.5, fontweight="bold", **kw)


def strip_spines(ax, keep=("left", "bottom")):
    for s in ("left", "right", "top", "bottom"):
        ax.spines[s].set_visible(s in keep)


# ------------------------------------------------------------------ data
def load_cells(derived, c):
    return pd.read_csv(os.path.join(derived, f"{c}_cells_plus.csv"))


def load_ts(derived, c):
    z = np.load(os.path.join(derived, f"{c}_ts.npz"))
    return dict(U=z["U"], E=z["E"], dates=z["dates"].astype("datetime64[D]"),
                cell_id=z["cell_id"], lon_c=z["lon_c"], lat_c=z["lat_c"])


def load_geojson(name):
    path = os.path.join(BASE, name)
    if not os.path.exists(path):
        return []
    return json.load(open(path))["features"]


def rings_of(geom):
    g = geom
    if g["type"] == "Polygon":
        return [np.asarray(r)[:, :2] for r in g["coordinates"]]
    if g["type"] == "MultiPolygon":
        return [np.asarray(r)[:, :2] for p in g["coordinates"] for r in p]
    if g["type"] in ("LineString",):
        return [np.asarray(g["coordinates"])[:, :2]]
    if g["type"] == "MultiLineString":
        return [np.asarray(r)[:, :2] for r in g["coordinates"]]
    return []


THANAS = None


def thanas(city):
    global THANAS
    if THANAS is None:
        THANAS = load_geojson("thana_3city.geojson")
    return [(f["properties"]["name"], rings_of(f["geometry"]))
            for f in THANAS if f["properties"]["city"] == city]


# ------------------------------------------------------------- rasters
class CityGrid:
    """Exact 100 m raster of the decomposition cells, in lon/lat."""

    def __init__(self, df, lon_c, lat_c):
        self.df = df
        self.r0, self.r1 = df.row.min(), df.row.max()
        self.c0, self.c1 = df.col.min(), df.col.max()
        self.lon_c, self.lat_c = lon_c, lat_c
        self.ri = (df.row - self.r0).to_numpy()
        self.ci = (df.col - self.c0).to_numpy()
        self.shape = (self.r1 - self.r0 + 1, self.c1 - self.c0 + 1)
        self.extent = (lon_c[0] + lon_c[1] * self.c0, lon_c[0] + lon_c[1] * (self.c1 + 1),
                       lat_c[0] + lat_c[1] * self.r0, lat_c[0] + lat_c[1] * (self.r1 + 1))
        self.latm = np.deg2rad(df.lat.mean())

    def raster(self, values, mask=None):
        a = np.full(self.shape, np.nan)
        v = np.asarray(values, float).copy()
        if mask is not None:
            v[~np.asarray(mask)] = np.nan
        a[self.ri, self.ci] = v
        return a

    def km_per_deg(self):
        return 111.195 * np.cos(self.latm), 111.195

    def xy_to_lonlat(self, x_m, y_m):
        lon = self.lon_c[0] + self.lon_c[1] * (np.asarray(x_m) / 100.0)
        lat = self.lat_c[0] + self.lat_c[1] * (np.asarray(y_m) / 100.0)
        return lon, lat


def map_limits(grid, pad_km=0.6):
    kx, ky = grid.km_per_deg()
    e = grid.extent
    return (e[0] - pad_km / kx, e[1] + pad_km / kx, e[2] - pad_km / ky, e[3] + pad_km / ky)


def map_size_km(grid, pad_km=0.6):
    kx, ky = grid.km_per_deg()
    l = map_limits(grid, pad_km)
    return (l[1] - l[0]) * kx, (l[3] - l[2]) * ky


def setup_map(ax, grid, lims=None, ticks=True, nticks=3, land=True, bounds=True,
              bound_lw=0.3, bound_col="0.35", city=None):
    lims = lims or getattr(ax, "map_lims", None) or map_limits(grid)
    ax.set_xlim(lims[0], lims[1])
    ax.set_ylim(lims[2], lims[3])
    ax.set_aspect(1 / np.cos(grid.latm))
    ax.set_facecolor("white")
    if city is not None and (land or bounds):
        for nm, rings in thanas(city):
            for r in rings:
                if land:
                    ax.add_patch(Polygon(r, closed=True, fc="#f1f1ef", ec="none", zorder=0))
                if bounds:
                    ax.plot(r[:, 0], r[:, 1], color=bound_col, lw=bound_lw, zorder=3,
                            solid_joinstyle="round")
    if ticks:
        geo_ticks(ax, lims, nticks)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.6)


def _nice_step(span, n):
    raw = span / n
    for s in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
        if s >= raw:
            return s
    return 1.0


def geo_ticks(ax, lims, n=3, dec=2):
    sx = _nice_step(lims[1] - lims[0], n)
    sy = _nice_step(lims[3] - lims[2], n)
    mx, my = 0.07 * (lims[1] - lims[0]), 0.07 * (lims[3] - lims[2])
    xt = np.arange(np.ceil(lims[0] / sx) * sx, lims[1], sx)
    yt = np.arange(np.ceil(lims[2] / sy) * sy, lims[3], sy)
    xt = xt[(xt > lims[0] + mx) & (xt < lims[1] - mx)]
    yt = yt[(yt > lims[2] + my) & (yt < lims[3] - my)]
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_xticklabels([f"{v:.{dec}f}°E" for v in xt])
    ax.set_yticklabels([f"{v:.{dec}f}°N" for v in yt], rotation=90, va="center")
    ax.tick_params(axis="both", labelsize=5.8, pad=1.5, length=2)


def scalebar(ax, grid, length_km=5, loc=(0.06, 0.05), lw=1.6, fs=6, dark=False):
    kx, _ = grid.km_per_deg()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xs = x0 + loc[0] * (x1 - x0)
    ys = y0 + loc[1] * (y1 - y0)
    dl = length_km / kx
    fg, bg = ("white", "#1a1a1a") if dark else ("k", "white")
    halo = [pe.withStroke(linewidth=2.2, foreground=bg)]
    # alternating bar (white/black on a dark map, black/white on a light one)
    nseg = 2
    for i in range(nseg):
        ax.add_patch(Rectangle((xs + i * dl / nseg, ys), dl / nseg, (y1 - y0) * 0.008,
                               fc=fg if i % 2 == 0 else bg, ec=fg, lw=0.4, zorder=20))
    ax.text(xs + dl / 2, ys + (y1 - y0) * 0.018, f"{length_km:g} km", ha="center", va="bottom",
            fontsize=fs, zorder=21, color=fg, path_effects=halo)


def north_arrow(ax, loc=(0.92, 0.88), size=0.06, fs=6.5, dark=False):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x = x0 + loc[0] * (x1 - x0)
    y = y0 + loc[1] * (y1 - y0)
    h = size * (y1 - y0)
    w = h * 0.35 * (x1 - x0) / (y1 - y0) * (ax.bbox.height / ax.bbox.width)
    fg, bg = ("white", "#1a1a1a") if dark else ("k", "white")
    ax.add_patch(Polygon([[x, y + h], [x - w, y], [x, y + h * 0.28]], fc=fg, ec=fg, lw=0.4, zorder=20))
    ax.add_patch(Polygon([[x, y + h], [x + w, y], [x, y + h * 0.28]], fc=bg, ec=fg, lw=0.4, zorder=20))
    ax.text(x, y + h * 1.08, "N", ha="center", va="bottom", fontsize=fs, fontweight="bold", zorder=21,
            color=fg, path_effects=[pe.withStroke(linewidth=2, foreground=bg)])


def city_tag(ax, c, loc=(0.03, 0.97), fs=7.5, extra=None, dark=False):
    txt = CITY_NAME[c] if extra is None else f"{CITY_NAME[c]}\n{extra}"
    fg, bg = ("white", "#1a1a1a") if dark else ("k", "white")
    ax.text(*loc, txt, transform=ax.transAxes, ha="left", va="top", fontsize=fs,
            fontweight="bold", color=fg, zorder=25, linespacing=1.15,
            path_effects=[pe.withStroke(linewidth=2.4, foreground=bg)])


def map_scale_height(grids, total_w, gap):
    """Height (in) needed when maps fill total_w at a common km scale."""
    sizes = [map_size_km(g) for g in grids]
    s = (total_w - gap * (len(grids) - 1)) / sum(x[0] for x in sizes)
    return s * max(x[1] for x in sizes)


def same_scale_axes(fig, grids, left, bottom, total_w, height_max, gap):
    """Place one map axes per grid at an identical km scale (in figure inches)."""
    sizes = [map_size_km(g) for g in grids]
    sx = (total_w - gap * (len(grids) - 1)) / sum(s[0] for s in sizes)   # in/km
    sy = height_max / max(s[1] for s in sizes)
    s = min(sx, sy)
    fw, fh = fig.get_size_inches()
    axes = []
    x = left
    for (wkm, hkm) in sizes:
        w, h = wkm * s, hkm * s
        ax = fig.add_axes([x / fw, (bottom + height_max - h) / fh, w / fw, h / fh])
        axes.append(ax)
        x += w + gap
    return axes, s


def monsoon_shading(ax, years=range(2017, 2024), color="#dbe9f6", alpha=1.0, months=(6, 10)):
    for y in years:
        a = np.datetime64(f"{y}-{months[0]:02d}-01")
        b = np.datetime64(f"{y}-{months[1]:02d}-01")
        ax.axvspan(a, b, color=color, alpha=alpha, lw=0, zorder=0)


def add_cbar(fig, rect, mappable, label, ticks=None, extend="both", orientation="horizontal"):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(mappable, cax=cax, orientation=orientation, extend=extend,
                      extendfrac=0.04, ticks=ticks)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, width=0.4, length=2, pad=1.5)
    cb.set_label(label, fontsize=6.5, labelpad=2)
    return cb
