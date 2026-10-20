"""The 100 m analysis lattice and the tools that put fine rasters onto it.

The 100 m cells of the three-city decomposition sit on a regular lon/lat lattice:
    lon = a + b*col,   lat = c + d*row      (row increases northwards)
Every factor is aggregated exactly onto that lattice by assigning each fine pixel
(30 m terrain/geology, 10 m land cover/OSM) to the cell that contains its centre.
"""
import json, time
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from PIL import Image
import tifffile
from geo import ll2utm, utm2ll, Grid, read_gpkg_table
from paths import INSAR, FACTORS, PARTS, CITY_MAPS

Image.MAX_IMAGE_PIXELS = None
CITY = {'dha': dict(name='Dhaka', lc='Dhaka', dem='n23_e090', boundary='Dhaka_city_bound.gpkg'),
        'ctg': dict(name='Chattogram', lc='Chittagong', dem='n22_e091', boundary='Chattogram_city_boundary.gpkg')}
AOI_PAD_CELLS = 10          # prediction area: cells within 1 km of a PS cell
BUF = 5000.                 # work-grid buffer beyond the prediction area (m)
DIST_CAP = 2000.            # all distance factors are capped here (m); beyond it a distance acts as a
                            # regional gradient (a coordinate) rather than a local effect
MIN_VALID = 0.5             # a cell needs this share of valid fine pixels, else NaN


def city_rings(c):
    """City corporation outline (lon/lat rings) from analysis/City_Maps; [] if the file is missing."""
    import sqlite3
    from geo import gpkg_geom
    fn = CITY_MAPS / CITY[c]['boundary']
    if not fn.exists():
        return []
    db = sqlite3.connect(str(fn))
    for t, in db.execute("select table_name from gpkg_contents where data_type = 'features'"):
        col = db.execute("select column_name from gpkg_geometry_columns where table_name = ?", (t,)).fetchone()[0]
        rows = [r[0] for r in db.execute(f'select "{col}" from "{t}"') if r[0]]
        if rows:
            return [p[:, :2] for b in rows for p in gpkg_geom(b)[1]]
    return []


class Lattice:
    def __init__(self, cells, rings=()):
        A = np.c_[np.ones(len(cells)), cells.col]
        self.a, self.b = np.linalg.lstsq(A, cells.lon, rcond=None)[0]
        B = np.c_[np.ones(len(cells)), cells.row]
        self.c, self.d = np.linalg.lstsq(B, cells.lat, rcond=None)[0]
        self.r0, self.r1 = int(cells.row.min()) - AOI_PAD_CELLS, int(cells.row.max()) + AOI_PAD_CELLS
        self.c0, self.c1 = int(cells.col.min()) - AOI_PAD_CELLS, int(cells.col.max()) + AOI_PAD_CELLS
        if len(rings):                       # also cover the whole city outline (+2 cells)
            p = np.vstack(rings)
            r, c = self.cell_of(p[:, 0], p[:, 1])
            self.r0, self.r1 = min(self.r0, int(r.min()) - 2), max(self.r1, int(r.max()) + 2)
            self.c0, self.c1 = min(self.c0, int(c.min()) - 2), max(self.c1, int(c.max()) + 2)
        self.nr, self.nc = self.r1 - self.r0 + 1, self.c1 - self.c0 + 1

    def lonlat(self, row, col):
        return self.a + self.b * np.asarray(col), self.c + self.d * np.asarray(row)

    def cell_of(self, lon, lat):
        col = np.rint((np.asarray(lon) - self.a) / self.b).astype(np.int64)
        row = np.rint((np.asarray(lat) - self.c) / self.d).astype(np.int64)
        return row, col

    def flat(self, row, col):
        """Index into the 2-D array (north-up: array row 0 = lattice row r1); -1 outside."""
        i = self.r1 - row; j = col - self.c0
        ok = (i >= 0) & (i < self.nr) & (j >= 0) & (j < self.nc)
        return np.where(ok, i * self.nc + j, -1)

    def meta(self):
        return dict(a=self.a, b=self.b, c=self.c, d=self.d, r0=self.r0, r1=self.r1, c0=self.c0, c1=self.c1,
                    nrows=self.nr, ncols=self.nc,
                    west=self.a + self.b * (self.c0 - .5), east=self.a + self.b * (self.c1 + .5),
                    south=self.c + self.d * (self.r0 - .5), north=self.c + self.d * (self.r1 + .5),
                    crs='EPSG:4326',
                    crsTransform=[self.b, 0, self.a + self.b * (self.c0 - .5), 0, -self.d, self.c + self.d * (self.r1 + .5)],
                    note='cell centre lon = a + b*col, lat = c + d*row; array row 0 = lattice row r1')


class Agg:
    """Aggregate a fine grid onto the lattice by pixel-centre assignment."""
    def __init__(self, lat_, grid):
        X, Y = grid.centers()
        lon, lat = utm2ll(X, Y)
        r, c = lat_.cell_of(lon, lat)
        self.k = lat_.flat(r, c).ravel()
        self.ok = self.k >= 0
        self.n = lat_.nr * lat_.nc
        self.shape = (lat_.nr, lat_.nc)
        self.cnt = np.bincount(self.k[self.ok], minlength=self.n).astype(float)

    def mean(self, a, min_valid=0.0):
        """Cell mean of the finite fine pixels; NaN where fewer than min_valid of the cell's
        pixels are finite (min_valid = 0: any finite pixel is enough)."""
        a = np.asarray(a, float).ravel()
        m = self.ok & np.isfinite(a)
        s = np.bincount(self.k[m], a[m], minlength=self.n)
        c = np.bincount(self.k[m], minlength=self.n)
        good = (c > 0) & (c >= min_valid * self.cnt)
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(good, s / np.maximum(c, 1), np.nan).reshape(self.shape)

    def majority(self, a, ncls):
        a = np.asarray(a).ravel().astype(np.int64)
        m = self.ok & (a >= 0)
        h = np.bincount(self.k[m] * ncls + a[m], minlength=self.n * ncls).reshape(self.n, ncls)
        out = h.argmax(1).astype(float); out[h.sum(1) == 0] = np.nan
        return out.reshape(self.shape)


class Ctx:
    """Everything a factor script needs for one city: lattice, prediction area, work grids."""
    def __init__(self, c):
        self.c, self.cfg, self.t0 = c, CITY[c], time.time()
        self.cells = pd.read_csv(INSAR / 'target' / f'{c}_cells_plus.csv')
        self.rings = city_rings(c)
        L = self.L = Lattice(self.cells, self.rings)
        self.R, self.C = np.meshgrid(np.arange(L.r1, L.r0 - 1, -1), np.arange(L.c0, L.c1 + 1), indexing='ij')
        self.LON, self.LAT = L.lonlat(self.R, self.C)
        hp = np.zeros(L.nr * L.nc, bool); hp[L.flat(self.cells.row.values, self.cells.col.values)] = True
        self.has_ps = hp.reshape(L.nr, L.nc)
        # prediction area: within 1 km of a PS cell, plus every cell of the city corporation area
        # (e.g. northern Chattogram, in a Sentinel-1 frame that was not processed: predicted, never trained on)
        self.near_ps = ndi.distance_transform_edt(~self.has_ps) <= AOI_PAD_CELLS
        self.in_city = np.zeros_like(self.has_ps)
        if self.rings:
            from matplotlib.path import Path
            ok = np.zeros(self.LON.size, bool)
            for r in self.rings:
                ok ^= Path(r).contains_points(np.c_[self.LON.ravel(), self.LAT.ravel()])
            self.in_city = ok.reshape(self.LON.shape)
        self.aoi = self.near_ps | self.in_city
        self.E, self.N = ll2utm(self.LON, self.LAT)
        E, N, aoi = self.E, self.N, self.aoi
        self.x0, self.y1 = E[aoi].min() - BUF, N[aoi].max() + BUF
        self.nx30 = int((E[aoi].max() - E[aoi].min() + 2 * BUF) / 30) + 1
        self.ny30 = int((N[aoi].max() - N[aoi].min() + 2 * BUF) / 30) + 1
        self.g30 = Grid(self.x0, self.y1, 30., self.nx30, self.ny30)
        self._A30 = self._A10 = None

    @property
    def A30(self):
        if self._A30 is None:
            self._A30 = Agg(self.L, self.g30)
        return self._A30

    @property
    def g10(self):   # 10 m grid: prediction area + 3 km
        return Grid(self.x0 + BUF - 3000, self.y1 - BUF + 3000, 10.,
                    int((self.nx30 * 30 - 2 * BUF + 6000) / 10), int((self.ny30 * 30 - 2 * BUF + 6000) / 10))

    @property
    def A10(self):
        if self._A10 is None:
            self._A10 = Agg(self.L, self.g10)
        return self._A10

    def log(self, *s):
        print(self.c, f'{time.time() - self.t0:5.0f}s', *s, flush=True)


def dist_m(mask, res, cap=DIST_CAP):
    if not mask.any():
        return np.full(mask.shape, cap)
    return np.minimum(ndi.distance_transform_edt(~mask) * res, cap)


def dist_censored(feature, mapped, res, cap=DIST_CAP):
    """Distance to `feature` for a source that covers only `mapped` pixels.

    A distance is kept only where it is known: the pixel is mapped and the nearest feature is
    closer than the nearest unmapped pixel (otherwise a closer feature could lie in the gap).
    If neither a feature nor a gap lies within `cap`, the value is `cap`. Everything else is NaN,
    so unmapped ground is never read as 'far from the feature'."""
    d = ndi.distance_transform_edt(~feature) * res if feature.any() else np.full(feature.shape, np.inf)
    gap = ndi.distance_transform_edt(mapped) * res if (~mapped).any() else np.full(feature.shape, np.inf)
    ok = mapped & ((d <= gap) | (gap >= cap))
    return np.where(ok, np.minimum(d, cap), np.nan)


def coverage_share(A, cls, cover, r_cells):
    """Share of class pixels among the mapped pixels within a disk of r_cells lattice cells.
    NaN where less than MIN_VALID of the disk is mapped."""
    import hydro
    num = hydro.focal_mean(A.mean((cls & cover).astype(float)), r_cells)
    den = hydro.focal_mean(A.mean(cover.astype(float)), r_cells)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den >= MIN_VALID, np.clip(num / den, 0, 1), np.nan)


def read_tif(fn):
    """Single-band GeoTIFF (LZW allowed) -> array, x/y of the upper-left corner, pixel size."""
    a = np.array(Image.open(fn)).astype(np.float32)
    with tifffile.TiffFile(fn) as t:
        tg = t.pages[0].tags
        tie = tg['ModelTiepointTag'].value; sc = tg['ModelPixelScaleTag'].value
        nd = tg['GDAL_NODATA'].value if 'GDAL_NODATA' in tg else None
    if nd is not None:
        a[a == float(nd)] = np.nan
    return a, tie[3], tie[4], sc[0], sc[1]


def polys_utm(rings_ll):
    out = []
    for r in rings_ll:
        r = np.asarray(r)[:, :2]; e, n = ll2utm(r[:, 0], r[:, 1]); out.append(np.c_[e, n])
    return out


def lines_utm(db, table, where=None, cols=()):
    out = []
    for r in read_gpkg_table(db, table, cols):
        kind, parts = r[0]
        if where is not None and not where(r[1:]):
            continue
        for p in parts:
            e, n = ll2utm(p[:, 0], p[:, 1]); out.append(np.c_[e, n])
    return out


def save_part(c, name, F, **extra):
    np.savez_compressed(PARTS / f'{c}_{name}.npz', **{k: np.asarray(v, np.float64) for k, v in F.items()},
                        _keys=json.dumps(list(F)), **{f'_x_{k}': json.dumps(v) for k, v in extra.items()})


def load_part(c, name):
    p = PARTS / f'{c}_{name}.npz'
    if not p.exists():
        return None, {}
    z = np.load(p, allow_pickle=False)
    keys = json.loads(str(z['_keys']))
    extra = {k[3:]: json.loads(str(z[k])) for k in z.files if k.startswith('_x_')}
    return {k: z[k] for k in keys}, extra
