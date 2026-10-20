"""Small geospatial helpers (no GDAL available): UTM 46N, shapefile, GPKG WKB,
polygon/line rasterisation, SRTM tile reading."""
import struct, sqlite3, numpy as np, tifffile
from matplotlib.path import Path

# ---------------------------------------------------------------- UTM (Krueger series)
_a, _f = 6378137.0, 1 / 298.257223563
_k0, _lon0, _FE = 0.9996, np.radians(93.0), 500000.0  # UTM zone 46N
_n = _f / (2 - _f)
_A = _a / (1 + _n) * (1 + _n**2 / 4 + _n**4 / 64)
_al = [_n/2 - 2*_n**2/3 + 5*_n**3/16, 13*_n**2/48 - 3*_n**3/5, 61*_n**3/240]
_be = [_n/2 - 2*_n**2/3 + 37*_n**3/96, _n**2/48 + _n**3/15, 17*_n**3/480]
_de = [2*_n - 2*_n**2/3 - 2*_n**3, 7*_n**2/3 - 8*_n**3/5, 56*_n**3/15]


def ll2utm(lon, lat):
    lon, lat = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
    e = np.sqrt(_f * (2 - _f))
    t = np.sinh(np.arctanh(np.sin(lat)) - e * np.arctanh(e * np.sin(lat)))
    dl = lon - _lon0
    xi, eta = np.arctan2(t, np.cos(dl)), np.arctanh(np.sin(dl) / np.sqrt(1 + t**2))
    x, y = eta.copy(), xi.copy()
    for j, a in enumerate(_al, 1):
        x += a * np.cos(2*j*xi) * np.sinh(2*j*eta)
        y += a * np.sin(2*j*xi) * np.cosh(2*j*eta)
    return _FE + _k0 * _A * x, _k0 * _A * y


def utm2ll(E, N):
    xi = np.asarray(N, float) / (_k0 * _A)
    eta = (np.asarray(E, float) - _FE) / (_k0 * _A)
    xp, ep = xi.copy(), eta.copy()
    for j, b in enumerate(_be, 1):
        xp -= b * np.sin(2*j*xi) * np.cosh(2*j*eta)
        ep -= b * np.cos(2*j*xi) * np.sinh(2*j*eta)
    chi = np.arcsin(np.sin(xp) / np.cosh(ep))
    lat = chi.copy()
    for j, d in enumerate(_de, 1):
        lat += d * np.sin(2*j*chi)
    lon = _lon0 + np.arctan2(np.sinh(ep), np.cos(xp))
    return np.degrees(lon), np.degrees(lat)


# ---------------------------------------------------------------- shapefile
def read_dbf(fn):
    b = open(fn, 'rb').read()
    n = struct.unpack('<I', b[4:8])[0]; hl, rl = struct.unpack('<HH', b[8:12])
    flds = []
    for i in range(32, hl - 1, 32):
        nm = b[i:i+11].split(b'\0')[0].decode(); flds.append((nm, chr(b[i+11]), b[i+16]))
    rows = []
    for r in range(n):
        o = hl + r * rl + 1; rec = {}
        for nm, t, ln in flds:
            v = b[o:o+ln].decode('latin1').strip(); o += ln
            rec[nm] = float(v) if t in 'NF' and v else v
        rows.append(rec)
    return rows


def read_shp_polygons(fn):
    """Return list of lists of rings (each ring an (n,2) array)."""
    b = open(fn, 'rb').read(); o = 100; out = []
    while o < len(b):
        _, clen = struct.unpack('>II', b[o:o+8]); o += 8
        st = struct.unpack('<i', b[o:o+4])[0]
        if st in (5, 15, 25):
            npart, npt = struct.unpack('<ii', b[o+36:o+44])
            parts = np.frombuffer(b, '<i4', npart, o+44)
            pts = np.frombuffer(b, '<f8', 2*npt, o+44+4*npart).reshape(-1, 2)
            idx = list(parts) + [npt]
            out.append([pts[idx[i]:idx[i+1]] for i in range(npart)])
        else:
            out.append([])
        o += clen * 2
    return out


# ---------------------------------------------------------------- GPKG WKB
def _wkb(b, o=0):
    bo = '<' if b[o] == 1 else '>'; t = struct.unpack(bo+'I', b[o+1:o+5])[0] % 1000; o += 5
    if t == 1:
        return ('pt', [np.frombuffer(b, bo+'f8', 2, o).reshape(1, 2)]), o + 16
    if t == 2:
        n = struct.unpack(bo+'I', b[o:o+4])[0]
        return ('ln', [np.frombuffer(b, bo+'f8', 2*n, o+4).reshape(-1, 2)]), o + 4 + 16*n
    if t == 3:
        nr = struct.unpack(bo+'I', b[o:o+4])[0]; o += 4; rings = []
        for _ in range(nr):
            n = struct.unpack(bo+'I', b[o:o+4])[0]
            rings.append(np.frombuffer(b, bo+'f8', 2*n, o+4).reshape(-1, 2)); o += 4 + 16*n
        return ('pg', rings), o
    if t in (4, 5, 6, 7):
        n = struct.unpack(bo+'I', b[o:o+4])[0]; o += 4; kind, parts = None, []
        for _ in range(n):
            (k, p), o = _wkb(b, o); kind = k; parts += p
        return (kind, parts), o
    raise ValueError(t)


def gpkg_geom(blob):
    b = bytes(blob); flags = b[3]; env = (flags >> 1) & 7
    off = 8 + {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]
    return _wkb(b, off)[0]


def read_gpkg_table(fn, table, cols=()):
    c = sqlite3.connect(fn)
    q = 'select geom' + ''.join(f',"{k}"' for k in cols) + f' from "{table}"'
    out = []
    for r in c.execute(q):
        if r[0] is None: continue
        out.append((gpkg_geom(r[0]),) + tuple(r[1:]))
    return out


# ---------------------------------------------------------------- grids
class Grid:
    """North-up UTM grid: x0,y0 = upper-left corner, res metres."""
    def __init__(self, x0, y0, res, nx, ny):
        self.x0, self.y0, self.res, self.nx, self.ny = x0, y0, res, nx, ny

    def centers(self):
        xs = self.x0 + (np.arange(self.nx) + .5) * self.res
        ys = self.y0 - (np.arange(self.ny) + .5) * self.res
        return np.meshgrid(xs, ys)

    def ij(self, x, y):
        return ((self.y0 - np.asarray(y)) / self.res).astype(int), ((np.asarray(x) - self.x0) / self.res).astype(int)


def rasterize_polys(grid, polys, values, fill=0, dtype=np.int16):
    """polys: list of ring-lists in grid CRS; even-odd rule within each feature."""
    out = np.full((grid.ny, grid.nx), fill, dtype)
    X, Y = grid.centers()
    for rings, v in zip(polys, values):
        if not len(rings): continue
        allp = np.vstack(rings)
        xmn, ymn = allp.min(0); xmx, ymx = allp.max(0)
        j0, j1 = max(int((xmn - grid.x0) / grid.res), 0), min(int((xmx - grid.x0) / grid.res) + 1, grid.nx)
        i0, i1 = max(int((grid.y0 - ymx) / grid.res), 0), min(int((grid.y0 - ymn) / grid.res) + 1, grid.ny)
        if j1 <= j0 or i1 <= i0: continue
        pts = np.c_[X[i0:i1, j0:j1].ravel(), Y[i0:i1, j0:j1].ravel()]
        inside = np.zeros(len(pts), bool)
        for r in rings:
            inside ^= Path(r).contains_points(pts)
        sub = out[i0:i1, j0:j1].reshape(-1); sub[inside] = v
        out[i0:i1, j0:j1] = sub.reshape(i1 - i0, j1 - j0)
    return out


def rasterize_lines(grid, lines, step=None):
    """Burn polylines (grid CRS) into a boolean mask by dense sampling."""
    step = step or grid.res / 3
    m = np.zeros((grid.ny, grid.nx), bool)
    for ln in lines:
        if len(ln) < 2:
            seg = ln
        else:
            d = np.r_[0, np.cumsum(np.hypot(*np.diff(ln, axis=0).T))]
            t = np.arange(0, d[-1] + step, step)
            seg = np.c_[np.interp(t, d, ln[:, 0]), np.interp(t, d, ln[:, 1])]
        i, j = grid.ij(seg[:, 0], seg[:, 1])
        k = (i >= 0) & (i < grid.ny) & (j >= 0) & (j < grid.nx)
        m[i[k], j[k]] = True
    return m


def read_srtm(fn):
    """SRTM 1-arcsec tile -> array, lon0 (west edge), lat0 (north edge), d (deg)."""
    with tifffile.TiffFile(fn) as t:
        p = t.pages[0]; a = p.asarray().astype(np.float32)
        tie = p.tags['ModelTiepointTag'].value; sc = p.tags['ModelPixelScaleTag'].value
    a[a < -1000] = np.nan
    return a, tie[3], tie[4], sc[0], sc[1]


def sample_regular(arr, lon0, lat0, dx, dy, lon, lat):
    """Bilinear sample of a lon/lat regular raster (pixel-is-area corner origin)."""
    from scipy.ndimage import map_coordinates
    c = (np.asarray(lon) - lon0) / dx - .5; r = (lat0 - np.asarray(lat)) / dy - .5
    return map_coordinates(arr, [r, c], order=1, mode='nearest')


def rasterize_polys_fast(grid, polys, values, fill=0, dtype=np.int16):
    """Scanline rasterisation with Pillow; even-odd rule within each feature
    (outer ring XOR holes). Later features overwrite earlier ones."""
    from PIL import Image, ImageDraw
    out = np.full((grid.ny, grid.nx), fill, dtype)
    for rings, v in zip(polys, values):
        if not len(rings): continue
        allp = np.vstack(rings)
        xmn, ymn = allp.min(0); xmx, ymx = allp.max(0)
        j0, j1 = max(int((xmn - grid.x0) / grid.res) - 1, 0), min(int((xmx - grid.x0) / grid.res) + 2, grid.nx)
        i0, i1 = max(int((grid.y0 - ymx) / grid.res) - 1, 0), min(int((grid.y0 - ymn) / grid.res) + 2, grid.ny)
        if j1 <= j0 or i1 <= i0: continue
        m = np.zeros((i1 - i0, j1 - j0), bool)
        for r in rings:
            if len(r) < 3: continue
            # pixel coordinates relative to the bbox; PIL samples pixel centres at +0.5
            px = (r[:, 0] - grid.x0) / grid.res - j0 - 0.5
            py = (grid.y0 - r[:, 1]) / grid.res - i0 - 0.5
            im = Image.new('1', (j1 - j0, i1 - i0), 0)
            ImageDraw.Draw(im).polygon(list(zip(px.tolist(), py.tolist())), fill=1)
            m ^= np.array(im, bool)
        sub = out[i0:i1, j0:j1]; sub[m] = v
    return out
