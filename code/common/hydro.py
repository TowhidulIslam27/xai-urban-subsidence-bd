"""Terrain and flow-routing helpers on a north-up metric grid (row 0 = north).

No GDAL/WhiteboxTools here, so the routines are written out:
priority-flood depression filling with an epsilon gradient (Barnes et al. 2014),
D8 receivers, flow accumulation, TWI, streams, HAND and terrain derivatives.
"""
import heapq
import numpy as np
from scipy import ndimage as ndi

D8 = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]


def fill_epsilon(z, eps=1e-4):
    """Priority-flood fill; every interior cell drains to the edge afterwards."""
    ny, nx = z.shape
    zf = z.astype(np.float64).ravel().copy()
    done = np.zeros(ny * nx, bool)
    heap = []
    edge = np.zeros((ny, nx), bool)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
    for k in np.flatnonzero(edge.ravel()):
        heap.append((zf[k], int(k)))
        done[k] = True
    heapq.heapify(heap)
    offs = [di * nx + dj for di, dj in D8]
    push, pop = heapq.heappush, heapq.heappop
    zl = zf.tolist(); dl = done.tolist()
    while heap:
        h, k = pop(heap)
        i, j = divmod(k, nx)
        for (di, dj), o in zip(D8, offs):
            ii, jj = i + di, j + dj
            if ii < 0 or jj < 0 or ii >= ny or jj >= nx:
                continue
            n = k + o
            if dl[n]:
                continue
            dl[n] = True
            v = zl[n]
            if v <= h:
                v = h + eps
                zl[n] = v
            push(heap, (v, n))
    return np.asarray(zl).reshape(ny, nx)


def d8_receivers(zf, res):
    """Index of the steepest downslope neighbour (-1 at outlets)."""
    ny, nx = zf.shape
    best = np.zeros((ny, nx)); rec = np.full((ny, nx), -1, np.int64)
    idx = np.arange(ny * nx).reshape(ny, nx)
    pad = np.pad(zf, 1, constant_values=np.inf)
    pidx = np.pad(idx, 1, constant_values=-1)
    for di, dj in D8:
        dist = res * (np.sqrt(2) if di and dj else 1.0)
        zn = pad[1 + di:1 + di + ny, 1 + dj:1 + dj + nx]
        drop = (zf - zn) / dist
        m = drop > best
        best[m] = drop[m]
        rec[m] = pidx[1 + di:1 + di + ny, 1 + dj:1 + dj + nx][m]
    return rec, best


def accumulate(zf, rec):
    """Upslope cell count (including the cell) by descending-elevation sweep."""
    order = np.argsort(-zf.ravel(), kind='stable')
    acc = np.ones(zf.size)
    r = rec.ravel()
    al = acc.tolist(); rl = r.tolist()
    for k in order.tolist():
        t = rl[k]
        if t >= 0:
            al[t] += al[k]
    return np.asarray(al).reshape(zf.shape), order


def hand(zf, z, rec, stream, order):
    """Height above the stream cell each cell drains to (ascending sweep)."""
    base = np.where(stream.ravel(), z.ravel(), np.nan)
    bl = base.tolist(); rl = rec.ravel().tolist(); sl = stream.ravel().tolist(); zl = z.ravel().tolist()
    for k in order[::-1].tolist():           # low to high: receiver first
        if sl[k]:
            continue
        t = rl[k]
        bl[k] = bl[t] if t >= 0 else zl[k]
    return z - np.asarray(bl).reshape(z.shape)


def derivatives(z, res):
    """Zevenbergen-Thorne first and second derivatives; x = east, y = north."""
    p = np.pad(z, 1, mode='edge')
    z1, z2, z3 = p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:]
    z4, z5, z6 = p[1:-1, :-2], p[1:-1, 1:-1], p[1:-1, 2:]
    z7, z8, z9 = p[2:, :-2], p[2:, 1:-1], p[2:, 2:]
    zx = (z6 - z4) / (2 * res)
    zy = (z2 - z8) / (2 * res)                # north minus south
    zxx = (z4 - 2 * z5 + z6) / res ** 2
    zyy = (z2 - 2 * z5 + z8) / res ** 2
    zxy = (-z1 + z3 + z7 - z9) / (4 * res ** 2)
    return zx, zy, zxx, zyy, zxy


def terrain(z, res):
    zx, zy, zxx, zyy, zxy = derivatives(z, res)
    g2 = zx ** 2 + zy ** 2
    slope = np.degrees(np.arctan(np.sqrt(g2)))
    asp = np.arctan2(-zx, -zy)                # downslope azimuth, clockwise from north
    flat = g2 < (np.tan(np.radians(0.05))) ** 2
    with np.errstate(invalid='ignore', divide='ignore'):
        prof = -(zx ** 2 * zxx + 2 * zx * zy * zxy + zy ** 2 * zyy) / (g2 * (1 + g2) ** 1.5)
        plan = -(zy ** 2 * zxx - 2 * zx * zy * zxy + zx ** 2 * zyy) / (g2 ** 1.5)
    prof[flat] = 0; plan[flat] = 0
    tri = np.zeros_like(z)
    p = np.pad(z, 1, mode='edge')
    for di, dj in D8:
        tri += (p[1 + di:1 + di + z.shape[0], 1 + dj:1 + dj + z.shape[1]] - z) ** 2
    tri = np.sqrt(tri)                        # Riley et al. 1999
    return dict(slope=slope, asp_sin=np.where(flat, 0, np.sin(asp)), asp_cos=np.where(flat, 0, np.cos(asp)),
                curv_prof=prof, curv_plan=plan, tri=tri)


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y <= r * r).astype(float)


def focal_mean(a, r, valid=None):
    from scipy.signal import fftconvolve
    m = np.isfinite(a) if valid is None else (valid & np.isfinite(a))
    zz = np.where(m, a, 0.)
    k = disk(r)
    num = fftconvolve(zz, k, 'same'); den = fftconvolve(m.astype(float), k, 'same')
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0.5, num / den, np.nan)


def err(z, size):
    mn = ndi.minimum_filter(z, size, mode='nearest'); mx = ndi.maximum_filter(z, size, mode='nearest')
    me = ndi.uniform_filter(z, size, mode='nearest')
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(mx - mn > 0.01, (me - mn) / (mx - mn), 0.5)
