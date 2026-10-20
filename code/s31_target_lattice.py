"""Section 3.1 / 4.1  InSAR target and the 100 m analysis lattice.

Reads the vertical-rate cells of the three-city decomposition (InSAR_data/target/<c>_cells_plus.csv),
fits the regular lon/lat lattice of those cells (extended to cover the city corporation outline in
analysis/City_Maps) and defines the prediction area: every cell within 1 km of a cell with PS in both
tracks, plus every cell inside the city outline.

Outputs  analysis/factors/<c>_grid.json   lattice definition (also used by the GEE export script)
         analysis/factors/parts/<c>_lattice.npz   prediction-area mask, PS mask, vU on the lattice
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import json, numpy as np
from lattice import Ctx, save_part
from paths import FACTORS, CITIES


def run(c):
    x = Ctx(c); L = x.L
    json.dump(L.meta(), open(FACTORS / f'{c}_grid.json', 'w'), indent=1)
    vU = np.full(L.nr * L.nc, np.nan); vU[L.flat(x.cells.row.values, x.cells.col.values)] = x.cells.vU
    save_part(c, 'lattice', dict(aoi=x.aoi, has_ps=x.has_ps, near_ps=x.near_ps, in_city=x.in_city,
                                 vU=vU.reshape(L.nr, L.nc)))
    x.log('lattice', L.nr, 'x', L.nc, '| prediction cells', int(x.aoi.sum()), '| near PS', int(x.near_ps.sum()),
          '| in city', int(x.in_city.sum()), '| in city, > 1 km from PS', int((x.in_city & ~x.near_ps).sum()),
          '| cells with vU', int(x.has_ps.sum()))


if __name__ == '__main__':
    for c in sys.argv[1:] or CITIES:
        run(c)
