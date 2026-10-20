"""Section 5.1, Fig. S (factor atlas)  Every candidate factor mapped for both cities in the layout
of Fig. 3: six factors per sheet (three rows of Dhaka | Chattogram pairs), clipped to the city
corporation boundaries, one map scale, north arrow and scale bar on every panel.
Sheets follow the factor groups in priority order; the land-cover check factors (partial
coverage) and the building factors that exist only at PS cells follow on their own sheets.
Output    manuscript/Figures/FigS_factor_atlas_<n>.png (numbered in that order)
"""
import os, sys, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import figlib
from factor_groups import PRIORITY, LC_CHECK, PS_ONLY
from paths import FIGS

PER_SHEET = 6
figlib.init()
cols = set(figlib.T['dha']) | set(figlib.T['ctg']) | set(figlib.M['dha']) | set(figlib.M['ctg'])
blocks = [figlib.display_order([k for k in PRIORITY if k in cols]),
          [k for k in LC_CHECK if k in cols], [k for k in PS_ONLY if k in cols]]
sheets = [b[i:i + PER_SHEET] for b in blocks for i in range(0, len(b), PER_SHEET)]
for f in glob.glob(str(FIGS / 'FigS_factor_atlas_*.png')):     # drop sheets of an earlier layout
    os.remove(f)
for n, keys in enumerate(sheets, 1):
    figlib.map_sheet(keys, f'FigS_factor_atlas_{n}', clip=True)
    print(n, keys)
