"""Section 3.2, Fig. S (DEM check)  SRTM ground proxy versus the GLO-30 derivatives of the
Chattogram kit: elevation, slope, TWI and profile curvature at the 100 m PS cells (hexbin,
Spearman rho). Uses whichever ground model s32a ran on (FABDEM once exported).
Output    manuscript/Figures/FigS_dem_check.png
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import figlib
figlib.init(); figlib.dem_fig()
