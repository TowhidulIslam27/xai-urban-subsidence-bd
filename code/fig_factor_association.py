"""Section 5.1, Fig. S (factor association)  Spearman rho of each factor with vU at PS cells, per
city; open markers = dropped in screening. A first look only, not a model result.
Output    manuscript/Figures/FigS_factor_association.png
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import figlib
figlib.init(); figlib.assoc_fig()
