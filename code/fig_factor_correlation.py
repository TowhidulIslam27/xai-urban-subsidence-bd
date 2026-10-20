"""Section 5.1, Fig. S (factor correlation)  Spearman correlation between candidate factors at
PS cells, Dhaka and Chattogram; grey labels = dropped in screening (needs screening.py).
Output    manuscript/Figures/FigS_factor_correlation.png
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import figlib
figlib.init(); figlib.corr_fig()
