"""Section 5.1, Fig. 3  Six conditioning factors, Dhaka beside Chattogram at one map scale:
ground elevation, HAND, signed distance to the Holocene contact, pond/ditch share, land under
development, road density, clipped to the city corporation boundaries (analysis/City_Maps) on a
dark grey canvas. Shared colour scale per factor; light grey = no source data.
Output    manuscript/Figures/Fig03_factor_maps.png
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import figlib
figlib.init()
figlib.map_sheet(['elev', 'hand', 'd_contact', 'pond_share_300', 'dev_share_300', 'road_dens_500'], 'Fig03_factor_maps',
                 clip=True)
