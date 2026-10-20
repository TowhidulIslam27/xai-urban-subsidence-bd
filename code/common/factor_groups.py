"""Candidate factors for screening, in priority order, with their group, model set, native
resolution and transfer status (shared by screening, catalogue and figures).

Rules behind the lists (v2, 18 Oct 2026)
  - Model candidates must cover every prediction cell. Factors from the 2024 city land-cover maps
    cover only part of each city, and coverage itself separates core from fringe (different vU),
    so they are kept as LC_CHECK (checks on the covered area), not as candidates.
  - Building factors of the three-city paper exist only at PS cells (PS_ONLY); the Open Buildings
    and GHSL layers from GEE replace them as candidates.
  - City-specific objects (MRT-6, expressway, port, OSM embankments) may enter the single-city
    models but never the transfer test (TRANSFER = False).
  - Fill (M3) means fill history (GSW, Landsat, Dynamic World, GHSL), not present-day ponds.
Factors missing from a table (e.g. GEE layers before the export) are skipped by present().
"""
# (id, group, set, native resolution in m, transfer)
_F = [
    # geology (M1); geo_shared is categorical and handled by the models directly
    ('d_contact', 'Geology', 'M1', 30, True),
    # fill history (M3)
    ('gsw_lost_share', 'Fill', 'M3', 30, True),
    ('gsw_lost_window', 'Fill', 'M3', 30, True),
    ('gsw_lost_after2010', 'Fill', 'M3', 30, True),
    ('fill_age_2023', 'Fill', 'M3', 30, True),
    ('water_loss_1990_2022', 'Fill', 'M3', 30, True),     # the fill process first,
    ('water_share_1990', 'Fill', 'M3', 30, True),         # then the state it started from
    ('ghs_urban_epoch', 'Fill', 'M3', 100, True),
    ('dw_built_change', 'Fill', 'M3', 10, True),
    # groundwater proxies (M3)
    ('pop_density_2020', 'Groundwater proxy', 'M3', 100, True),
    ('viirs_rad', 'Groundwater proxy', 'M3', 500, True),
    # urban load (M2)
    ('ob_bfrac_2023', 'Urban load', 'M2', 4, True),
    ('ob_height_2023', 'Urban load', 'M2', 4, True),
    ('ob_volume_2023', 'Urban load', 'M2', 4, True),
    ('ob_new_frac', 'Urban load', 'M2', 4, True),
    ('ob_dheight_existing', 'Urban load', 'M2', 4, True),
    ('ob_bld_per_ha', 'Urban load', 'M2', 4, True),
    ('ghs_bv_2020', 'Urban load', 'M2', 100, True),
    ('ghs_bh_2018', 'Urban load', 'M2', 100, True),
    ('ghs_bs_d2000_2020', 'Urban load', 'M2', 100, True),
    ('road_dens_500', 'Urban load', 'M2', 10, True),
    ('dist_mainroad', 'Urban load', 'M2', 10, True),
    ('dist_rail', 'Urban load', 'M2', 10, True),
    ('facility_dens_500', 'Urban load', 'M2', 10, True),
    ('dsm_minus_ground', 'Urban load', 'M2', 30, True),
    # land cover (M2)
    ('dw_built', 'Land cover', 'M2', 10, True),
    ('dw_trees', 'Land cover', 'M2', 10, True),
    ('dw_grass', 'Land cover', 'M2', 10, True),
    ('dw_crops', 'Land cover', 'M2', 10, True),
    ('dw_bare', 'Land cover', 'M2', 10, True),
    ('dw_water', 'Land cover', 'M2', 10, True),
    ('dw_flooded_vegetation', 'Land cover', 'M2', 10, True),
    ('s2_ndvi', 'Land cover', 'M2', 10, True),
    ('s2_ndbi', 'Land cover', 'M2', 20, True),
    ('s2_mndwi', 'Land cover', 'M2', 20, True),
    ('s2_bsi', 'Land cover', 'M2', 20, True),
    ('s2_ndvi_sd', 'Land cover', 'M2', 10, True),
    # elevation ranks before the hydrology and terrain derivatives: in a delta city it separates
    # hills, levees and fill, and the derivatives are partly built from it
    ('elev', 'Terrain', 'M0', 30, True),
    # hydrology (M0)
    ('hand', 'Hydrology', 'M0', 30, True),
    ('dist_water', 'Hydrology', 'M0', 30, True),
    ('merit_hand', 'Hydrology', 'M0', 90, True),
    ('twi', 'Hydrology', 'M0', 30, True),
    ('fill_depth', 'Hydrology', 'M0', 30, True),
    ('dist_stream', 'Hydrology', 'M0', 30, True),
    ('drain_dens_1k', 'Hydrology', 'M0', 30, True),
    # soil (M0)
    ('olm_clay', 'Soil', 'M0', 250, True),
    ('olm_sand', 'Soil', 'M0', 250, True),
    ('olm_bulkdens', 'Soil', 'M0', 250, True),
    ('olm_orgc', 'Soil', 'M0', 250, True),
    ('hsg_dual_share_300', 'Soil', 'M0', 250, True),
    ('hsg_poor_share_300', 'Soil', 'M0', 250, True),
    ('smap_sm_annual', 'Soil', 'M0', 9000, True),
    ('chirps_mm_yr', 'Soil', 'M0', 5500, True),
    # terrain (M0)
    ('slope', 'Terrain', 'M0', 30, True),
    ('tpi_300', 'Terrain', 'M0', 30, True),
    ('tpi_1000', 'Terrain', 'M0', 30, True),
    ('relief_1k', 'Terrain', 'M0', 30, True),
    ('tri', 'Terrain', 'M0', 30, True),
    ('err_450', 'Terrain', 'M0', 30, True),
    ('asp_sin', 'Terrain', 'M0', 30, True),
    ('asp_cos', 'Terrain', 'M0', 30, True),
    ('curv_prof', 'Terrain', 'M0', 30, True),
    ('curv_plan', 'Terrain', 'M0', 30, True),
    # city-specific (single-city models only)
    ('dist_embank', 'City-specific', 'M2', 10, False),
    ('dist_mrt', 'City-specific', 'M2', 10, False),
    ('dist_expressway', 'City-specific', 'M2', 10, False),
    ('dist_port', 'City-specific', 'M2', 10, False),
]
PRIORITY = [f[0] for f in _F]                      # process meaning first, generic terrain last
GROUP = {f[0]: f[1] for f in _F}
SET = {f[0]: f[2] for f in _F}
RES = {f[0]: f[3] for f in _F}
TRANSFER = {f[0]: f[4] for f in _F}

# kept in the tables for checks only (partial coverage or PS cells only)
LC_CHECK = ['lc_cover', 'pond_share_300', 'dev_share_300', 'built_share_300', 'green_share_300',
            'water_big_share_500', 'dist_pond', 'dist_dev', 'dist_lc_water']
PS_ONLY = ['bld_density_ps', 'bld_height_ps', 'h23_ob_ps', 'new_share_ps', 'bld_growth_ps']
# Dhaka borehole/well locations (wells_dhaka.py): a survey network with no type, depth or year,
# regular rather than clustered (Clark-Evans R = 1.16) and denser in the old core (rho 0.66 with road
# density); it marks where boreholes were logged, not where water is pumped -> check only
WELL_CHECK = ['dist_well', 'well_dens_1k', 'well_kde_500']
GROUP.update({k: 'Land cover (check)' for k in LC_CHECK}, **{k: 'Urban load (PS only)' for k in PS_ONLY},
             **{k: 'Groundwater (check)' for k in WELL_CHECK})

MAX_RES = 1000          # m; coarser layers have only a handful of values per city and act as location codes
TREND_R2_MAX = 0.7      # factors explained this well by a quadratic surface in E, N act as coordinates


def present(cols, keys=None):
    """Keys (default PRIORITY) that exist in a table's columns, in order."""
    cols = set(cols)
    return [k for k in (keys or PRIORITY) if k in cols]
