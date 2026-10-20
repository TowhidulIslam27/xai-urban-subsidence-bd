"""Section 3.2, Table 3  Factor catalogue: one row per conditioning factor, with source, status per city,
model set (M0-M3; '-' = check only or not modelled), transfer eligibility and flags.
v2 (18 Oct 2026): land-cover factors and PS-only building factors are check only; distances
capped at 2 km; embankment and other single-city objects are city-specific (no transfer).
Writes analysis/screening/factor_catalogue.csv."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import pandas as pd

B, G, LC, PS, NA, LOC = 'built', 'GEE', 'built (CTG kit)', 'PS cells only', 'not available', 'local file'
rows = [
 # id, group, set, factor, source, res, method, dha, ctg, transfer, flag
 ('elev', 'Terrain', 'M0', 'Ground elevation', 'FABDEM (GEE); SRTM 1" proxy until then', '30 m', 'FABDEM 3x3 mean; interim: SRTM 20th pct in 270 m', B, B, 'yes', 'SRTM proxy cuts CTG ridges (7 m low on slopes >= 10 deg): provisional'),
 ('slope', 'Terrain', 'M0', 'Slope', 'same DEM', '30 m', 'Zevenbergen-Thorne', B, B, 'yes', 'DEM-dependent: rho 0.38 vs GLO-30 at CTG PS cells'),
 ('asp_sin', 'Terrain', 'M0', 'Aspect, eastness (sin)', 'same DEM', '30 m', 'downslope azimuth; 0 where flat', B, B, 'yes', 'noise on flat ground'),
 ('asp_cos', 'Terrain', 'M0', 'Aspect, northness (cos)', 'same DEM', '30 m', 'as above', B, B, 'yes', 'noise on flat ground'),
 ('curv_prof', 'Terrain', 'M0', 'Profile curvature', 'same DEM', '30 m', 'Zevenbergen-Thorne', B, B, 'yes', 'poor agreement between DEMs (rho 0.15 at PS cells)'),
 ('curv_plan', 'Terrain', 'M0', 'Plan curvature', 'same DEM', '30 m', 'Zevenbergen-Thorne', B, B, 'yes', 'poor agreement between DEMs (rho 0.07 at PS cells)'),
 ('tri', 'Terrain', 'M0', 'Terrain ruggedness index', 'same DEM', '30 m', 'Riley et al. 1999', B, B, 'yes', ''),
 ('tpi_300', 'Terrain', 'M0', 'Topographic position, 300 m', 'same DEM', '30 m', 'z minus 300 m disk mean', B, B, 'yes', ''),
 ('tpi_1000', 'Terrain', 'M0', 'Topographic position, 1 km', 'same DEM', '30 m', 'z minus 1 km disk mean', B, B, 'yes', ''),
 ('relief_1k', 'Terrain', 'M0', 'Local relief, 1 km', 'same DEM', '30 m', 'max - min in 990 m window', B, B, 'yes', ''),
 ('err_450', 'Terrain', 'M0', 'Elevation-relief ratio', 'same DEM', '30 m', '(mean-min)/(max-min), 450 m', B, B, 'yes', ''),
 ('twi', 'Hydrology', 'M0', 'Topographic wetness index', 'same DEM', '30 m', 'priority-flood fill, D8, ln(a/tan b)', B, B, 'yes', 'rho 0.30 vs GLO-30 at PS cells; routing on a city surface is approximate'),
 ('hand', 'Hydrology', 'M0', 'Height above nearest drainage', 'same DEM', '30 m', 'D8 to streams >= 0.5 km2; clipped at 0', B, B, 'yes', 'derived streams, not the mapped khals; MERIT HAND as check'),
 ('merit_hand', 'Hydrology', 'M0', 'HAND (MERIT Hydro)', 'MERIT/Hydro/v1_0_1', '90 m', 'mean per cell', G, G, 'yes', ''),
 ('fill_depth', 'Hydrology', 'M0', 'Depth of closed depressions', 'same DEM', '30 m', 'filled minus ground', B, B, 'yes', ''),
 ('dist_stream', 'Hydrology', 'M0', 'Distance to derived stream', 'same DEM', '30 m', 'Euclidean', B, B, 'yes', ''),
 ('drain_dens_1k', 'Hydrology', 'M0', 'Drainage density', 'same DEM', '30 m', 'stream km per km2 in 1 km disk', B, B, 'yes', ''),
 ('dist_water', 'Hydrology', 'M0', 'Distance to permanent surface water', 'JRC GSW occurrence >= 50 % (GEE)', '30 m', 'Euclidean, capped 2 km', G, G, 'yes', 'replaces LC + Natural Earth distance to river'),
 ('dist_lc_water', 'Land cover (check)', '-', 'Distance to mapped river/canal/lake', 'City LC 2024', '10 m', 'censored at the map edge, capped 2 km', B, B, 'no', 'CHECK ONLY: partial coverage'),
 ('water_big_share_500', 'Land cover (check)', '-', 'Open-water share, 500 m', 'City LC 2024', '10 m', 'share of mapped area; NaN if < 50 % mapped', B, B, 'no', 'CHECK ONLY: map misses 28 % of Dhaka PS cells'),
 ('dist_coast', 'Hydrology', '-', 'Distance to coastline', 'CTG kit F24', '30 m', '', NA, LC, 'no', 'constant in Dhaka; CTG-only'),
 ('geo_unit', 'Geology', 'M1', 'Geological unit', 'GSB 1:1M (paper polygons)', '1:1M', 'majority per cell', B, B, 'shared classes', '99 % agreement with paper units'),
 ('geo_shared', 'Geology', 'M1', 'Shared geology class', 'derived', '-', 'holocene fine/marsh/sand, older, channel 1990', B, B, 'yes', ''),
 ('d_contact', 'Geology', 'M1', 'Signed distance to Holocene/older contact', 'GSB 1:1M', '30 m', '+ inside Holocene; capped +/-2 km', B, B, 'yes', 'uncapped it acted as a coordinate (trend R2 0.8)'),
 ('young_share', 'Geology', 'M1', 'Holocene share of cell', 'GSB 1:1M', '30 m', '', B, B, 'yes', ''),
 ('dist_fault', 'Geology', '-', 'Distance to active fault', 'GEM (CTG kit F14)', '30 m', '', NA, LC, 'no', 'drop: coordinate proxy (trend R2 0.97)'),
 ('olm_clay', 'Soil', 'M0', 'Clay, 0-200 cm mean', 'OpenLandMap v02', '250 m', 'trapezoid depth mean', G, G, 'yes', ''),
 ('olm_sand', 'Soil', 'M0', 'Sand, 0-200 cm mean', 'OpenLandMap v02', '250 m', '', G, G, 'yes', ''),
 ('olm_bulkdens', 'Soil', 'M0', 'Bulk density', 'OpenLandMap v02', '250 m', '', G, G, 'yes', ''),
 ('olm_orgc', 'Soil', 'M0', 'Organic carbon (peat proxy)', 'OpenLandMap v02', '250 m', '', G, G, 'yes', 'compressible organic soils'),
 ('olm_texture_30cm', 'Soil', 'M0', 'USDA texture class, 30 cm', 'OpenLandMap v02', '250 m', 'mode', G, G, 'yes', ''),
 ('clay', 'Soil', 'M0', 'Clay, weighted mean', 'SoilGrids (CTG kit F16)', '250 m', '', NA, LC, 'cross-check', 'replaced by OpenLandMap for both cities'),
 ('hsg_dual_share_300', 'Soil', 'M0', 'High-water-table (dual) soils and poorly drained soils, share within 300 m', 'HYSOGs250m (Dataset)', '250 m', 'nearest pixel per cell, share of valid pixels', B, B, 'yes', 'same source as CTG kit F15 (98.6 % agreement)'),
 ('bsr', 'Soil', '-', 'BSR (CTG kit F30)', 'CTG kit', '30 m', '', NA, LC, 'no', 'definition to confirm'),
 ('smap_sm_annual', 'Moisture/climate', 'M0', 'Soil moisture, annual/wet/dry/root zone', 'SMAP L4 (SPL4SMGP)', '9-11 km', '2018-2023 means', NA, NA, 'no', 'not exported: native resolution > 1 km'),
 ('sm_annual', 'Moisture/climate', 'M0', 'Soil moisture (CTG kit)', 'CTG kit', '500 m', '', NA, LC, 'cross-check', 'range 0.235-0.250: near constant'),
 ('chirps_mm_yr', 'Moisture/climate', '-', 'Annual rainfall', 'CHIRPS daily', '5 km', '2018-2022 mean', NA, NA, 'no', 'not exported: native resolution > 1 km'),
 ('lst_dry', 'Moisture/climate', '-', 'Dry-season land surface temperature', 'Landsat 8/9 C2 L2', '30 m (100 m TIR)', 'median', G, G, 'yes', 'context only, not a process factor'),
 ('built_share_300', 'Land cover (check)', '-', 'Built-up share, 300 m', 'City LC 2024', '10 m', 'NaN if < 50 % of the disk mapped', B, B, 'no', 'CHECK ONLY: partial coverage; dw_built is the candidate'),
 ('green_share_300', 'Land cover (check)', '-', 'Green/open share, 300 m', 'City LC 2024', '10 m', 'NaN if < 50 % of the disk mapped', B, B, 'no', 'CHECK ONLY: partial coverage'),
 ('dw_built', 'Land cover', 'M2', 'Dynamic World class probabilities (7)', 'GOOGLE/DYNAMICWORLD/V1', '10 m', '2017-10..2023-09 mean', G, G, 'yes', ''),
 ('s2_ndvi', 'Land cover', 'M2', 'NDVI, NDBI, MNDWI, BSI, NDVI s.d.', 'Sentinel-2 SR, Cloud Score+', '10-20 m', 'dry-season median 2018-2023', G, G, 'yes', ''),
 ('bld_density_ps', 'Urban load (PS only)', '-', 'Building footprint density', 'Open Buildings v3 (paper)', '0.5 m', '', PS, PS, 'yes', 'needs full grid: ob_bfrac_2023'),
 ('bld_height_ps', 'Urban load (PS only)', '-', 'Building height', 'GHSL 2018 (paper)', '100 m', '', PS, PS, 'yes', 'needs full grid: ghs_bh_2018'),
 ('h23_ob_ps', 'Urban load (PS only)', '-', 'Building height 2023', 'Open Buildings 2.5D Temporal (paper)', '4 m', '', PS, PS, 'yes', 'needs full grid: ob_height_2023'),
 ('new_share_ps', 'Urban load (PS only)', '-', 'New-construction share 2017-2023', 'Open Buildings 2.5D Temporal (paper)', '4 m', '', PS, PS, 'yes', 'needs full grid: ob_new_frac'),
 ('bld_growth_ps', 'Urban load (PS only)', '-', 'Built-load growth', 'Open Buildings 2.5D Temporal (paper)', '4 m', '', PS, PS, 'yes', ''),
 ('ob_bfrac_2023', 'Urban load', 'M2', 'Footprint share 2016/2023, new share, height, volume, height change of existing buildings, buildings/ha', 'GOOGLE/Research/open-buildings-temporal/v1', '4 m eff.', 'mean per cell', G, G, 'yes', 'replaces the PS-only columns'),
 ('ghs_bv_2020', 'Urban load', 'M2', 'Built volume 2020, height 2018, surface 1990/2020', 'GHSL P2023A', '100 m', '', G, G, 'yes', ''),
 ('road_dens_500', 'Urban load', 'M2', 'Road density, 500 m', 'OSM (paper extract)', '10 m', 'km per km2', B, B, 'yes', ''),
 ('dist_mainroad', 'Urban load', 'M2', 'Distance to main road', 'OSM', '10 m', 'capped 2 km', B, B, 'yes', ''),
 ('dist_rail', 'Urban load', 'M2', 'Distance to railway', 'OSM', '10 m', 'capped 2 km', B, B, 'yes', ''),
 ('facility_dens_500', 'Urban load', 'M2', 'Facility density (schools, health, power)', 'OSM', 'points', 'per km2, 500 m', B, B, 'yes', 'OSM completeness differs by city'),
 ('dist_embank', 'City-specific', 'M2', 'Distance to mapped embankment', 'OSM', '10 m', 'capped 2 km; railway-named features left out', B, B, 'no', 'not the same object: CTG = Outer Ring Road coastal embankment, DHA = a few beribadh roads'),
 ('dist_mrt', 'City-specific', 'M2', 'Distance to MRT line 6', 'OSM', '10 m', 'capped 2 km', B, NA, 'no', 'Dhaka only; built during the InSAR window'),
 ('dist_expressway', 'City-specific', 'M2', 'Distance to elevated expressway', 'OSM', '10 m', 'capped 2 km', B, NA, 'no', 'Dhaka only'),
 ('dist_port', 'City-specific', 'M2', 'Distance to port', 'OSM', '10 m', 'capped 2 km', NA, B, 'no', 'CTG only'),
 ('pond_share_300', 'Land cover (check)', '-', 'Pond/ditch share, 300 m (2024)', 'City LC 2024', '10 m', 'NaN if < 50 % of the disk mapped', B, B, 'no', 'CHECK ONLY: surviving ponds, not fill history; partial coverage'),
 ('dist_pond', 'Land cover (check)', '-', 'Distance to pond/ditch', 'City LC 2024', '10 m', 'censored at the map edge, capped 2 km', B, B, 'no', 'CHECK ONLY: partial coverage'),
 ('dev_share_300', 'Land cover (check)', '-', 'Land under development, 300 m', 'City LC 2024', '10 m', 'NaN if < 50 % of the disk mapped', B, B, 'no', 'CHECK ONLY: 2024 status, after the InSAR window; partial coverage'),
 ('dist_dev', 'Land cover (check)', '-', 'Distance to land under development', 'City LC 2024', '10 m', 'censored at the map edge, capped 2 km', B, B, 'no', 'CHECK ONLY: partial coverage'),
 ('gsw_fill_year', 'Fill', 'M3', 'Year water was lost (fill year) and fill age', 'JRC GSW 1.4 YearlyHistory', '30 m', 'last wet year + 1, where dry 2019-21', G, G, 'yes', 'starts 1984; lost water = drained or filled'),
 ('gsw_lost_window', 'Fill', 'M3', 'Share of cell that lost water in 2017-2019 (inside the InSAR window)', 'JRC GSW 1.4 YearlyHistory', '30 m', 'lost and last wet year 2016-2018', G, G, 'yes', 'later losses not visible: dry 2019-2021 required, GSW ends 2021'),
 ('gsw_lost_share', 'Fill', 'M3', 'Share of cell converted from water; after 2010; ever water; occurrence', 'JRC GSW 1.4', '30 m', '', G, G, 'yes', ''),
 ('water_share_1990', 'Fill', 'M3', 'Water share 1990, 2005, 2022 and loss 1990-2022 (Landsat 5/8 MNDWI)', 'LANDSAT C2 L2', '30 m', 'dry-season median, same method each epoch', G, G, 'yes', 'independent check on GSW; scene counts exported'),
 ('ghs_urban_epoch', 'Fill', 'M3', 'First epoch >= 20 % built (urbanisation age)', 'GHSL P2023A built surface', '100 m', '1975-2020, 5-yr; earliest per cell; never built = NaN', G, G, 'yes', 'with ghs_urban_by2020 (share urban by 2020)'),
 ('dw_built_change', 'Fill', 'M3', 'Built probability 2018 vs 2023', 'Dynamic World', '10 m', '', G, G, 'yes', ''),
 ('pop_density_2020', 'Groundwater proxy', 'M3', 'Population density (pumping demand)', 'WorldPop 100 m', '100 m', '', G, G, 'yes', 'PROXY - not groundwater'),
 ('viirs_rad', 'Groundwater proxy', 'M3', 'Night-light radiance (economic/industrial demand)', 'VIIRS DNB monthly', '500 m', '2017-2023 median', G, G, 'yes', 'PROXY - not groundwater'),
 ('gw_head', 'Groundwater', 'M3', 'Water-table depth / head decline', 'BWDB piezometers, DWASA/CWASA', '-', '', NA, NA, 'yes', 'NOT AVAILABLE: GRACE/GLDAS too coarse'),
 ('gw_wells', 'Groundwater', 'M3', 'Production-well density', 'DWASA/CWASA', '-', '', NA, NA, 'yes', 'NOT AVAILABLE: not public'),
 ('clay_thickness', 'Geology', 'M3', 'Soft clay / fill thickness, depth to Pleistocene', 'CDMP/UDMP boreholes', '-', '', NA, NA, 'yes', 'NOT AVAILABLE in GEE'),
 ('foundation', 'Urban load', '-', 'Foundation type/depth', '-', '-', '', NA, NA, '-', 'NOT AVAILABLE'),
 ('fill_pre1984', 'Fill', 'M3', 'Fill before 1984', 'old maps, CORONA', '-', '', NA, 'partial (1990 channel)', '-', 'NOT AVAILABLE in GEE'),
]
cols = ['id', 'group', 'set', 'factor', 'source', 'res', 'method', 'dha', 'ctg', 'transfer', 'flag']
CAT = pd.DataFrame(rows, columns=cols)

if __name__ == '__main__':
    from paths import SCREEN
    CAT.to_csv(SCREEN / 'factor_catalogue.csv', index=False)
    print('written', SCREEN / 'factor_catalogue.csv')
    print(CAT.groupby('group').size().to_string()); print(len(CAT))
    print(pd.crosstab(CAT.dha, CAT.ctg))
