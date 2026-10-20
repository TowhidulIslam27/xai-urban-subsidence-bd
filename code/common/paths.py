"""Folder layout of the paper (all paths relative to XAI_Paper).

    XAI_Paper/
      InSAR_data/            target (vU per 100 m cell), factors_from_paper, basedata, validation
      Dataset/               Chattogram GLO-30 kit, soil moisture, geology
      analysis/inputs/       OSM lines (sqlite), DEMs/ (SRTM 1"), 3_City_LC_2024/
      analysis/factors/      outputs: lattice, factor parts, <c>_factors.csv, <c>_maps.npz
      analysis/screening/    outputs: factor catalogue and screening tables
      analysis/gee/exports/  GEE exports returned from Google Drive
      manuscript/code/       these scripts
      manuscript/Figures/    figures (PNG, 600 dpi)

Override the root with the environment variable XAI_ROOT."""
import os
from pathlib import Path

ROOT = Path(os.environ.get('XAI_ROOT', Path(__file__).resolve().parents[3]))
INSAR = ROOT / 'InSAR_data'
DATASET = ROOT / 'Dataset'
INPUTS = ROOT / 'analysis' / 'inputs'
DEMS = INPUTS / 'DEMs'
LC = INPUTS / '3_City_LC_2024'
FACTORS = ROOT / 'analysis' / 'factors'
PARTS = FACTORS / 'parts'
SCREEN = ROOT / 'analysis' / 'screening'
GEE_EXPORTS = ROOT / 'analysis' / 'gee' / 'exports'
CITY_MAPS = ROOT / 'analysis' / 'City_Maps'      # city corporation boundaries and wards (GeoPackage)
CODE = ROOT / 'manuscript' / 'code'
FIGS = ROOT / 'manuscript' / 'Figures'
for d in (FACTORS, PARTS, SCREEN, GEE_EXPORTS, FIGS):
    d.mkdir(parents=True, exist_ok=True)
CITIES = ['dha', 'ctg']
