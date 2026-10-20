"""Run the whole factor pipeline in manuscript order (about 5 minutes on a laptop).

    python run_all.py            # everything
    python run_all.py --figures  # figures only (tables must exist)
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = ['make_grid.py', 'terrain_hydrology.py', 'geology.py', 'landcover_water.py',
        'osm_factors.py', 'ctg_kit_check.py', 'soil_hsg.py', 'wells_dhaka.py', 'build_table.py', 'make_gee_script.py',
        'merge_gee.py', 'factor_catalogue.py', 'screening.py', 'table4.py']
FIGS = ['fig1_study_area.py', 'fig_dem_check.py', 'fig3_factor_maps.py', 'fig_factor_atlas.py',
        'fig_factor_correlation.py', 'fig_factor_association.py']
steps = FIGS if '--figures' in sys.argv else DATA + FIGS
for s in steps:
    if s == 'merge_gee.py' and not any(f.endswith('.csv') for f in os.listdir(
            os.path.join(HERE, '..', '..', 'analysis', 'gee', 'exports'))):
        print('skip', s, '(no GEE exports yet)'); continue
    print('==', s, flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, s)], check=True)
