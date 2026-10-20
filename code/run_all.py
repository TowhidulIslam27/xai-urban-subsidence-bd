"""Run the whole factor pipeline in manuscript order (about 5 minutes on a laptop).

    python run_all.py            # everything
    python run_all.py --figures  # figures only (tables must exist)
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = ['s31_target_lattice.py', 's32a_terrain_hydrology.py', 's32b_geology.py', 's32c_landcover_fill.py',
        's32d_urban_osm.py', 's32e_ctg_kit_crosscheck.py', 's32e2_soil_hsg.py', 's32m1_wells_dha.py', 's32f_assemble_table.py', 's32g_make_gee_script.py',
        's32h_merge_gee.py', 's32i_factor_catalogue.py', 's32j_screening.py', 's32l_table4_tex.py']
FIGS = ['s21_fig01_study_area.py', 's32k_figS_dem_check.py', 's51_fig03_factor_maps.py', 's51_figS_factor_atlas.py',
        's51_figS_factor_correlation.py', 's51_figS_factor_association.py']
steps = FIGS if '--figures' in sys.argv else DATA + FIGS
for s in steps:
    if s == 's32h_merge_gee.py' and not any(f.endswith('.csv') for f in os.listdir(
            os.path.join(HERE, '..', '..', 'analysis', 'gee', 'exports'))):
        print('skip', s, '(no GEE exports yet)'); continue
    print('==', s, flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, s)], check=True)
