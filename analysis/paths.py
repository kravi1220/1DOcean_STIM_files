"""Locations of model output, observations and results used by the analysis scripts.

By default everything is relative to the package root; set KG6_ROOT to analyse a copy of
the experiments/ and data/ folders elsewhere.
"""
import os

ROOT = os.environ.get('KG6_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GOTM = os.path.join(ROOT, 'experiments', 'gotm') + '/'        # GOTM runs (all.nc, icepack_column.dat)
MPAS = os.path.join(ROOT, 'experiments', 'mpas_adc') + '/'    # MPAS-ADC cases (hourly/<site>_adc.nc)
LES = os.path.join(ROOT, 'experiments', 'les') + '/'          # Oceananigans output and LES input files
DATA = os.path.join(ROOT, 'data') + '/'                       # CTD casts, OSI SAF extracts
FIG = os.path.join(ROOT, 'results', 'figures') + '/'
TAB = os.path.join(ROOT, 'results', 'tables') + '/'
os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')
