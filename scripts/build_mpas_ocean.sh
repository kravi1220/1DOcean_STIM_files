#!/bin/bash
# Build the MPAS-Ocean executable with the ADC closure (models/mpas_adc -> models/mpas_adc/ocean_model).
# Needs: MPI Fortran/C compilers (mpif90, mpicc), NetCDF-C/Fortran, and ParallelIO 2.6.2 built without PnetCDF.
# Usage: NETCDF=<netcdf-fortran prefix> PIO=<pio install prefix> scripts/build_mpas_ocean.sh
set -e
: "${NETCDF:?set NETCDF to the netCDF-Fortran install prefix}"
: "${PIO:?set PIO to the ParallelIO install prefix}"
unset PNETCDF
cd "$(dirname "$0")/../models/mpas_adc"
# bash is required: the Makefile's PIO-version probe uses '&>' (dash runs it in the background).
# Serial make: parallel make races on ocn_vmix_adc.mod.
make SHELL=/bin/bash gfortran CORE=ocean USE_PIO2=true
ls -l ocean_model
