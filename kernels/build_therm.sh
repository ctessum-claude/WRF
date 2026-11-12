#!/bin/bash
# Build kernels/mosaic_therm_driver against WRF's own compiled objects (main/libwrflib.a,
# the chem170 build), with the model's own compiler flags.  Run inside wrf-chem-build.sif:
#   apptainer exec -B /projects:/projects ../../apptainer/wrf-chem-build.sif bash build_therm.sh
set -e
W=$(cd "$(dirname "$0")/.." && pwd)
B=$W/kernels/build_therm
mkdir -p $B
FL="-O2 -ftree-vectorize -funroll-loops -w -ffree-form -ffree-line-length-none -fconvert=big-endian -frecord-marker=4 -fallow-argument-mismatch -fallow-invalid-boz"
INC="-I$W/main -I$W/frame -I$W/share -I$W/phys -I$W/chem -I$W/inc -I$W/external/esmf_time_f90 -I/opt/netcdf/include"
cd $B
gfortran $FL $INC -J$B -c $W/kernels/mosaic_therm_driver.F90 -o mosaic_therm_driver.o
gfortran $FL -c $W/kernels/stub_therm_link.F90 -o stub_therm_link.o
gfortran $FL -o mosaic_therm_driver mosaic_therm_driver.o stub_therm_link.o $W/main/libwrflib.a \
  $W/external/fftpack/fftpack5/libfftpack.a $W/external/io_grib1/libio_grib1.a \
  $W/external/io_grib_share/libio_grib_share.a $W/external/io_int/libwrfio_int.a \
  $W/external/esmf_time_f90/libesmf_time.a $W/frame/module_internal_header_util.o $W/frame/pack_utils.o \
  -L$W/external/io_netcdf -lwrfio_nf -L/opt/netcdf/lib -lnetcdff -lnetcdf \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lm -lz
ls -la mosaic_therm_driver
