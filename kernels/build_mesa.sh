#!/bin/bash
# Build kernels/mosaic_therm_driver WITH the MESA tape (module_esm_mesa_tape and
# the hooked module_mosaic_therm.F), against WRF's own compiled objects for
# everything else (main/libwrflib.a, the chem170 build), with the model's own
# preprocessing and compiler flags.  The two rebuilt objects are linked ahead
# of the archive, so the archive's module_mosaic_therm.o is never pulled in.
# Run inside wrf-chem-build.sif:
#   apptainer exec -B /projects:/projects ../../apptainer/wrf-chem-build.sif bash build_mesa.sh
set -e
W=$(cd "$(dirname "$0")/.." && pwd)
B=$W/kernels/build_mesa
mkdir -p $B
FL="-O2 -ftree-vectorize -funroll-loops -w -ffree-form -ffree-line-length-none -fconvert=big-endian -frecord-marker=4 -fallow-argument-mismatch -fallow-invalid-boz"
CPPF="-DEM_CORE= -DNMM_CORE= -DNMM_MAX_DIM=2600 -DDA_CORE= -DWRFPLUS= -DIWORDSIZE=4 -DDWORDSIZE=8 -DRWORDSIZE=4 -DLWORDSIZE=4 -DNONSTANDARD_SYSTEM_SUBR -DWRF_USE_CLM -DNETCDF -DHDF5 -DLANDREAD_STUB=1 -DUSE_ALLOCATABLES -Dwrfmodel -DGRIB1 -DINTIO -DKEEP_INT_AROUND -DLIMIT_ARGS -DBUILD_RRTMG_FAST=0 -DBUILD_RRTMK=0 -DBUILD_SBM_FAST=1 -DSHOW_ALL_VARS_USED=0 -DCONFIG_BUF_LEN=65536 -DMAX_DOMAINS_F=21 -DMAX_HISTORY=25 -DNMM_NEST= -DWRF_CHEM -DBUILD_CHEM=1 -DWRF_KPP -I$W/chem -traditional-cpp"
INC="-I$B -I$W/main -I$W/frame -I$W/share -I$W/phys -I$W/chem -I$W/inc -I$W/external/esmf_time_f90 -I/opt/netcdf/include"
cd $B
for f in module_esm_mesa_tape module_mosaic_therm; do
  sed -e "s/^\!.*'.*//" -e "s/^ *\!.*'.*//" $W/chem/$f.F > $f.G
  /lib/cpp -P -nostdinc -I$W/inc $CPPF $f.G > $f.bb
  $W/tools/standard.exe $f.bb | /lib/cpp -P -nostdinc -traditional-cpp > $f.f90
  rm -f $f.G $f.bb
  gfortran $FL $INC -J$B -c $f.f90 -o $f.o
done
gfortran $FL $INC -J$B -c $W/kernels/mosaic_therm_driver.F90 -o mosaic_therm_driver.o
gfortran $FL -c $W/kernels/stub_therm_link.F90 -o stub_therm_link.o
gfortran $FL -o mosaic_therm_driver mosaic_therm_driver.o module_esm_mesa_tape.o module_mosaic_therm.o \
  stub_therm_link.o $W/main/libwrflib.a \
  $W/external/fftpack/fftpack5/libfftpack.a $W/external/io_grib1/libio_grib1.a \
  $W/external/io_grib_share/libio_grib_share.a $W/external/io_int/libwrfio_int.a \
  $W/external/esmf_time_f90/libesmf_time.a $W/frame/module_internal_header_util.o $W/frame/pack_utils.o \
  -L$W/external/io_netcdf -lwrfio_nf -L/opt/netcdf/lib -lnetcdff -lnetcdf \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lm -lz
ls -la mosaic_therm_driver
