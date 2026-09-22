# Standalone kernel drivers for the EarthSciML EqWeFiC project.
#   make            real64 drivers (kind_phys = real64; WRF-native modules
#                   promoted with -fdefault-real-8)
#   make PREC=      real32 drivers (WRF's default kinds)
# Each driver replays one esm dump (flattened by EqWeFiC tools/esm_dump.py
# --flat) and writes an esm dump of the outputs; see the header of each
# *_driver.F90.  rrtm has no driver (rrtminit needs WRF's I/O and error
# modules); use the real32 in-model dumps for it.
FC = gfortran
PREC ?= -DDOUBLE_PRECISION
PROMO = $(if $(PREC),-fdefault-real-8,)
FFLAGS ?= -O2 -ffree-form -ffree-line-length-none -cpp $(PREC)
B = build
P = ../phys
M = ../phys/physics_mmm

DRIVERS = $(B)/ysu_driver $(B)/sfclayrev_driver $(B)/wsm6_driver $(B)/sw_driver $(B)/slab_driver $(B)/ntiedtke_driver $(B)/noah_driver
all: $(DRIVERS)

$(B):
	mkdir -p $(B)

# shared
$(B)/ccpp_kind_types.o: $(P)/ccpp_kind_types.F | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/module_esm_dump.o: $(P)/module_esm_dump.F | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/esm_flat_input.o: esm_flat_input.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
COMMON = $(B)/ccpp_kind_types.o $(B)/module_esm_dump.o $(B)/esm_flat_input.o

# CCPP kernels (kind_phys)
$(B)/bl_ysu.o: $(M)/bl_ysu.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/sf_sfclayrev.o: $(M)/sf_sfclayrev.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/module_libmassv.o: $(M)/module_libmassv.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/mp_radar.o: $(M)/mp_radar.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/mp_wsm6.o: $(M)/mp_wsm6.F90 $(B)/ccpp_kind_types.o $(B)/module_libmassv.o $(B)/mp_radar.o $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) -DESM_DUMP -J$(B) -I$(B) -c $< -o $@

# New Tiedtke: a build-time copy of the kernel with ntiedtke_esm_inner.inc
# spliced in after the cumastrn call (the physics_mmm source stays untouched)
$(B)/cu_ntiedtke_esm.F90: $(M)/cu_ntiedtke.F90 ntiedtke_esm_inner.inc | $(B)
	sed -e 's/^ *use cu_ntiedtke_common$$/&\n use module_esm_dump/' \
	    -e '/^     &     scale_fac, scale_fac2)$$/{n;s/^!$$/#include "ntiedtke_esm_inner.inc"\n!/}' $< > $@
	grep -q 'ntiedtke_esm_inner.inc' $@
$(B)/cu_ntiedtke.o: $(B)/cu_ntiedtke_esm.F90 $(B)/ccpp_kind_types.o $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) -I. -J$(B) -I$(B) -c $< -o $@

# WRF-native modules (default REAL): promoted in the real64 build
$(B)/module_ra_sw.o: $(P)/module_ra_sw.F $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
$(B)/module_model_constants.o: ../share/module_model_constants.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/stub_wrf_error.o: stub_wrf_error.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/module_sf_noahlsm.o: $(P)/module_sf_noahlsm.F $(B)/module_model_constants.o $(B)/stub_wrf_error.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
$(B)/module_sf_slab.o: $(P)/module_sf_slab.F $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@

$(B)/ysu_driver: ysu_driver.F90 $(COMMON) $(B)/bl_ysu.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/sfclayrev_driver: sfclayrev_driver.F90 $(COMMON) $(B)/sf_sfclayrev.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/wsm6_driver: wsm6_driver.F90 $(COMMON) $(B)/module_libmassv.o $(B)/mp_radar.o $(B)/mp_wsm6.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/sw_driver: sw_driver.F90 $(COMMON) $(B)/module_ra_sw.o
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@
$(B)/slab_driver: slab_driver.F90 $(COMMON) $(B)/module_sf_slab.o
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@

$(B)/ntiedtke_driver: ntiedtke_driver.F90 $(COMMON) $(B)/cu_ntiedtke.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@

$(B)/noah_driver: noah_driver.F90 $(COMMON) $(B)/module_model_constants.o $(B)/stub_wrf_error.o $(B)/module_sf_noahlsm.o
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@

clean:
	rm -rf $(B)
