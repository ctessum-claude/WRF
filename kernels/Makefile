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

DRIVERS = $(B)/gwdo_driver $(B)/gwdo_diag_driver $(B)/ysu_driver $(B)/sfclayrev_driver $(B)/wsm6_driver $(B)/sw_driver $(B)/slab_driver $(B)/solar_driver $(B)/cldfra_driver
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
$(B)/bl_gwdo.o: $(M)/bl_gwdo.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/module_libmassv.o: $(M)/module_libmassv.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/mp_radar.o: $(M)/mp_radar.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/mp_wsm6.o: $(M)/mp_wsm6.F90 $(B)/ccpp_kind_types.o $(B)/module_libmassv.o $(B)/mp_radar.o $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) -DESM_DUMP -J$(B) -I$(B) -c $< -o $@

# WRF-native modules (default REAL): promoted in the real64 build
$(B)/module_ra_sw.o: $(P)/module_ra_sw.F $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
$(B)/module_sf_slab.o: $(P)/module_sf_slab.F $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@

$(B)/ysu_driver: ysu_driver.F90 $(COMMON) $(B)/bl_ysu.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/gwdo_driver: gwdo_driver.F90 $(COMMON) $(B)/bl_gwdo.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
# gwdo_diag_driver: same kernel with the per-column branch record exposed, for
# comparing the real32 and real64 branch decisions of one dumped row.
$(B)/bl_gwdo_diag.F90: $(M)/bl_gwdo.F90 gwdo_diag_patch.py | $(B)
	python3 gwdo_diag_patch.py $< $@
$(B)/bl_gwdo_diag.o: $(B)/bl_gwdo_diag.F90 $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/gwdo_diag_driver: gwdo_diag_driver.F90 $(COMMON) $(B)/bl_gwdo_diag.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/sfclayrev_driver: sfclayrev_driver.F90 $(COMMON) $(B)/sf_sfclayrev.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/wsm6_driver: wsm6_driver.F90 $(COMMON) $(B)/module_libmassv.o $(B)/mp_radar.o $(B)/mp_wsm6.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@
$(B)/sw_driver: sw_driver.F90 $(COMMON) $(B)/module_ra_sw.o
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@
$(B)/slab_driver: slab_driver.F90 $(COMMON) $(B)/module_sf_slab.o
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@

clean:
	rm -rf $(B)

# solar geometry: radconst + calc_coszen extracted verbatim from the radiation driver
$(B)/solar_kernels.F90: $(P)/module_radiation_driver.F | $(B)
	( echo 'module solar_kernels'; echo 'contains'; \
	  awk '/^   SUBROUTINE radconst\(/,/END SUBROUTINE radconst/' $<; \
	  awk '/^   SUBROUTINE calc_coszen\(/,/END SUBROUTINE calc_coszen/' $<; \
	  echo 'end module solar_kernels' ) > $@
$(B)/wrf_error_stub.o: wrf_error_stub.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/solar_kernels.o: $(B)/solar_kernels.F90 $(B)/wrf_error_stub.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/solar_driver: solar_driver.F90 $(B)/solar_kernels.o $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) $< $(B)/solar_kernels.o $(B)/wrf_error_stub.o $(B)/module_esm_dump.o -o $@

# grid-scale cloud fraction: cal_cldfra1 (icloud = 1) extracted verbatim from the
# radiation driver, against a stub of the four module_state_description ids it reads
$(B)/cldfra_kernels.F90: $(P)/module_radiation_driver.F | $(B)
	( echo 'module cldfra_kernels'; echo 'contains'; \
	  awk '/^   SUBROUTINE cal_cldfra1\(/,/END SUBROUTINE cal_cldfra1/' $<; \
	  echo 'end module cldfra_kernels' ) > $@
$(B)/cldfra_state_stub.o: cldfra_state_stub.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/cldfra_kernels.o: $(B)/cldfra_kernels.F90 $(B)/cldfra_state_stub.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
$(B)/cldfra_driver: cldfra_driver.F90 $(COMMON) $(B)/cldfra_state_stub.o $(B)/cldfra_kernels.o
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@
