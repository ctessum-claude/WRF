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

DRIVERS = $(B)/ysu_driver $(B)/sfclayrev_driver $(B)/wsm6_driver $(B)/sw_driver $(B)/slab_driver $(B)/ntiedtke_driver $(B)/noah_driver $(B)/cbmz_driver $(B)/mosaic_drydep_driver $(B)/mosaic_subproc_driver
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

# KPP CBM-Z: the WRF build's own preprocessed sources (already kind dp = real64)
C = ../chem
KPPSRC = $(C)/module_kpp_cbmz_mosaic_Precision.f90 $(C)/module_kpp_cbmz_mosaic_Parameters.f90 \
         $(C)/module_kpp_cbmz_mosaic_JacobianSP.f90 $(C)/module_kpp_cbmz_mosaic_Jacobian.f90 \
         $(C)/module_kpp_cbmz_mosaic_Integr.f90 $(C)/module_kpp_cbmz_mosaic_Update_Rconst.f90
$(B)/kpp_cbmz.stamp: $(KPPSRC) | $(B)
	for f in $(KPPSRC); do $(FC) $(FFLAGS) -J$(B) -I$(B) -c $$f -o $(B)/`basename $$f .f90`.o || exit 1; done
	touch $@
KPPOBJ = $(B)/module_kpp_cbmz_mosaic_Precision.o $(B)/module_kpp_cbmz_mosaic_Parameters.o \
         $(B)/module_kpp_cbmz_mosaic_JacobianSP.o $(B)/module_kpp_cbmz_mosaic_Jacobian.o \
         $(B)/module_kpp_cbmz_mosaic_Integr.o $(B)/module_kpp_cbmz_mosaic_Update_Rconst.o

# MOSAIC dry deposition: mosaic_drydep_1clm and aerosol_depvel_2 extracted
# VERBATIM from chem/module_mosaic_drydep.F at build time (lines 260-542), with
# the two WRF-only use statements (module_configure, module_mosaic_driver)
# dropped -- neither routine calls anything from them.  The data modules they
# read are WRF's own and have no dependencies of their own.
$(B)/mosaic_drydep_kernel.F90: $(C)/module_mosaic_drydep.F | $(B)
	( echo 'module mosaic_drydep_kernel'; echo 'contains'; \
	  sed -n '260,542p' $< | sed -e '/use module_configure/d' -e '/use module_mosaic_driver/d'; \
	  echo 'end module mosaic_drydep_kernel' ) > $@
	grep -q 'subroutine aerosol_depvel_2' $@
$(B)/module_peg_util.o: $(C)/module_peg_util.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/stub_wrf_msg.o: stub_wrf_msg.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/module_data_mosaic_asect.o: $(C)/module_data_mosaic_asect.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/module_data_mosaic_other.o: $(C)/module_data_mosaic_other.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/mosaic_drydep_kernel.o: $(B)/mosaic_drydep_kernel.F90 $(B)/module_data_mosaic_asect.o $(B)/module_data_mosaic_other.o $(B)/module_peg_util.o
	$(FC) $(FFLAGS) $(PROMO) -ffixed-line-length-none -J$(B) -I$(B) -c $< -o $@
$(B)/module_state_description.o: ../frame/module_state_description.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/module_esm_mosaic_rates.o: $(C)/module_esm_mosaic_rates.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/module_mosaic_newnuc.o: $(B)/module_esm_mosaic_rates.o $(C)/module_mosaic_newnuc.F $(B)/module_data_mosaic_asect.o $(B)/module_data_mosaic_other.o $(B)/module_peg_util.o $(B)/module_state_description.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
$(B)/module_mosaic_coag.o: $(B)/module_esm_mosaic_rates.o $(C)/module_mosaic_coag.F $(B)/module_data_mosaic_asect.o $(B)/module_data_mosaic_other.o $(B)/module_peg_util.o $(B)/module_state_description.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
MOSSUB = $(B)/module_esm_mosaic_rates.o $(B)/stub_wrf_msg.o $(B)/module_data_mosaic_asect.o $(B)/module_data_mosaic_other.o \
         $(B)/module_peg_util.o $(B)/module_state_description.o \
         $(B)/module_mosaic_newnuc.o $(B)/module_mosaic_coag.o

MOSDD = $(B)/stub_wrf_msg.o $(B)/module_data_mosaic_asect.o $(B)/module_data_mosaic_other.o $(B)/module_peg_util.o $(B)/mosaic_drydep_kernel.o

# WRF-native modules (default REAL): promoted in the real64 build
$(B)/module_ra_sw.o: $(P)/module_ra_sw.F $(B)/module_esm_dump.o
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -I$(B) -c $< -o $@
$(B)/module_model_constants.o: ../share/module_model_constants.F | $(B)
	$(FC) $(FFLAGS) $(PROMO) -J$(B) -c $< -o $@
$(B)/stub_wrf_error.o: stub_wrf_error.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
# Noah: a build-time copy of the module with the three noah_esm_*.inc splices,
# publishing SFLX's and NOPAC's internal DF1/YY/ZZ1 and linearised-Penman
# quantities as module variables (the phys source stays untouched)
$(B)/module_sf_noahlsm_esm.F: $(P)/module_sf_noahlsm.F noah_esm_state.inc noah_esm_penman.inc noah_esm_nopac.inc | $(B)
	sed -e 's|^CONTAINS$$|#include "noah_esm_state.inc"\n&|' \
	    -e 's|^\t\t       ALBEDO,SOLDN,FVB,GAMA,STC(1),ETPN,FLX4,UA_PHYS)$$|&\n#include "noah_esm_penman.inc"|' \
	    -e 's|^      ZZ1 = DF1 / ( -0.5 \* ZSOIL (1) \* RCH \* RR ) + 1.0$$|&\n#include "noah_esm_nopac.inc"|' $< > $@
	grep -q 'noah_esm_state.inc' $@ && grep -q 'noah_esm_penman.inc' $@ && grep -q 'noah_esm_nopac.inc' $@
$(B)/module_sf_noahlsm.o: $(B)/module_sf_noahlsm_esm.F $(B)/module_model_constants.o $(B)/stub_wrf_error.o
	$(FC) $(FFLAGS) $(PROMO) -I. -J$(B) -I$(B) -c $< -o $@
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

$(B)/cbmz_driver: cbmz_driver.F90 $(COMMON) $(B)/kpp_cbmz.stamp
	$(FC) $(FFLAGS) -I$(B) -J$(B) cbmz_driver.F90 $(COMMON) $(KPPOBJ) -o $@

$(B)/mosaic_drydep_driver: mosaic_drydep_driver.F90 $(COMMON) $(MOSDD)
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@

$(B)/mosaic_subproc_driver: mosaic_subproc_driver.F90 $(COMMON) $(MOSSUB)
	$(FC) $(FFLAGS) $(PROMO) -I$(B) -J$(B) $^ -o $@

clean:
	rm -rf $(B)
