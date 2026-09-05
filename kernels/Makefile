# Standalone CCPP-kernel drivers for the EarthSciML EqWeFiC project.
#   make            real64 (kind_phys = real64) drivers
#   make PREC=      real32 drivers (WRF's default kind_phys)
FC = gfortran
PREC ?= -DDOUBLE_PRECISION
FFLAGS ?= -O2 -ffree-form -ffree-line-length-none -cpp $(PREC)
B = build
SRC_KIND = ../phys/ccpp_kind_types.F
SRC_DUMP = ../phys/module_esm_dump.F
SRC_YSU  = ../phys/physics_mmm/bl_ysu.F90

all: $(B)/ysu_driver

$(B):
	mkdir -p $(B)

$(B)/ccpp_kind_types.o: $(SRC_KIND) | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/module_esm_dump.o: $(SRC_DUMP) | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/esm_flat_input.o: esm_flat_input.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/bl_ysu.o: $(SRC_YSU) $(B)/ccpp_kind_types.o
	$(FC) $(FFLAGS) -J$(B) -I$(B) -c $< -o $@
$(B)/ysu_driver: ysu_driver.F90 $(B)/ccpp_kind_types.o $(B)/module_esm_dump.o $(B)/esm_flat_input.o $(B)/bl_ysu.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) $^ -o $@

clean:
	rm -rf $(B)
