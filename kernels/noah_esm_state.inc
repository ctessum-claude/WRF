! noah_esm_state.inc -- EqWeFiC instrumentation, spliced into a BUILD-TIME COPY
! of phys/module_sf_noahlsm.F by kernels/Makefile (the WRF source is untouched).
! SUBROUTINE SFLX and SUBROUTINE NOPAC compute several quantities that WRF never
! stores in a Registry array and that the driver-level `noah` dump therefore
! cannot carry: the soil thermal conductivity DF1, the surface-energy-balance
! closure pair YY/ZZ1, the linearised Penman quantities RCH/RR/EPSCA/T24, and
! the value of SSOIL that PENMAN is called with (which is the PREVIOUS skin
! temperature's ground heat flux, not the one SHFLX returns).  These module
! variables publish them for kernels/noah_driver.
      REAL :: ESM_DF1 = 0.0, ESM_YY = 0.0, ESM_ZZ1 = 0.0
      REAL :: ESM_SSOIL_PRE = 0.0, ESM_RCH = 0.0, ESM_RR = 0.0
      REAL :: ESM_EPSCA = 0.0, ESM_T24 = 0.0, ESM_FDOWN = 0.0
