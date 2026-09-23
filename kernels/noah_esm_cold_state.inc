! noah_esm_cold_state.inc -- EqWeFiC instrumentation for the COLD-SEASON paths,
! spliced into a BUILD-TIME COPY of phys/module_sf_noahlsm.F by kernels/Makefile
! (the WRF source is untouched).  SNOPAC's surface-energy-balance closure and
! the frozen-soil sink (HRT -> TBND/TMPAVG/SNKSRC -> FRH2O) compute quantities
! that WRF never stores in a Registry array, so the driver-level `noah` dump
! cannot carry them.  These module variables publish them for
! kernels/noah_driver, which zeroes them before each SFLX call.
      REAL :: ESM_SNO_T12 = 0.0, ESM_SNO_T12A = 0.0, ESM_SNO_T12B = 0.0
      REAL :: ESM_SNO_DENOM = 0.0, ESM_SNO_DTOT = 0.0, ESM_SNO_DSOIL = 0.0
      REAL :: ESM_SNO_DF1 = 0.0, ESM_SNO_YY = 0.0, ESM_SNO_ZZ1 = 0.0
      REAL :: ESM_SNO_ETANRG = 0.0, ESM_SNO_ESNOW1 = 0.0, ESM_SNO_ESNOW2 = 0.0
      REAL :: ESM_SNO_EX = 0.0
      REAL :: ESM_SNO_PRCP1 = 0.0, ESM_SNO_SNDENS = 0.0, ESM_SNO_ESD = 0.0
      REAL :: ESM_SNO_SNOMLT = 0.0, ESM_SNO_FLX1 = 0.0, ESM_SNO_FLX3 = 0.0
!     0 = SNOPAC not called, 1 = sub-freezing block (T12 <= TFREEZ), 2 = melt block
      INTEGER :: ESM_SNO_BRANCH = 0
!     frozen-soil sink, one entry per soil layer
      INTEGER, PARAMETER :: ESM_NSOIL_MAX = 10
      REAL :: ESM_FRZ_TAVG(ESM_NSOIL_MAX)  = 0.0
      REAL :: ESM_FRZ_TBND(ESM_NSOIL_MAX)  = 0.0
      REAL :: ESM_FRZ_TSNSR(ESM_NSOIL_MAX) = 0.0
      REAL :: ESM_FRZ_FREE(ESM_NSOIL_MAX)  = 0.0
      REAL :: ESM_FRZ_XH2O(ESM_NSOIL_MAX)  = 0.0
      REAL :: ESM_FRZ_QTOT(ESM_NSOIL_MAX)  = 0.0
      REAL :: ESM_FRZ_SH2O_IN(ESM_NSOIL_MAX) = 0.0
      INTEGER :: ESM_FRZ_CALLED(ESM_NSOIL_MAX) = 0
      INTEGER :: ESM_FRZ_NLOG(ESM_NSOIL_MAX)   = 0
!     FRH2O path: 1 = T above T0-1e-3 (FREE = SMC), 2 = Newton converged,
!                 3 = Flerchinger explicit fallback (Newton hit its 10-iteration cap)
      INTEGER :: ESM_FRZ_PATH(ESM_NSOIL_MAX)   = 0
      REAL :: ESM_FRZ_TSURF = 0.0
!     scratch written by FRH2O, copied into the per-layer slots by SNKSRC
      INTEGER :: ESM_FRH_NLOG = 0, ESM_FRH_PATH = 0
      REAL :: ESM_SNO_EX_RAW = 0.0, ESM_SNO_ESD_IN = 0.0
      INTEGER :: ESM_SNO_MELT_BLOCK = 0
      REAL, DIMENSION(ESM_NSOIL_MAX) :: ESM_FRZ_XH2O_RAW = 0.0
      REAL, DIMENSION(ESM_NSOIL_MAX) :: ESM_FRZ_SH2O_PRE = 0.0
      REAL, DIMENSION(ESM_NSOIL_MAX) :: ESM_FRZ_DZ = 0.0
      INTEGER, DIMENSION(ESM_NSOIL_MAX) :: ESM_FRZ_CLAMP = -1
      REAL, DIMENSION(ESM_NSOIL_MAX) :: ESM_FRZ_SICE = 0.0
      REAL, DIMENSION(ESM_NSOIL_MAX) :: ESM_FRZ_SH2O_HRT = 0.0
