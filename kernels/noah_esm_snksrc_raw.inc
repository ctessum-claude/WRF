! Captured BEFORE SNKSRC's clamps, so the sink is referenceable without
! differencing: XH2O here is the raw energy-driven value, SH2O is still the
! pre-update liquid water, and DZ is the layer thickness the conversion uses.
! FORTRAN_BUGS N97: without these, tsnsr can only be recovered from
! (XH2O - SH2O), which cancellation destroys at the small dt a dt -> 0
! reference needs, while a larger dt lets the clamps bind instead.
      IF (K >= 1 .AND. K <= ESM_NSOIL_MAX) THEN
         ESM_FRZ_XH2O_RAW(K) = XH2O ; ESM_FRZ_SH2O_PRE(K) = SH2O
         ESM_FRZ_DZ(K) = DZ
      END IF
