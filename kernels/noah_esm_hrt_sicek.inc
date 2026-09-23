! The same for layers 2..NSOIL.
      IF (K >= 1 .AND. K <= ESM_NSOIL_MAX) THEN
         ESM_FRZ_SICE(K) = SICE ; ESM_FRZ_SH2O_HRT(K) = SH2O (K)
      END IF
