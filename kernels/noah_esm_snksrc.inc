      IF (K >= 1 .AND. K <= ESM_NSOIL_MAX) THEN
         ESM_FRZ_TAVG(K) = TAVG ; ESM_FRZ_TSNSR(K) = TSNSR
         ESM_FRZ_FREE(K) = FREE ; ESM_FRZ_XH2O(K) = XH2O
         ESM_FRZ_QTOT(K) = QTOT ; ESM_FRZ_CALLED(K) = 1
         ESM_FRZ_NLOG(K) = ESM_FRH_NLOG ; ESM_FRZ_PATH(K) = ESM_FRH_PATH
         ESM_FRZ_XH2O(K) = XH2O
         IF (XH2O == ESM_FRZ_XH2O_RAW(K)) THEN
            ESM_FRZ_CLAMP(K) = 0                ! no clamp: tsnsr = -qtot exactly
         ELSE IF (XH2O == ESM_FRZ_SH2O_PRE(K)) THEN
            ESM_FRZ_CLAMP(K) = 1                ! held at the incoming sh2o: no phase change
         ELSE IF (XH2O == FREE) THEN
            ESM_FRZ_CLAMP(K) = 2                ! stopped at the equilibrium free water
         ELSE IF (XH2O == 0.0) THEN
            ESM_FRZ_CLAMP(K) = 3                ! floored at zero liquid water
         ELSE
            ESM_FRZ_CLAMP(K) = 4                ! capped at smc
         END IF
      END IF
