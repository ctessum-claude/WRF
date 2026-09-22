      ESM_FRH_NLOG = NLOG
      IF (TKELV > (T0- 1.E-3)) THEN
         ESM_FRH_PATH = 1
      ELSE IF (KCOUNT > 0) THEN
         ESM_FRH_PATH = 2
      ELSE
         ESM_FRH_PATH = 3
      END IF
