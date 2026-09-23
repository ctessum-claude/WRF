! The compacted density and the state behind it.  SNDENS is still the entry
! value here (WRF assigns SNDENS = DSX only below), so recomputing SNDENS*PEXP
! reproduces the unclamped DSX bit for bit.  ESM_SNP_DSX_CLAMP: 0 none,
! 1 capped at 0.40, 2 floored at 0.05.  (0.40 > 0.05, so the second test can
! only fire on a raw value below 0.05 and the two are mutually exclusive.)
      ESM_SNP_TAVGC = TAVGC ; ESM_SNP_ESDCX = ESDCX ; ESM_SNP_BFAC = BFAC
      ESM_SNP_PEXP = PEXP
      ESM_SNP_DSX_RAW = SNDENS * (PEXP)
      ESM_SNP_DSX = DSX
      IF (ESM_SNP_DSX_RAW > 0.40) THEN
         ESM_SNP_DSX_CLAMP = 1
      ELSE IF (ESM_SNP_DSX_RAW < 0.05) THEN
         ESM_SNP_DSX_CLAMP = 2
      ELSE
         ESM_SNP_DSX_CLAMP = 0
      END IF
