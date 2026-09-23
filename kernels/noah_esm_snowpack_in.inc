! SNOWPACK's ARGUMENTS at entry, published because they are not recoverable
! from the driver-level dump: SFLX has already reset SNDENS to SNEQV/SNOWH and,
! on a snowing step, SNOW_NEW has moved both of them again, so neither the
! dumped SNOWH nor SNEQV/SNOWH is the density this routine receives.  ESD here
! is the pack AFTER sublimation and melt have been taken out of it.
      ESM_SNP_CALLED = 1
      ESM_SNP_SNOWH_IN = SNOWH ; ESM_SNP_SNDENS_IN = SNDENS
      ESM_SNP_ESD = ESD ; ESM_SNP_TSNOW = TSNOW ; ESM_SNP_TSOIL = TSOIL
      ESM_SNP_DTSEC = DTSEC
