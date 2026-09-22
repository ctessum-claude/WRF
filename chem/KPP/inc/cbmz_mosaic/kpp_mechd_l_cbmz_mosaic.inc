!
! EarthSciML instrumentation (column i = its, j = jts; level index k)
  LOGICAL :: esm_dumping
  REAL( KIND = dp ), ALLOCATABLE, DIMENSION(:,:) :: esm_var_in, esm_var_out, esm_fix, esm_rconst, esm_jv
  REAL( KIND = dp ), ALLOCATABLE, DIMENSION(:)   :: esm_temp, esm_cm, esm_ch2o, esm_p, esm_rho, esm_qv
