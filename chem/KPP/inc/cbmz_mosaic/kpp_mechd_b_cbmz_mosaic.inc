!
! EarthSciML instrumentation: decide once per call whether this call is dumped
  esm_dumping = esm_dump_want('cbmz_kpp')
  IF ( esm_dumping ) THEN
    ALLOCATE( esm_var_in(NVAR,kts:kte), esm_var_out(NVAR,kts:kte), esm_fix(NFIX,kts:kte), &
              esm_rconst(NREACT,kts:kte), esm_jv(njv,kts:kte), esm_temp(kts:kte), esm_cm(kts:kte), &
              esm_ch2o(kts:kte), esm_p(kts:kte), esm_rho(kts:kte), esm_qv(kts:kte) )
    esm_var_in = 0._dp; esm_var_out = 0._dp; esm_fix = 0._dp; esm_rconst = 0._dp; esm_jv = 0._dp
    esm_temp = 0._dp; esm_cm = 0._dp; esm_ch2o = 0._dp; esm_p = 0._dp; esm_rho = 0._dp; esm_qv = 0._dp
  END IF
