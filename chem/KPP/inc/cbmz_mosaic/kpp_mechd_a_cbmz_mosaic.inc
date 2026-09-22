!
! EarthSciML instrumentation: write the column dump
  IF ( esm_dumping ) THEN
    CALL esm_dump_open('cbmz_kpp')
    CALL esm_dump_var('its', its);   CALL esm_dump_var('ite', ite)
    CALL esm_dump_var('jts', jts);   CALL esm_dump_var('jte', jte)
    CALL esm_dump_var('kts', kts);   CALL esm_dump_var('kte', kte)
    CALL esm_dump_var('nvar', NVAR); CALL esm_dump_var('nfix', NFIX)
    CALL esm_dump_var('nreact', NREACT); CALL esm_dump_var('njv', njv)
    CALL esm_dump_var('dtstepc', dtstepc)
    CALL esm_dump_var('time_start', TIME_START); CALL esm_dump_var('time_end', TIME_END)
    CALL esm_dump_var('atol', ATOL); CALL esm_dump_var('rtol', RTOL)
    CALL esm_dump_var('icntrl', ICNTRL); CALL esm_dump_var('rcntrl', RCNTRL)
    CALL esm_dump_var('temp', esm_temp); CALL esm_dump_var('c_m', esm_cm); CALL esm_dump_var('c_h2o', esm_ch2o)
    CALL esm_dump_var('p_phy', esm_p); CALL esm_dump_var('rho_phy', esm_rho); CALL esm_dump_var('qv', esm_qv)
    CALL esm_dump_var('jv', esm_jv)
    CALL esm_dump_var('fix', esm_fix)
    CALL esm_dump_var('rconst', esm_rconst)
    CALL esm_dump_var('var_in', esm_var_in)
    CALL esm_dump_var('var_out', esm_var_out)
    CALL esm_dump_close()
    DEALLOCATE( esm_var_in, esm_var_out, esm_fix, esm_rconst, esm_jv, esm_temp, esm_cm, esm_ch2o, esm_p, esm_rho, esm_qv )
  END IF
