!
! EarthSciML instrumentation: state leaving INTEGRATE
  IF ( esm_dumping .AND. i == its .AND. j == jts ) THEN
    esm_var_out(:,k) = var(:)
  END IF
