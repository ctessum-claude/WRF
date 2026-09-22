!
! EarthSciML instrumentation: state entering INTEGRATE (after Update_Rconst)
  IF ( esm_dumping .AND. i == its .AND. j == jts ) THEN
    esm_var_in(:,k) = var(:)
    esm_fix(:,k)    = fix(:)
    esm_rconst(:,k) = RCONST(:)
    esm_jv(:,k)     = jv(:)
    esm_temp(k)     = TEMP
    esm_cm(k)       = C_M
    esm_ch2o(k)     = fix(indf_H2O)
    esm_p(k)        = REAL( p_phy(i,k,j), KIND=dp )
    esm_rho(k)      = REAL( rho_phy(i,k,j), KIND=dp )
    esm_qv(k)       = REAL( moist(i,k,j,P_QV), KIND=dp )
  END IF
