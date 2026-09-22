! EarthSciML instrumentation, spliced into a build-time copy of
! physics_mmm/cu_ntiedtke.F90 right after the cumastrn call (kernels/Makefile):
! the convection type (1 deep, 2 shallow, 3 mid-level), cloud base and top
! (kernel levels, counted from the model top), and cumastrn's own tendencies
! and mass fluxes, before cu_ntiedtke_run applies them over dt.
      if (esm_dump_inner) then
        call esm_dump_var('ktype', ktype); call esm_dump_var('kcbot', icbot); call esm_dump_var('kctop', ictop)
        call esm_dump_var('ldcum', locum)
        call esm_dump_var('ztt', ztt); call esm_dump_var('zqq', zqq)
        call esm_dump_var('ptte', ptte); call esm_dump_var('pqte', pqte); call esm_dump_var('pcte', pcte)
        call esm_dump_var('pvom', pvom); call esm_dump_var('pvol', pvol)
        call esm_dump_var('pmfu', zmfu); call esm_dump_var('pmfd', zmfd)
        call esm_dump_var('prsfc', prsfc); call esm_dump_var('pssfc', pssfc)
        call esm_dump_var('scale_fac', scale_fac); call esm_dump_var('scale_fac2', scale_fac2)
      endif
