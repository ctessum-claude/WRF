! solar_driver: evaluate WRF's radconst + calc_coszen for synthetic regimes that
! the em_scm_xy reference run never reaches (the julian < 80 branch of radconst,
! radt > 0, multi-day xtime, southern/eastern hemisphere, polar day), writing one
! esm dump per regime. The two subroutines are EXTRACTED VERBATIM from
! phys/module_radiation_driver.F at build time (Makefile rule solar_kernels.F90),
! so the driver cannot drift from the model code. DEGRAD and DPD are the
! share/module_model_constants.F:72-74 definitions. Build real32 (PREC=) to
! reproduce WRF's default REAL.
!   solar_driver <output_dir>
program solar_driver
  use module_esm_dump
  use solar_kernels, only: radconst, calc_coszen
  implicit none
  integer, parameter :: nreg = 6
  real, parameter :: piconst = 3.1415926535897932384626433
  real, parameter :: degrad = piconst/180., dpd = 360./365.
  character(len=32) :: names(nreg)
  real :: julian(nreg), xtime(nreg), radt(nreg), gmt(nreg), lat(nreg), lon(nreg)
  real :: declin, solcon, xt
  real :: coszen(1,1), hrang(1,1), xlat(1,1), xlon(1,1)
  character(len=1024) :: outdir
  integer :: n
  names = [character(len=32) :: 'winter_branch', 'equinox_edge', 'southern_hemisphere', &
           'radt_offset', 'multi_day_wrap', 'polar_day']
  julian = [35.25, 79.9, 172.5, 295.5, 297.2, 172.25]
  xtime  = [0., 30., 120., 700., 4000., 0.]
  radt   = [0., 0., 0., 30., 10., 0.]
  gmt    = [6., 12., 0., 19., 19., 6.]
  lat    = [37.6, 51.5, -33.9, 37.6, 37.6, 89.9]
  lon    = [-96.7, 0.0, 151.2, -96.7, -96.7, 0.0]
  call get_command_argument(1, outdir)
  do n = 1, nreg
    xlat = lat(n); xlon = lon(n)
    call radconst(xtime(n), declin, solcon, julian(n), degrad, dpd)
    xt = xtime(n) + radt(n)*0.5
    call calc_coszen(1,1,1,1,1,1,1,1, julian(n), xt, gmt(n), declin, degrad, xlon, xlat, coszen, hrang)
    call esm_dump_open_file('solar', trim(outdir)//'/esm_dump_solar_'//trim(names(n))//'.json')
    call esm_dump_var('regime', trim(names(n)))
    call esm_dump_var('julian', julian(n)); call esm_dump_var('xtime', xtime(n)); call esm_dump_var('radt', radt(n))
    call esm_dump_var('gmt', gmt(n)); call esm_dump_var('xlat', lat(n)); call esm_dump_var('xlong', lon(n))
    call esm_dump_var('degrad', degrad); call esm_dump_var('dpd', dpd)
    call esm_dump_var('declin', declin); call esm_dump_var('solcon', solcon)
    call esm_dump_var('coszen', coszen(1,1)); call esm_dump_var('hrang', hrang(1,1))
    call esm_dump_close()
  end do
end program solar_driver
