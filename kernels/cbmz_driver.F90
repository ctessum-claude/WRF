! cbmz_driver: replay one dumped column (esm dump scheme 'cbmz_kpp', written by
! the chem/KPP/inc/cbmz_mosaic hooks) through the KPP CBM-Z mechanism and write
! the results as an esm dump.
!   cbmz_driver <input.flat> <output.json>     env ESM_DT overrides dtstepc
! For every dumped level it
!   (a) recomputes RCONST from TEMP, C_M, C_H2O and jv with
!       cbmz_mosaic_Update_RCONST, to check the dumped rate coefficients,
!   (b) evaluates the right-hand side cbmz_mosaic_Fun(var_in, fix, RCONST) --
!       the INSTANTANEOUS production-minus-loss rates, which WRF never dumps --
!       both with the dumped RCONST and with the recomputed one, and
!   (c) re-integrates the step with cbmz_mosaic_INTEGRATE, to check var_out.
! KPP is already double precision inside WRF (kind dp = real64), so this is not
! a precision upgrade: it is an independent replay of the same arithmetic, and
! (b) is the quantity an .esm reaction system must reproduce.
program cbmz_driver
  use esm_flat_input
  use module_esm_dump
  use cbmz_mosaic_Precision, only: dp
  use cbmz_mosaic_Parameters, only: NVAR, NFIX, NREACT, NSPEC, ind_CH3O2, ind_ETHP, ind_RO2, &
       ind_C2O3, ind_ANO2, ind_NAP, ind_ISOPP, ind_ISOPN, ind_ISOPO2, ind_XO2
  use cbmz_mosaic_UpdateRconstWRF, only: cbmz_mosaic_Update_RCONST
  use cbmz_mosaic_Integrator, only: cbmz_mosaic_INTEGRATE, cbmz_mosaic_Fun
  implicit none
  integer, parameter :: njv = 52
  character(len=1024) :: inpath, outpath, dtstr
  integer :: kts, kte, nk, k, stat, n, ierr
  real(dp) :: dtstepc, tstart, tend
  real(dp), allocatable :: var_in(:,:), var_out(:,:), fixd(:,:), rconst(:,:), jv(:,:)
  real(dp), allocatable :: temp(:), c_m(:), c_h2o(:)
  real(dp), allocatable :: rconst_re(:,:), vdot(:,:), vdot_re(:,:), var_re(:,:)
  real(dp) :: atol(NSPEC), rtol(NSPEC), irr(NREACT), y(NVAR), f(NFIX), rc(NREACT)
  integer :: icntrl(20)
  real(dp) :: rcntrl(20)
  ! Pj_* indices, in the order the generated interface assigns them
  integer :: p(njv)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  kts = flat_i0('kts'); kte = flat_i0('kte'); nk = kte - kts + 1
  dtstepc = flat_r0('dtstepc')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) dtstepc
  allocate(var_in(NVAR,nk), var_out(NVAR,nk), fixd(NFIX,nk), rconst(NREACT,nk), jv(njv,nk))
  allocate(temp(nk), c_m(nk), c_h2o(nk))
  allocate(rconst_re(NREACT,nk), vdot(NVAR,nk), vdot_re(NVAR,nk), var_re(NVAR,nk))
  call flat_r2('var_in', var_in); call flat_r2('var_out', var_out); call flat_r2('fix', fixd)
  call flat_r2('rconst', rconst); call flat_r2('jv', jv)
  call flat_r1('temp', temp); call flat_r1('c_m', c_m); call flat_r1('c_h2o', c_h2o)
  call flat_r1('atol', atol); call flat_r1('rtol', rtol)
  call flat_i1('icntrl', icntrl); call flat_r1('rcntrl', rcntrl)
  do n = 1, njv
    p(n) = n
  end do

  do k = 1, nk
    y = var_in(:,k); f = fixd(:,k)
    ! (a) rate coefficients from the dumped state
    call cbmz_mosaic_Update_RCONST( &
         y(ind_CH3O2), y(ind_ETHP), y(ind_RO2), y(ind_C2O3), y(ind_ANO2), &
         y(ind_NAP), y(ind_ISOPP), y(ind_ISOPN), y(ind_ISOPO2), y(ind_XO2), &
         jv(:,k), njv, rconst_re(:,k), &
         p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), &
         p(11), p(12), p(13), p(14), p(15), p(16), p(17), p(18), p(19), p(20), &
         p(21), p(22), p(23), p(24), p(25), p(26), p(27), p(28), p(29), p(30), &
         p(31), p(32), p(33), p(34), p(35), p(36), p(37), p(38), p(39), p(40), &
         p(41), p(42), p(43), p(44), p(45), p(46), p(47), p(48), p(49), p(50), &
         p(51), p(52), c_m(k), c_h2o(k), temp(k) )
    ! (b) instantaneous rates, with the dumped and with the recomputed RCONST
    call cbmz_mosaic_Fun(y, f, rconst(:,k), vdot(:,k))
    call cbmz_mosaic_Fun(y, f, rconst_re(:,k), vdot_re(:,k))
    ! (c) re-integrate the operator-split step
    rc = rconst(:,k); irr = 0._dp; tstart = 0._dp; tend = dtstepc
    call cbmz_mosaic_INTEGRATE(tstart, tend, f, y, rc, atol, rtol, irr, &
         ICNTRL_U=icntrl, RCNTRL_U=rcntrl, IERR_U=ierr)
    var_re(:,k) = y
  end do

  call esm_dump_open_file('cbmz_kpp', trim(outpath))
  call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('dtstepc', dtstepc); call esm_dump_var('dp_bytes', storage_size(dtstepc)/8)
  call esm_dump_var('rconst_replay', rconst_re)
  call esm_dump_var('vdot', vdot)
  call esm_dump_var('vdot_rconst_replay', vdot_re)
  call esm_dump_var('var_out_replay', var_re)
  call esm_dump_var('ierr', ierr)
  call esm_dump_close()
end program cbmz_driver
