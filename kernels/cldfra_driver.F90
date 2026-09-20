! cldfra_driver: replay one dumped WRF tile (one j-row) through cal_cldfra1, the
! icloud = 1 grid-scale cloud-fraction diagnostic of phys/module_radiation_driver.F,
! and write CLDFRA and the branch flag as an esm dump.
!   cldfra_driver <input.flat> <output.json>
! cal_cldfra1 is EXTRACTED VERBATIM from module_radiation_driver.F at build time
! (Makefile rule cldfra_kernels.F90), so the driver cannot drift from the model
! code. It is a WRF-native default-REAL routine, promoted by -fdefault-real-8 in
! the real64 build; the scheme is a pure diagnostic (no time step), so no dt
! sensitivity applies. The saturation difference RHGRID*QVS_WEIGHT - QV that the
! Randall (1994) branch raises to the power GAMMA is a cancellation, which is why
! the real32 in-model values are not an adequate reference near saturation
! (EqWeFiC PLAN.md section 3.1): the real64 replay is.
program cldfra_driver
  use esm_flat_input
  use module_esm_dump
  use cldfra_kernels, only: cal_cldfra1
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: its, ite, kts, kte, mp_physics, k
  logical :: f_qv, f_qc, f_qi, f_qs
  real, allocatable, dimension(:,:,:) :: cldfra, qv, qc, qi, qs, t_phy, p_phy, f_ice_phy, f_rain_phy
  integer, allocatable, dimension(:,:,:) :: cldfra1_flag

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite'); kts = flat_i0('kts'); kte = flat_i0('kte')
  mp_physics = flat_i0('mp_physics')
  f_qv = flat_l0('f_qv'); f_qc = flat_l0('f_qc'); f_qi = flat_l0('f_qi'); f_qs = flat_l0('f_qs')
  allocate(cldfra(its:ite,kts:kte,1), qv(its:ite,kts:kte,1), qc(its:ite,kts:kte,1), &
       qi(its:ite,kts:kte,1), qs(its:ite,kts:kte,1), t_phy(its:ite,kts:kte,1), &
       p_phy(its:ite,kts:kte,1), f_ice_phy(its:ite,kts:kte,1), f_rain_phy(its:ite,kts:kte,1), &
       cldfra1_flag(its:ite,kts:kte,1))
  call flat_r2('qv3d', qv(:,:,1)); call flat_r2('qc3d', qc(:,:,1)); call flat_r2('qi3d', qi(:,:,1))
  call flat_r2('qs3d', qs(:,:,1)); call flat_r2('t3d', t_phy(:,:,1)); call flat_r2('p3d', p_phy(:,:,1))
  ! WSM6 (mp_physics = 6) leaves F_ICE_PHY/F_RAIN_PHY out of the actual argument
  ! list; they are OPTIONAL and only the mp = 5 branch reads them.
  f_ice_phy = 0.; f_rain_phy = 0.; cldfra = 0.; cldfra1_flag = 0

  call esm_dump_open_file('cldfra', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite)
  call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('mp_physics', mp_physics)
  call esm_dump_var('real_bytes', storage_size(cldfra(its,kts,1))/8)
  call cal_cldfra1(cldfra, qv, qc, qi, qs, f_qv, f_qc, f_qi, f_qs, t_phy, p_phy, &
       f_ice_phy, f_rain_phy, mp_physics, cldfra1_flag, &
       its,ite, 1,1, kts,kte,  its,ite, 1,1, kts,kte,  its,ite, 1,1, kts,kte)
  call esm_dump_var('cldfra3d_out', cldfra(:,:,1))
  call esm_dump_var('cldfra1_flag_out', real(cldfra1_flag(:,:,1)))
  call esm_dump_close()
end program cldfra_driver
