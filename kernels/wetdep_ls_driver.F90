! wetdep_ls_driver: replay one dumped column (esm dump scheme 'wetdep_ls',
! written by the hook around the wetdep_ls call in chem/chem_driver.F,
! wetscav_onoff < 0) through chem/module_wetdep_ls.F, compiled unchanged,
! and write the scavenged column as an esm dump.
!   wetdep_ls_driver <input.flat> <output.json>
! Built real64 by default (-fdefault-real-8 promotes module_wetdep_ls and
! module_state_description), PREC= for a real32 build, so the two can be
! compared.  The flat must carry p_sulf (see below).  The column is i = j = 1; moist carries only qc (wetdep_ls reads
! no other moisture species).
program wetdep_ls_driver
  use esm_flat_input
  use module_esm_dump
  use module_state_description, only: p_qc, p_sulf
  use module_wetdep_ls, only: wetdep_ls
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: kts, kte, num_chem, numgas, chem_opt, nmoist
  real :: dt
  real, allocatable :: chem(:,:,:,:), moist(:,:,:,:), rho(:,:,:), dz8w(:,:,:), vvel(:,:,:)
  real, allocatable :: chem_in(:,:), col(:)
  real :: rain(1,1)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  kts = flat_i0('kts'); kte = flat_i0('kte')
  num_chem = flat_i0('num_chem'); numgas = flat_i0('numgas'); chem_opt = flat_i0('chem_opt')
  dt = flat_r0('dt'); rain(1,1) = flat_r0('rain')
  ! WRF sets the species pointers at run time (set_scalar_indices_from_config);
  ! standalone they keep their default of 1, which would make wetdep_ls skip
  ! sulf.  The flat carries the run's p_sulf (chem_opt 170: 3, from the
  ! registry package order); every other pointer wetdep_ls reads (bc/oc/dms,
  ! sea salt) is absent from chem_opt 170 and so is 1 in the model as well.
  p_sulf = flat_i0('p_sulf')
  nmoist = max(p_qc, 2)
  allocate(chem(1,kts:kte,1,num_chem), moist(1,kts:kte,1,nmoist))
  allocate(rho(1,kts:kte,1), dz8w(1,kts:kte,1), vvel(1,kts:kte,1))
  allocate(chem_in(kts:kte,num_chem), col(kts:kte))
  moist = 0.0
  call flat_r1('qc', col);   moist(1,:,1,p_qc) = col
  call flat_r1('rho', col);  rho(1,:,1) = col
  call flat_r1('dz8w', col); dz8w(1,:,1) = col
  call flat_r1('vvel', col); vvel(1,:,1) = col
  call flat_r2('chem_in', chem_in)
  chem(1,:,1,:) = chem_in

  call wetdep_ls(dt, chem, rain, moist, rho, nmoist,              &
                 num_chem, numgas, dz8w, vvel, chem_opt,          &
                 1,1, 1,1, kts,kte,                               &
                 1,1, 1,1, kts,kte,                               &
                 1,1, 1,1, kts,kte)

  call esm_dump_open_file('wetdep_ls', trim(outpath))
  call esm_dump_var('real_bytes', storage_size(dt)/8)
  call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('num_chem', num_chem); call esm_dump_var('numgas', numgas)
  call esm_dump_var('chem_out', chem(1,:,1,:))
  call esm_dump_close()
end program wetdep_ls_driver
