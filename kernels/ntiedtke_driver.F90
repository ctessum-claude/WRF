! ntiedtke_driver: replay one dumped WRF row (esm dump scheme 'ntiedtke', written
! by phys/module_cu_ntiedtke.F) through cu_ntiedtke_run (real64 when built with
! -DDOUBLE_PRECISION) and write the outputs as an esm dump.
!   ntiedtke_driver <input.flat> <output.json>      env ESM_DT overrides delt
! The kernel reads the state in its own top-down ordering (the dump's pu, pv,
! pt, pqv, pqc, pqi, pqvf, ptf, poz, pzz, pomg, pap, paph, which
! cu_ntiedtke_pre_run built), updates it in place as x + (dx/dt)*delt and
! returns the surface convective precipitation over delt.  The driver then
! forms the WRF tendencies exactly as cu_ntiedtke_post_run does (bottom-up,
! theta tendency = dT/dt / exner).  The kernel is compiled from a copy of
! physics_mmm/cu_ntiedtke.F90 that includes ntiedtke_esm_inner.inc right after
! the cumastrn call (see the Makefile), so the cloud type, base and top and the
! cumastrn tendencies land in the output as well.
program ntiedtke_driver
  use ccpp_kind_types, only: kind_phys
  use esm_flat_input
  use module_esm_dump
  use cu_ntiedtke, only: cu_ntiedtke_init, cu_ntiedtke_run
  implicit none
  character(len=1024) :: inpath, outpath, dtstr
  character(len=512) :: errmsg
  integer :: its, ite, kts, kte, im, kx, kx1, stat, errflg, i, k, zz
  real(kind_phys) :: delt, cp, rd, rv, xlv, xls, xlf, grav
  integer, allocatable :: slimsk(:)
  real(kind_phys), allocatable, dimension(:) :: dx, hfx, qfx, zprecc
  real(kind_phys), allocatable, dimension(:,:) :: pu, pv, pt, pqv, pqc, pqi, pqvf, ptf, poz, pomg, pap, pzz, paph
  real(kind_phys), allocatable, dimension(:,:) :: pu0, pv0, pt0, pqv0, pqc0, pqi0, pi3d
  real(kind_phys), allocatable, dimension(:,:) :: rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite'); kts = flat_i0('kts'); kte = flat_i0('kte')
  im = flat_i0('im'); kx = flat_i0('kx'); kx1 = flat_i0('kx1')
  delt = flat_r0('delt')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) delt
  cp = flat_r0('cp'); rd = flat_r0('rd'); rv = flat_r0('rv')
  xlv = flat_r0('xlv'); xls = flat_r0('xls'); xlf = flat_r0('xlf'); grav = flat_r0('grav')

  allocate(slimsk(im), dx(im), hfx(im), qfx(im), zprecc(im))
  allocate(pu(im,kx), pv(im,kx), pt(im,kx), pqv(im,kx), pqc(im,kx), pqi(im,kx), pqvf(im,kx), ptf(im,kx), &
           poz(im,kx), pomg(im,kx), pap(im,kx), pzz(im,kx1), paph(im,kx1), pi3d(im,kx))
  allocate(rthcuten(im,kx), rqvcuten(im,kx), rqccuten(im,kx), rqicuten(im,kx), rucuten(im,kx), rvcuten(im,kx))
  call flat_i1('slimsk', slimsk)
  call flat_r1('dx', dx); call flat_r1('hfx', hfx); call flat_r1('qfx', qfx)
  call flat_r2('pu', pu); call flat_r2('pv', pv); call flat_r2('pt', pt)
  call flat_r2('pqv', pqv); call flat_r2('pqc', pqc); call flat_r2('pqi', pqi)
  call flat_r2('pqvf', pqvf); call flat_r2('ptf', ptf)
  call flat_r2('poz', poz); call flat_r2('pzz', pzz); call flat_r2('pomg', pomg)
  call flat_r2('pap', pap); call flat_r2('paph', paph)
  call flat_r2('pi3d', pi3d)
  pu0 = pu; pv0 = pv; pt0 = pt; pqv0 = pqv; pqc0 = pqc; pqi0 = pqi
  zprecc = 0

  call cu_ntiedtke_init(con_cp=cp, con_rd=rd, con_rv=rv, con_xlv=xlv, con_xls=xls, con_xlf=xlf, &
                        con_grav=grav, errmsg=errmsg, errflg=errflg)

  call esm_dump_open_file('ntiedtke', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite); call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('delt', delt); call esm_dump_var('kind_phys_bytes', storage_size(delt)/8)
  esm_dump_inner = .true.
  call cu_ntiedtke_run(pu=pu, pv=pv, pt=pt, pqv=pqv, pqc=pqc, pqi=pqi, pqvf=pqvf, ptf=ptf, &
       poz=poz, pzz=pzz, pomg=pomg, pap=pap, paph=paph, evap=qfx, hfx=hfx, zprecc=zprecc, &
       lndj=slimsk, lq=im, km=kx, km1=kx1, dt=delt, dx=dx, errmsg=errmsg, errflg=errflg)
  esm_dump_inner = .false.

  ! cu_ntiedtke_post_run: kernel level zz = kte-(k-kts) is WRF level k
  do k = kts, kte
    zz = kte - (k - kts)
    do i = 1, im
      rthcuten(i,k-kts+1) = (pt(i,zz) - pt0(i,zz)) / pi3d(i,k-kts+1) / delt
      rqvcuten(i,k-kts+1) = (pqv(i,zz) - pqv0(i,zz)) / delt
      rqccuten(i,k-kts+1) = (pqc(i,zz) - pqc0(i,zz)) / delt
      rqicuten(i,k-kts+1) = (pqi(i,zz) - pqi0(i,zz)) / delt
      rucuten(i,k-kts+1)  = (pu(i,zz) - pu0(i,zz)) / delt
      rvcuten(i,k-kts+1)  = (pv(i,zz) - pv0(i,zz)) / delt
    end do
  end do
  call esm_dump_var('pu_out', pu); call esm_dump_var('pv_out', pv); call esm_dump_var('pt_out', pt)
  call esm_dump_var('pqv_out', pqv); call esm_dump_var('pqc_out', pqc); call esm_dump_var('pqi_out', pqi)
  call esm_dump_var('zprecc', zprecc); call esm_dump_var('pratec', zprecc / delt)
  call esm_dump_var('rthcuten', rthcuten); call esm_dump_var('rqvcuten', rqvcuten)
  call esm_dump_var('rqccuten', rqccuten); call esm_dump_var('rqicuten', rqicuten)
  call esm_dump_var('rucuten', rucuten); call esm_dump_var('rvcuten', rvcuten)
  call esm_dump_var('errflg', errflg)
  call esm_dump_close()
end program ntiedtke_driver
