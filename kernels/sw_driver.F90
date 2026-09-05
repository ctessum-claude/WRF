! sw_driver: replay one dumped WRF tile (one j-row) through SWRAD/SWPARA
! (module_ra_sw.F compiled with -fdefault-real-8 for the real64 build) and
! write gsw and the theta tendency, plus SWPARA's column internals for the
! first column (esm_dump_inner), as an esm dump.
!   sw_driver <input.flat> <output.json>
! Dudhia SW is diagnostic (no time step).
program sw_driver
  use, intrinsic :: iso_fortran_env, only: real64
  use esm_flat_input
  use module_esm_dump
  use module_ra_sw, only: swrad, swinit
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: its, ite, kts, kte, kms, kme, julday, ghg_input, icloud, i, k
  real :: dt, gmt, r, cp, g, xtime, declin, solcon, radfrq, degrad, cssca, julian
  logical :: warm_rain, f_qv, f_qc, f_qr, f_qi, f_qs, f_qg
  real, allocatable, dimension(:,:,:) :: rthraten, rho_phy, t3d, qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, p3d, pi3d, dz8w
  real, allocatable, dimension(:,:) :: gsw, xlat, xlong, albedo, coszen, obscur, tmp2

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite'); kts = flat_i0('kts'); kte = flat_i0('kte')
  kms = flat_i0('kms'); kme = flat_i0('kme')
  dt = flat_r0('dt'); gmt = flat_r0('gmt'); r = flat_r0('r'); cp = flat_r0('cp'); g = flat_r0('g')
  julday = flat_i0('julday'); ghg_input = flat_i0('ghg_input'); xtime = flat_r0('xtime'); declin = flat_r0('declin')
  solcon = flat_r0('solcon'); radfrq = flat_r0('radfrq'); icloud = flat_i0('icloud'); degrad = flat_r0('degrad')
  warm_rain = flat_l0('warm_rain'); cssca = flat_r0('cssca')
  julian = 0.; if (flat_has('julian')) julian = flat_r0('julian')
  f_qv = flat_has('f_qv') .and. flat_has('qv3d'); if (f_qv) f_qv = flat_l0('f_qv')
  f_qc = flat_has('f_qc') .and. flat_has('qc3d'); if (f_qc) f_qc = flat_l0('f_qc')
  f_qr = flat_has('f_qr') .and. flat_has('qr3d'); if (f_qr) f_qr = flat_l0('f_qr')
  f_qi = flat_has('f_qi') .and. flat_has('qi3d'); if (f_qi) f_qi = flat_l0('f_qi')
  f_qs = flat_has('f_qs') .and. flat_has('qs3d'); if (f_qs) f_qs = flat_l0('f_qs')
  f_qg = flat_has('f_qg') .and. flat_has('qg3d'); if (f_qg) f_qg = flat_l0('f_qg')

  allocate(rthraten(its:ite,kms:kme,1), rho_phy(its:ite,kms:kme,1), t3d(its:ite,kms:kme,1), qv3d(its:ite,kms:kme,1), &
           qc3d(its:ite,kms:kme,1), qr3d(its:ite,kms:kme,1), qi3d(its:ite,kms:kme,1), qs3d(its:ite,kms:kme,1), &
           qg3d(its:ite,kms:kme,1), p3d(its:ite,kms:kme,1), pi3d(its:ite,kms:kme,1), dz8w(its:ite,kms:kme,1))
  allocate(gsw(its:ite,1), xlat(its:ite,1), xlong(its:ite,1), albedo(its:ite,1), coszen(its:ite,1), obscur(its:ite,1), &
           tmp2(its:ite,kts:kte))
  rthraten = 0; rho_phy = 0; t3d = 0; qv3d = 0; qc3d = 0; qr3d = 0; qi3d = 0; qs3d = 0; qg3d = 0; p3d = 0; pi3d = 1; dz8w = 0
  call flat_r2('t3d', tmp2); t3d(:,kts:kte,1) = tmp2
  call flat_r2('p3d', tmp2); p3d(:,kts:kte,1) = tmp2
  call flat_r2('pi3d', tmp2); pi3d(:,kts:kte,1) = tmp2
  call flat_r2('rho_phy', tmp2); rho_phy(:,kts:kte,1) = tmp2
  call flat_r2('dz8w', tmp2); dz8w(:,kts:kte,1) = tmp2
  if (flat_has('qv3d')) then; call flat_r2('qv3d', tmp2); qv3d(:,kts:kte,1) = tmp2; end if
  if (flat_has('qc3d')) then; call flat_r2('qc3d', tmp2); qc3d(:,kts:kte,1) = tmp2; end if
  if (flat_has('qr3d')) then; call flat_r2('qr3d', tmp2); qr3d(:,kts:kte,1) = tmp2; end if
  if (flat_has('qi3d')) then; call flat_r2('qi3d', tmp2); qi3d(:,kts:kte,1) = tmp2; end if
  if (flat_has('qs3d')) then; call flat_r2('qs3d', tmp2); qs3d(:,kts:kte,1) = tmp2; end if
  if (flat_has('qg3d')) then; call flat_r2('qg3d', tmp2); qg3d(:,kts:kte,1) = tmp2; end if
  call flat_r1('xlat', xlat(:,1)); call flat_r1('xlong', xlong(:,1)); call flat_r1('albedo', albedo(:,1))
  call flat_r1('coszen', coszen(:,1)); call flat_r1('obscur', obscur(:,1))

  ! swinit stores cssca = swrad_scat*1e-5; the index arguments are unused.
  call swinit(cssca*1.e5, .true., its,ite,1,1,kts,kte, its,ite,1,1,kms,kme, its,ite,1,1,kts,kte)

  call esm_dump_open_file('sw', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite); call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('real_bytes', storage_size(dt)/8)
  ! SWRAD sets esm_dump_inner for column its because a dump file is open.
  call swrad(dt=dt,rthraten=rthraten,gsw=gsw,xlat=xlat,xlong=xlong,albedo=albedo,rho_phy=rho_phy,t3d=t3d, &
       qv3d=qv3d,qc3d=qc3d,qr3d=qr3d,qi3d=qi3d,qs3d=qs3d,qg3d=qg3d,p3d=p3d,pi3d=pi3d,dz8w=dz8w,gmt=gmt, &
       r=r,cp=cp,g=g,julday=julday,ghg_input=ghg_input,xtime=xtime,declin=declin,solcon=solcon, &
       f_qv=f_qv,f_qc=f_qc,f_qr=f_qr,f_qi=f_qi,f_qs=f_qs,f_qg=f_qg,radfrq=radfrq,icloud=icloud,degrad=degrad, &
       warm_rain=warm_rain,ids=its,ide=ite,jds=1,jde=1,kds=kts,kde=kte,ims=its,ime=ite,jms=1,jme=1,kms=kms,kme=kme, &
       its=its,ite=ite,jts=1,jte=1,kts=kts,kte=kte,coszen=coszen,julian=julian,obscur=obscur)
  call esm_dump_var('gsw', gsw(:,1)); call esm_dump_var('rthraten_sw', rthraten(:,kts:kte,1))
  call esm_dump_close()
end program sw_driver
