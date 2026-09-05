! ysu_driver: replay one dumped WRF tile through bl_ysu_run (real64 when built
! with -DDOUBLE_PRECISION) and write the outputs as an esm dump.
!   ysu_driver <input.flat> <output.json>      env ESM_DT overrides the time step
! The instantaneous tendency of the implicit scheme is approximated by running
! with a small ESM_DT (e.g. 0.01 s); see EqWeFiC PLAN.md section 6.
program ysu_driver
  use, intrinsic :: iso_fortran_env, only: real64
  use ccpp_kind_types, only: kind_phys
  use esm_flat_input
  use module_esm_dump
  use bl_ysu, only: bl_ysu_run
  implicit none
  character(len=1024) :: inpath, outpath, dtstr
  integer :: its, ite, kts, kte, kme, nmix, stat, errflg
  real(kind_phys) :: dt, cp, g, rovcp, rd, rovg, xlv, rv, ep1, ep2, karman
  logical :: f_qc, f_qi, topdown, flag_bep
  real(kind_phys), allocatable, dimension(:,:) :: ux, vx, tx, qvx, qcx, qix, p2d, pi2d, dz8w2d, rthraten, &
       utnp, vtnp, ttnp, qvtnp, qctnp, qitnp, exch_hx, exch_mx, p2di, &
       a_u, a_v, a_t, a_q, a_e, b_u, b_v, b_t, b_q, b_e, sfk, vlk, dlu, dlg
  real(kind_phys), allocatable, dimension(:,:,:) :: qmix, qmixtnp
  real(kind_phys), allocatable, dimension(:) :: psfcpa, znt, ust, hpbl, dusfc, dvsfc, dtsfc, dqsfc, psim, psih, &
       xland, hfx, qfx, wspd, br, wstar, delta, u10, v10, uox, vox, ctopo, ctopo2, frcurb
  integer, allocatable :: kpbl1d(:)
  character(len=512) :: errmsg

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite'); kts = flat_i0('kts'); kte = flat_i0('kte'); kme = flat_i0('kme')
  nmix = 0
  dt = flat_r0('dt')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) dt
  cp = flat_r0('cp'); g = flat_r0('g'); rovcp = flat_r0('rovcp'); rd = flat_r0('rd'); rovg = flat_r0('rovg')
  xlv = flat_r0('xlv'); rv = flat_r0('rv'); ep1 = flat_r0('ep1'); ep2 = flat_r0('ep2'); karman = flat_r0('karman')
  f_qc = flat_l0('f_qc'); f_qi = flat_l0('f_qi'); topdown = flat_l0('ysu_topdown_pblmix'); flag_bep = flat_l0('flag_bep')

  allocate(ux(its:ite,kts:kte), vx(its:ite,kts:kte), tx(its:ite,kts:kte), qvx(its:ite,kts:kte), &
           qcx(its:ite,kts:kte), qix(its:ite,kts:kte), p2d(its:ite,kts:kte), pi2d(its:ite,kts:kte), &
           dz8w2d(its:ite,kts:kte), rthraten(its:ite,kts:kte), p2di(its:ite,kts:kte+1), &
           utnp(its:ite,kts:kte), vtnp(its:ite,kts:kte), ttnp(its:ite,kts:kte), qvtnp(its:ite,kts:kte), &
           qctnp(its:ite,kts:kte), qitnp(its:ite,kts:kte), exch_hx(its:ite,kts:kte), exch_mx(its:ite,kts:kte), &
           a_u(its:ite,kts:kte), a_v(its:ite,kts:kte), a_t(its:ite,kts:kte), a_q(its:ite,kts:kte), a_e(its:ite,kts:kte), &
           b_u(its:ite,kts:kte), b_v(its:ite,kts:kte), b_t(its:ite,kts:kte), b_q(its:ite,kts:kte), b_e(its:ite,kts:kte), &
           sfk(its:ite,kts:kte), vlk(its:ite,kts:kte), dlu(its:ite,kts:kte), dlg(its:ite,kts:kte), &
           qmix(its:ite,kts:kte,nmix), qmixtnp(its:ite,kts:kte,nmix))
  allocate(psfcpa(its:ite), znt(its:ite), ust(its:ite), hpbl(its:ite), dusfc(its:ite), dvsfc(its:ite), dtsfc(its:ite), &
           dqsfc(its:ite), psim(its:ite), psih(its:ite), xland(its:ite), hfx(its:ite), qfx(its:ite), wspd(its:ite), &
           br(its:ite), wstar(its:ite), delta(its:ite), u10(its:ite), v10(its:ite), uox(its:ite), vox(its:ite), &
           ctopo(its:ite), ctopo2(its:ite), frcurb(its:ite), kpbl1d(its:ite))
  call flat_r2('ux', ux); call flat_r2('vx', vx); call flat_r2('tx', tx); call flat_r2('qvx', qvx)
  call flat_r2('qcx', qcx); call flat_r2('qix', qix); call flat_r2('p2d', p2d); call flat_r2('p2di', p2di)
  call flat_r2('pi2d', pi2d); call flat_r2('dz8w2d', dz8w2d); call flat_r2('rthraten', rthraten)
  call flat_r1('psfcpa', psfcpa); call flat_r1('znt', znt); call flat_r1('ust', ust); call flat_r1('psim', psim)
  call flat_r1('psih', psih); call flat_r1('xland', xland); call flat_r1('hfx', hfx); call flat_r1('qfx', qfx)
  call flat_r1('wspd', wspd); call flat_r1('br', br); call flat_r1('u10_in', u10); call flat_r1('v10_in', v10)
  call flat_r1('uox', uox); call flat_r1('vox', vox); call flat_r1('ctopo', ctopo); call flat_r1('ctopo2', ctopo2)
  call flat_r1('frcurb', frcurb)
  a_u = 0; a_v = 0; a_t = 0; a_q = 0; a_e = 0; b_u = 0; b_v = 0; b_t = 0; b_q = 0; b_e = 0
  sfk = 0; vlk = 0; dlu = 0; dlg = 0

  call bl_ysu_run(ux=ux,vx=vx,tx=tx,qvx=qvx,qcx=qcx,qix=qix,nmix=nmix,qmix=qmix,p2d=p2d,p2di=p2di,pi2d=pi2d, &
       f_qc=f_qc,f_qi=f_qi,utnp=utnp,vtnp=vtnp,ttnp=ttnp,qvtnp=qvtnp,qctnp=qctnp,qitnp=qitnp,qmixtnp=qmixtnp, &
       cp=cp,g=g,rovcp=rovcp,rd=rd,rovg=rovg,ep1=ep1,ep2=ep2,karman=karman,xlv=xlv,rv=rv, &
       dz8w2d=dz8w2d,psfcpa=psfcpa,znt=znt,ust=ust,hpbl=hpbl,dusfc=dusfc,dvsfc=dvsfc,dtsfc=dtsfc,dqsfc=dqsfc, &
       psim=psim,psih=psih,xland=xland,hfx=hfx,qfx=qfx,wspd=wspd,br=br,dt=dt,kpbl1d=kpbl1d, &
       exch_hx=exch_hx,exch_mx=exch_mx,wstar=wstar,delta=delta,u10=u10,v10=v10,uox=uox,vox=vox, &
       rthraten=rthraten,ysu_topdown_pblmix=topdown,ctopo=ctopo,ctopo2=ctopo2, &
       a_u=a_u,a_v=a_v,a_t=a_t,a_q=a_q,a_e=a_e,b_u=b_u,b_v=b_v,b_t=b_t,b_q=b_q,b_e=b_e, &
       sfk=sfk,vlk=vlk,dlu=dlu,dlg=dlg,frcurb=frcurb,flag_bep=flag_bep,its=its,ite=ite,kte=kte,kme=kme, &
       errmsg=errmsg,errflg=errflg)

  call esm_dump_open_file('ysu', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite); call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('dt', dt); call esm_dump_var('kind_phys_bytes', storage_size(dt)/8)
  call esm_dump_var('utnp', utnp); call esm_dump_var('vtnp', vtnp); call esm_dump_var('ttnp', ttnp)
  call esm_dump_var('qvtnp', qvtnp); call esm_dump_var('qctnp', qctnp); call esm_dump_var('qitnp', qitnp)
  call esm_dump_var('exch_hx', exch_hx); call esm_dump_var('exch_mx', exch_mx)
  call esm_dump_var('hpbl', hpbl); call esm_dump_var('kpbl1d', kpbl1d); call esm_dump_var('wstar', wstar)
  call esm_dump_var('delta', delta); call esm_dump_var('u10', u10); call esm_dump_var('v10', v10)
  call esm_dump_var('dusfc', dusfc); call esm_dump_var('dvsfc', dvsfc); call esm_dump_var('dtsfc', dtsfc); call esm_dump_var('dqsfc', dqsfc)
  call esm_dump_var('errflg', errflg)
  call esm_dump_close()
end program ysu_driver
