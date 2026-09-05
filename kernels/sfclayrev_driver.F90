! sfclayrev_driver: replay one dumped WRF tile through sf_sfclayrev_run
! (real64 when built with -DDOUBLE_PRECISION) and write the outputs as an
! esm dump with the same names as the WRF wrapper (<name>_out for inout).
!   sfclayrev_driver <input.flat> <output.json>
! The scheme is diagnostic (no time step), so no ESM_DT handling is needed.
program sfclayrev_driver
  use ccpp_kind_types, only: kind_phys
  use esm_flat_input
  use module_esm_dump
  use sf_sfclayrev, only: sf_sfclayrev_init, sf_sfclayrev_run
  implicit none
  character(len=1024) :: inpath, outpath
  character(len=512) :: errmsg
  integer :: its, ite, errflg, isftcflx, iz0tlnd
  real(kind_phys) :: cp, g, rovcp, r, xlv, svp1, svp2, svp3, svpt0, ep1, ep2, karman, p1000mb, tofd_factor
  logical :: isfflx, if_kim_tofd, shalwater_z0, scm_force_flux
  real(kind_phys), allocatable, dimension(:) :: ux, vx, t1d, qv1d, p1d, dz8w1d, psfcpa, mavail, pblh, tsk, &
       xland, lakemask, water_depth, dx, varf, chs, chs2, cqs2, cpm, rmol, znt, ust, zol, mol, regime, psim, &
       psih, fm, fh, hfx, qfx, qsfc, gz1oz0, wspd, br, flhc, flqc, qgh, ustm, lh, u10, v10, th2, t2, q2, &
       ck, cka, cd, cda

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite')
  cp = flat_r0('cp'); g = flat_r0('g'); rovcp = flat_r0('rovcp'); r = flat_r0('r'); xlv = flat_r0('xlv')
  svp1 = flat_r0('svp1'); svp2 = flat_r0('svp2'); svp3 = flat_r0('svp3'); svpt0 = flat_r0('svpt0')
  ep1 = flat_r0('ep1'); ep2 = flat_r0('ep2'); karman = flat_r0('karman'); p1000mb = flat_r0('p1000mb')
  tofd_factor = flat_r0('tofd_factor'); if_kim_tofd = flat_l0('if_kim_tofd')
  isfflx = flat_l0('isfflx'); shalwater_z0 = flat_l0('shalwater_z0'); scm_force_flux = flat_l0('scm_force_flux')
  isftcflx = flat_i0('isftcflx'); iz0tlnd = flat_i0('iz0tlnd')

  allocate(ux(its:ite), vx(its:ite), t1d(its:ite), qv1d(its:ite), p1d(its:ite), dz8w1d(its:ite), psfcpa(its:ite), &
       mavail(its:ite), pblh(its:ite), tsk(its:ite), xland(its:ite), lakemask(its:ite), water_depth(its:ite), &
       dx(its:ite), varf(its:ite), chs(its:ite), chs2(its:ite), cqs2(its:ite), cpm(its:ite), rmol(its:ite), &
       znt(its:ite), ust(its:ite), zol(its:ite), mol(its:ite), regime(its:ite), psim(its:ite), psih(its:ite), &
       fm(its:ite), fh(its:ite), hfx(its:ite), qfx(its:ite), qsfc(its:ite), gz1oz0(its:ite), wspd(its:ite), &
       br(its:ite), flhc(its:ite), flqc(its:ite), qgh(its:ite), ustm(its:ite), lh(its:ite), u10(its:ite), &
       v10(its:ite), th2(its:ite), t2(its:ite), q2(its:ite), ck(its:ite), cka(its:ite), cd(its:ite), cda(its:ite))
  call flat_r1('ux', ux); call flat_r1('vx', vx); call flat_r1('t1d', t1d); call flat_r1('qv1d', qv1d)
  call flat_r1('p1d', p1d); call flat_r1('dz8w1d', dz8w1d); call flat_r1('psfcpa', psfcpa); call flat_r1('mavail', mavail)
  call flat_r1('pblh', pblh); call flat_r1('tsk', tsk); call flat_r1('xland', xland); call flat_r1('lakemask', lakemask)
  call flat_r1('water_depth', water_depth); call flat_r1('dx', dx); call flat_r1('varf', varf)
  call flat_r1('chs', chs); call flat_r1('chs2', chs2); call flat_r1('cqs2', cqs2); call flat_r1('cpm', cpm)
  call flat_r1('rmol', rmol); call flat_r1('znt', znt); call flat_r1('ust', ust); call flat_r1('zol', zol)
  call flat_r1('mol', mol); call flat_r1('regime', regime); call flat_r1('psim', psim); call flat_r1('psih', psih)
  call flat_r1('fm', fm); call flat_r1('fh', fh); call flat_r1('hfx', hfx); call flat_r1('qfx', qfx)
  call flat_r1('qsfc', qsfc); call flat_r1('gz1oz0', gz1oz0); call flat_r1('wspd', wspd); call flat_r1('br', br)
  call flat_r1('flhc', flhc); call flat_r1('flqc', flqc); call flat_r1('qgh', qgh); call flat_r1('ustm', ustm)
  ck = 0; cka = 0; cd = 0; cda = 0

  call sf_sfclayrev_init(errmsg, errflg)
  call sf_sfclayrev_run(ux=ux,vx=vx,t1d=t1d,qv1d=qv1d,p1d=p1d,dz8w1d=dz8w1d,cp=cp,g=g,rovcp=rovcp,r=r,xlv=xlv, &
       psfcpa=psfcpa,chs=chs,chs2=chs2,cqs2=cqs2,cpm=cpm,pblh=pblh,rmol=rmol,znt=znt,ust=ust,mavail=mavail, &
       zol=zol,mol=mol,regime=regime,psim=psim,psih=psih,fm=fm,fh=fh,xland=xland,lakemask=lakemask, &
       hfx=hfx,qfx=qfx,tsk=tsk,u10=u10,varf=varf,if_kim_tofd=if_kim_tofd,tofd_factor=tofd_factor, &
       v10=v10,th2=th2,t2=t2,q2=q2,flhc=flhc,flqc=flqc,qgh=qgh,qsfc=qsfc,lh=lh,gz1oz0=gz1oz0,wspd=wspd,br=br, &
       isfflx=isfflx,dx=dx,svp1=svp1,svp2=svp2,svp3=svp3,svpt0=svpt0,ep1=ep1,ep2=ep2,karman=karman, &
       p1000mb=p1000mb,shalwater_z0=shalwater_z0,water_depth=water_depth,its=its,ite=ite,errmsg=errmsg,errflg=errflg, &
       isftcflx=isftcflx,iz0tlnd=iz0tlnd,scm_force_flux=scm_force_flux,ustm=ustm,ck=ck,cka=cka,cd=cd,cda=cda)

  call esm_dump_open_file('sfclayrev', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite); call esm_dump_var('kind_phys_bytes', storage_size(cp)/8)
  call esm_dump_var('chs_out', chs); call esm_dump_var('chs2_out', chs2); call esm_dump_var('cqs2_out', cqs2)
  call esm_dump_var('cpm_out', cpm); call esm_dump_var('rmol_out', rmol); call esm_dump_var('znt_out', znt)
  call esm_dump_var('ust_out', ust); call esm_dump_var('zol_out', zol); call esm_dump_var('mol_out', mol)
  call esm_dump_var('regime_out', regime); call esm_dump_var('psim_out', psim); call esm_dump_var('psih_out', psih)
  call esm_dump_var('fm_out', fm); call esm_dump_var('fh_out', fh); call esm_dump_var('hfx_out', hfx)
  call esm_dump_var('qfx_out', qfx); call esm_dump_var('qsfc_out', qsfc); call esm_dump_var('gz1oz0_out', gz1oz0)
  call esm_dump_var('wspd_out', wspd); call esm_dump_var('br_out', br); call esm_dump_var('flhc_out', flhc)
  call esm_dump_var('flqc_out', flqc); call esm_dump_var('qgh_out', qgh); call esm_dump_var('ustm_out', ustm)
  call esm_dump_var('lh', lh); call esm_dump_var('u10', u10); call esm_dump_var('v10', v10)
  call esm_dump_var('th2', th2); call esm_dump_var('t2', t2); call esm_dump_var('q2', q2)
  call esm_dump_var('ck', ck); call esm_dump_var('cka', cka); call esm_dump_var('cd', cd); call esm_dump_var('cda', cda)
  call esm_dump_var('errflg', errflg)
  call esm_dump_close()
end program sfclayrev_driver
