! slab_driver: replay one dumped WRF tile (one j-row) through SLAB
! (module_sf_slab.F compiled with -fdefault-real-8 for the real64 build) and
! write the outputs (<name>_out) and SLAB1D's diagnostics as an esm dump.
!   slab_driver <input.flat> <output.json>      env ESM_DT overrides deltsm
! SLAB is an in-place step (TSK, TSLB advance by deltsm); the surface tendency
! dthgdt is dt-independent only when the soil sub-stepping nsoil == 1, so use
! a small ESM_DT for instantaneous rates (EqWeFiC PLAN.md section 6).
program slab_driver
  use esm_flat_input
  use module_esm_dump
  use module_sf_slab, only: slab
  implicit none
  character(len=1024) :: inpath, outpath, dtstr
  integer :: its, ite, kts, kte, nsl, ifsnow, stat
  logical :: radiation
  real :: deltsm, rovcp, xlv, dtmin, svp1, svp2, svp3, svpt0, ep2, karman, eomeg, stbolt, p1000mb
  real, allocatable, dimension(:) :: zs, dzs
  real, allocatable, dimension(:,:,:) :: t3d, qv3d, p3d, tslb
  real, allocatable, dimension(:,:) :: flhc, flqc, psfc, xland, tmn, hfx, qfx, lh, tsk, qsfc, chklowq, gsw, glw, &
       capg, thc, snowc, emiss, mavail

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite'); kts = flat_i0('kts'); kte = flat_i0('kte')
  nsl = flat_i0('num_soil_layers'); radiation = flat_l0('radiation'); ifsnow = flat_i0('ifsnow')
  deltsm = flat_r0('deltsm')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) deltsm
  rovcp = flat_r0('rovcp'); xlv = flat_r0('xlv'); dtmin = flat_r0('dtmin')
  svp1 = flat_r0('svp1'); svp2 = flat_r0('svp2'); svp3 = flat_r0('svp3'); svpt0 = flat_r0('svpt0')
  ep2 = flat_r0('ep2'); karman = flat_r0('karman'); eomeg = flat_r0('eomeg'); stbolt = flat_r0('stbolt')
  p1000mb = flat_r0('p1000mb')
  allocate(zs(nsl), dzs(nsl)); call flat_r1('zs', zs); call flat_r1('dzs', dzs)
  allocate(t3d(its:ite,1:1,1), qv3d(its:ite,1:1,1), p3d(its:ite,1:1,1), tslb(its:ite,nsl,1))
  allocate(flhc(its:ite,1), flqc(its:ite,1), psfc(its:ite,1), xland(its:ite,1), tmn(its:ite,1), hfx(its:ite,1), &
       qfx(its:ite,1), lh(its:ite,1), tsk(its:ite,1), qsfc(its:ite,1), chklowq(its:ite,1), gsw(its:ite,1), &
       glw(its:ite,1), capg(its:ite,1), thc(its:ite,1), snowc(its:ite,1), emiss(its:ite,1), mavail(its:ite,1))
  call flat_r1('t1d', t3d(:,1,1)); call flat_r1('qv1d', qv3d(:,1,1)); call flat_r1('p1d', p3d(:,1,1))
  call flat_r2('tslb', tslb(:,:,1))
  call flat_r1('flhc', flhc(:,1)); call flat_r1('flqc', flqc(:,1)); call flat_r1('psfc', psfc(:,1))
  call flat_r1('xland', xland(:,1)); call flat_r1('tmn', tmn(:,1)); call flat_r1('hfx', hfx(:,1))
  call flat_r1('qfx', qfx(:,1)); call flat_r1('lh', lh(:,1)); call flat_r1('tsk', tsk(:,1)); call flat_r1('qsfc', qsfc(:,1))
  call flat_r1('chklowq', chklowq(:,1)); call flat_r1('gsw', gsw(:,1)); call flat_r1('glw', glw(:,1))
  call flat_r1('capg', capg(:,1)); call flat_r1('thc', thc(:,1)); call flat_r1('snowc', snowc(:,1))
  call flat_r1('emiss', emiss(:,1)); call flat_r1('mavail', mavail(:,1))

  call esm_dump_open_file('slab', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite); call esm_dump_var('num_soil_layers', nsl)
  call esm_dump_var('deltsm', deltsm); call esm_dump_var('real_bytes', storage_size(deltsm)/8)
  esm_dump_inner = .true.
  call slab(t3d,qv3d,p3d,flhc,flqc,psfc,xland,tmn,hfx,qfx,lh,tsk,qsfc,chklowq,gsw,glw,capg,thc,snowc,emiss,mavail, &
       deltsm,rovcp,xlv,dtmin,ifsnow,svp1,svp2,svp3,svpt0,ep2,karman,eomeg,stbolt,tslb,zs,dzs,nsl,radiation,p1000mb, &
       its,ite,1,1,1,1, its,ite,1,1,1,1, its,ite,1,1,kts,kte)
  esm_dump_inner = .false.
  call esm_dump_var('hfx_out', hfx(:,1)); call esm_dump_var('qfx_out', qfx(:,1)); call esm_dump_var('lh_out', lh(:,1))
  call esm_dump_var('tsk_out', tsk(:,1)); call esm_dump_var('qsfc_out', qsfc(:,1)); call esm_dump_var('chklowq_out', chklowq(:,1))
  call esm_dump_var('capg_out', capg(:,1)); call esm_dump_var('flhc_out', flhc(:,1)); call esm_dump_var('flqc_out', flqc(:,1))
  call esm_dump_var('tslb_out', tslb(:,:,1))
  call esm_dump_close()
end program slab_driver
