! wsm6_driver: replay one dumped WRF tile through mp_wsm6_run (real64 when
! built with -DDOUBLE_PRECISION; the kernel is compiled with -DESM_DUMP so its
! process-rate dumps land in the output file) and write the outputs as an esm
! dump with the wrapper's names (<name>_out for inout).
!   wsm6_driver <input.flat> <output.json>      env ESM_DT overrides delt
! WSM6 is an in-place step: the state change over delt divided by delt is a
! dt-dependent tendency; run with a small ESM_DT to approximate the
! instantaneous rates (EqWeFiC PLAN.md section 6).
program wsm6_driver
  use ccpp_kind_types, only: kind_phys
  use esm_flat_input
  use module_esm_dump
  use mp_wsm6, only: mp_wsm6_init, mp_wsm6_run
  implicit none
  character(len=1024) :: inpath, outpath, dtstr
  character(len=512) :: errmsg
  integer :: its, ite, kts, kte, stat, errflg, hail_opt
  real(kind_phys) :: delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, xlv0, xlf0, den0, denr, cliq, cice, psat
  real(kind_phys), parameter :: dens = 100.0_kind_phys   ! rhosnow, as in module_physics_init
  logical :: has_snow, has_graupel
  real(kind_phys), allocatable, dimension(:,:) :: t, q, qc, qi, qr, qs, qg, den, p, delz
  real(kind_phys), allocatable, dimension(:) :: rain, rainncv, sr, snow, snowncv, graupel, graupelncv

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  its = flat_i0('its'); ite = flat_i0('ite'); kts = flat_i0('kts'); kte = flat_i0('kte')
  delt = flat_r0('delt')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) delt
  g = flat_r0('g'); cpd = flat_r0('cpd'); cpv = flat_r0('cpv'); rd = flat_r0('rd'); rv = flat_r0('rv')
  t0c = flat_r0('t0c'); ep1 = flat_r0('ep1'); ep2 = flat_r0('ep2'); qmin = flat_r0('qmin')
  xls = flat_r0('xls'); xlv0 = flat_r0('xlv0'); xlf0 = flat_r0('xlf0'); den0 = flat_r0('den0'); denr = flat_r0('denr')
  cliq = flat_r0('cliq'); cice = flat_r0('cice'); psat = flat_r0('psat')
  has_snow = flat_l0('has_snow'); has_graupel = flat_l0('has_graupel')
  hail_opt = 0
  if (flat_has('hail_opt')) hail_opt = flat_i0('hail_opt')

  allocate(t(its:ite,kts:kte), q(its:ite,kts:kte), qc(its:ite,kts:kte), qi(its:ite,kts:kte), qr(its:ite,kts:kte), &
           qs(its:ite,kts:kte), qg(its:ite,kts:kte), den(its:ite,kts:kte), p(its:ite,kts:kte), delz(its:ite,kts:kte))
  allocate(rain(its:ite), rainncv(its:ite), sr(its:ite), snow(its:ite), snowncv(its:ite), graupel(its:ite), graupelncv(its:ite))
  call flat_r2('t', t); call flat_r2('q', q); call flat_r2('qc', qc); call flat_r2('qi', qi); call flat_r2('qr', qr)
  call flat_r2('qs', qs); call flat_r2('qg', qg); call flat_r2('den', den); call flat_r2('p', p); call flat_r2('delz', delz)
  call flat_r1('rain', rain)
  snow = 0; graupel = 0
  if (has_snow) call flat_r1('snow', snow)
  if (has_graupel) call flat_r1('graupel', graupel)
  rainncv = 0; snowncv = 0; graupelncv = 0; sr = 0

  call mp_wsm6_init(den0, denr, dens, cliq, cpv, hail_opt, errmsg, errflg)

  call esm_dump_open_file('wsm6', trim(outpath))
  call esm_dump_var('its', its); call esm_dump_var('ite', ite); call esm_dump_var('kts', kts); call esm_dump_var('kte', kte)
  call esm_dump_var('delt', delt); call esm_dump_var('kind_phys_bytes', storage_size(delt)/8)
  esm_dump_inner = .true.
  call mp_wsm6_run(t=t,q=q,qc=qc,qi=qi,qr=qr,qs=qs,qg=qg,den=den,p=p,delz=delz,delt=delt,g=g,cpd=cpd,cpv=cpv, &
       rd=rd,rv=rv,t0c=t0c,ep1=ep1,ep2=ep2,qmin=qmin,xls=xls,xlv0=xlv0,xlf0=xlf0,den0=den0,denr=denr, &
       cliq=cliq,cice=cice,psat=psat,rain=rain,rainncv=rainncv,sr=sr,snow=snow,snowncv=snowncv, &
       graupel=graupel,graupelncv=graupelncv,its=its,ite=ite,kts=kts,kte=kte,errmsg=errmsg,errflg=errflg)
  esm_dump_inner = .false.
  call esm_dump_var('t_out', t); call esm_dump_var('q_out', q); call esm_dump_var('qc_out', qc)
  call esm_dump_var('qi_out', qi); call esm_dump_var('qr_out', qr); call esm_dump_var('qs_out', qs)
  call esm_dump_var('qg_out', qg); call esm_dump_var('rain_out', rain); call esm_dump_var('rainncv', rainncv)
  call esm_dump_var('sr', sr); call esm_dump_var('snow_out', snow); call esm_dump_var('snowncv', snowncv)
  call esm_dump_var('graupel_out', graupel); call esm_dump_var('graupelncv', graupelncv)
  call esm_dump_var('errflg', errflg)
  call esm_dump_close()
end program wsm6_driver
