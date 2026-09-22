! noah_driver: replay one dumped Noah point (esm dump scheme 'noah', written by
! phys/module_sf_noahdrv.F around its SFLX call for point (its, jts)) through
! SFLX (module_sf_noahlsm.F compiled with -fdefault-real-8 for the real64
! build) and write the outputs as an esm dump.
!   noah_driver <input.flat> <output.json>      env ESM_DT overrides dt
! The dump carries the parameter tables REDPRM reads (VEGPARM/SOILPARM/GENPARM
! as WRF loaded them); the driver copies them into module_sf_noahlsm, so no
! table file is read.  SFLX is an in-place step: CMC, T1, STC, SMC, SH2O,
! SNOWH, SNEQV advance over dt, and the fluxes are the ones used for that step.
! The state tendency (x_out - x)/dt converges to the instantaneous rate as dt
! shrinks (EqWeFiC PLAN.md section 6).  The Noah driver passes an
! uninitialised local (DUMMY) for SFCSPD, COSZ, PRCPRAIN, SOLARDIRECT and CM;
! SFLX reads none of them, and the replay passes zero.
program noah_driver
  use esm_flat_input
  use module_esm_dump
  use module_sf_noahlsm
  implicit none
  character(len=1024) :: inpath, outpath, dtstr
  integer :: nsoil, stat, vegtyp, soiltyp, slopetyp, isurban, opt_thcnd, fasdas, nroot, iiloc, jjloc
  logical :: local, rdlai2d, usemonalb, ua_phys
  real :: dt, ffrozp, zlvl, lwdn, soldn, solnet, sfcprs, prcp, sfctmp, q2k, th2, q2sat, dqsdt2
  real :: shdfac, shdmin, shdmax, albbrd, snoalb, tbot, z0brd, z0, emissi, embrd
  real :: cmc, t1, snowh, sneqv, albedo, ch, cm, xlai, snotime1, ribb, sfcheadrt, infxsrt, etpnd1, aoasis
  real :: sfcspd, cosz, prcprain, solardirect
  real :: eta, sheat, eta_kinematic, fdown, ec, edir, ett, esnow, drip, dew, beta, etp, ssoil
  real :: flx1, flx2, flx3, flx4, fvb, fbur, fgsn, snomlt, sncovr, runoff1, runoff2, runoff3
  real :: rc, pc, rsmin, rcs, rct, rcq, rcsoil, soilw, soilm, q1, smcwlt, smcdry, smcref, smcmax
  real :: xsda_qfx, hfx_phy, qfx_phy, xqnorm, hcpct_fasdas, irrigation_channel
  real, allocatable, dimension(:) :: sldpth, stc, smc, sh2o, et, smav
  character(len=256) :: llanduse, lsoil

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  nsoil = flat_i0('nsoil'); iiloc = flat_i0('i'); jjloc = flat_i0('j')
  dt = flat_r0('dt')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) dt
  allocate(sldpth(nsoil), stc(nsoil), smc(nsoil), sh2o(nsoil), et(nsoil), smav(nsoil))

  ! parameter tables (module_sf_noahlsm variables REDPRM reads)
  call flat_i1('tbl_nrotbl', NROTBL)
  call flat_r1('tbl_snuptbl', SNUPTBL); call flat_r1('tbl_rstbl', RSTBL); call flat_r1('tbl_rgltbl', RGLTBL)
  call flat_r1('tbl_hstbl', HSTBL); call flat_r1('tbl_shdtbl', SHDTBL); call flat_r1('tbl_maxalb', MAXALB)
  call flat_r1('tbl_emissmintbl', EMISSMINTBL); call flat_r1('tbl_emissmaxtbl', EMISSMAXTBL)
  call flat_r1('tbl_laimintbl', LAIMINTBL); call flat_r1('tbl_laimaxtbl', LAIMAXTBL)
  call flat_r1('tbl_z0mintbl', Z0MINTBL); call flat_r1('tbl_z0maxtbl', Z0MAXTBL)
  call flat_r1('tbl_albedomintbl', ALBEDOMINTBL); call flat_r1('tbl_albedomaxtbl', ALBEDOMAXTBL)
  call flat_r1('tbl_ztopvtbl', ZTOPVTBL); call flat_r1('tbl_zbotvtbl', ZBOTVTBL)
  TOPT_DATA = flat_r0('topt_data'); CMCMAX_DATA = flat_r0('cmcmax_data')
  CFACTR_DATA = flat_r0('cfactr_data'); RSMAX_DATA = flat_r0('rsmax_data')
  call flat_r1('tbl_bb', BB); call flat_r1('tbl_drysmc', DRYSMC); call flat_r1('tbl_f11', F11)
  call flat_r1('tbl_maxsmc', MAXSMC); call flat_r1('tbl_refsmc', REFSMC); call flat_r1('tbl_satpsi', SATPSI)
  call flat_r1('tbl_satdk', SATDK); call flat_r1('tbl_satdw', SATDW); call flat_r1('tbl_wltsmc', WLTSMC)
  call flat_r1('tbl_qtz', QTZ); call flat_r1('tbl_slope_data', SLOPE_DATA)
  SBETA_DATA = flat_r0('sbeta_data'); FXEXP_DATA = flat_r0('fxexp_data'); CSOIL_DATA = flat_r0('csoil_data')
  SALP_DATA = flat_r0('salp_data'); REFDK_DATA = flat_r0('refdk_data'); REFKDT_DATA = flat_r0('refkdt_data')
  FRZK_DATA = flat_r0('frzk_data'); ZBOT_DATA = flat_r0('zbot_data'); SMLOW_DATA = flat_r0('smlow_data')
  SMHIGH_DATA = flat_r0('smhigh_data'); CZIL_DATA = flat_r0('czil_data'); LVCOEF_DATA = flat_r0('lvcoef_data')
  LUCATS = flat_i0('lucats'); BARE = flat_i0('bare'); NATURAL = flat_i0('natural')
  SLCATS = flat_i0('slcats'); SLPCATS = flat_i0('slpcats')

  ! point inputs
  ffrozp = flat_r0('ffrozp'); isurban = flat_i0('isurban'); zlvl = flat_r0('zlvl')
  call flat_r1('sldpth', sldpth); local = flat_l0('local')
  llanduse = 'USGS'; lsoil = 'STAS'   ! strings are not carried by the flat format; the SCM run uses these
  lwdn = flat_r0('lwdn'); soldn = flat_r0('soldn'); solnet = flat_r0('solnet'); sfcprs = flat_r0('sfcprs')
  prcp = flat_r0('prcp'); sfctmp = flat_r0('sfctmp'); q2k = flat_r0('q2k'); th2 = flat_r0('th2')
  q2sat = flat_r0('q2sat'); dqsdt2 = flat_r0('dqsdt2')
  vegtyp = flat_i0('vegtyp'); soiltyp = flat_i0('soiltyp'); slopetyp = flat_i0('slopetyp')
  shdfac = flat_r0('shdfac'); shdmin = flat_r0('shdmin'); shdmax = flat_r0('shdmax')
  albbrd = flat_r0('albbrd'); snoalb = flat_r0('snoalb'); tbot = flat_r0('tbot')
  z0brd = flat_r0('z0brd'); z0 = flat_r0('z0'); emissi = flat_r0('emissi'); embrd = flat_r0('embrd')
  cmc = flat_r0('cmc'); t1 = flat_r0('t1'); call flat_r1('stc', stc); call flat_r1('smc', smc)
  call flat_r1('sh2o', sh2o); snowh = flat_r0('snowh'); sneqv = flat_r0('sneqv')
  albedo = flat_r0('albedo'); ch = flat_r0('ch'); xlai = flat_r0('xlai')
  rdlai2d = flat_l0('rdlai2d'); usemonalb = flat_l0('usemonalb'); ua_phys = flat_l0('ua_phys')
  snotime1 = flat_r0('snotime1'); ribb = flat_r0('ribb'); opt_thcnd = flat_i0('opt_thcnd')
  fasdas = flat_i0('fasdas'); irrigation_channel = flat_r0('irrigation_channel')
  sfcheadrt = flat_r0('sfcheadrt'); aoasis = flat_r0('aoasis')
  infxsrt = 0; etpnd1 = 0; xsda_qfx = 0; hfx_phy = 0; qfx_phy = 0; xqnorm = 0
  sfcspd = 0; cosz = 0; prcprain = 0; solardirect = 0; cm = 0

  ! EqWeFiC cold-season instrumentation: the SNOPAC and frozen-soil publishers
  ! are module variables, so clear them before the call -- a layer that SNKSRC
  ! does not visit must read back as "not called", not as the previous call's value.
  ESM_SNO_BRANCH = 0
  ESM_SNO_T12 = 0.0; ESM_SNO_T12A = 0.0; ESM_SNO_T12B = 0.0
  ESM_SNO_DENOM = 0.0; ESM_SNO_DTOT = 0.0; ESM_SNO_DSOIL = 0.0
  ESM_SNO_DF1 = 0.0; ESM_SNO_YY = 0.0; ESM_SNO_ZZ1 = 0.0
  ESM_SNO_ETANRG = 0.0; ESM_SNO_ESNOW1 = 0.0; ESM_SNO_ESNOW2 = 0.0
  ESM_SNO_EX = 0.0; ESM_SNO_PRCP1 = 0.0; ESM_SNO_SNDENS = 0.0
  ESM_SNO_ESD = 0.0; ESM_SNO_SNOMLT = 0.0; ESM_SNO_FLX1 = 0.0; ESM_SNO_FLX3 = 0.0
  ESM_FRZ_TAVG = 0.0; ESM_FRZ_TBND = 0.0; ESM_FRZ_TSNSR = 0.0
  ESM_FRZ_FREE = 0.0; ESM_FRZ_XH2O = 0.0; ESM_FRZ_QTOT = 0.0
  ESM_FRZ_CALLED = 0; ESM_FRZ_NLOG = 0; ESM_FRZ_PATH = 0; ESM_FRZ_TSURF = 0.0
  ESM_FRZ_SH2O_IN(1:nsoil) = sh2o(1:nsoil)

  call SFLX(IILOC=iiloc, JJLOC=jjloc, FFROZP=ffrozp, ISURBAN=isurban, DT=dt, ZLVL=zlvl, NSOIL=nsoil, &
       SLDPTH=sldpth, LOCAL=local, LLANDUSE=llanduse, LSOIL=lsoil, LWDN=lwdn, SOLDN=soldn, SOLNET=solnet, &
       SFCPRS=sfcprs, PRCP=prcp, SFCTMP=sfctmp, Q2=q2k, SFCSPD=sfcspd, COSZ=cosz, PRCPRAIN=prcprain, &
       SOLARDIRECT=solardirect, TH2=th2, Q2SAT=q2sat, DQSDT2=dqsdt2, VEGTYP=vegtyp, SOILTYP=soiltyp, &
       SLOPETYP=slopetyp, SHDFAC=shdfac, SHDMIN=shdmin, SHDMAX=shdmax, ALB=albbrd, SNOALB=snoalb, &
       TBOT=tbot, Z0BRD=z0brd, Z0=z0, EMISSI=emissi, EMBRD=embrd, CMC=cmc, T1=t1, STC=stc, SMC=smc, &
       SH2O=sh2o, SNOWH=snowh, SNEQV=sneqv, ALBEDO=albedo, CH=ch, CM=cm, ETA=eta, SHEAT=sheat, &
       ETA_KINEMATIC=eta_kinematic, FDOWN=fdown, EC=ec, EDIR=edir, ET=et, ETT=ett, ESNOW=esnow, DRIP=drip, &
       DEW=dew, BETA=beta, ETP=etp, SSOIL=ssoil, FLX1=flx1, FLX2=flx2, FLX3=flx3, FLX4=flx4, FVB=fvb, &
       FBUR=fbur, FGSN=fgsn, UA_PHYS=ua_phys, SNOMLT=snomlt, SNCOVR=sncovr, RUNOFF1=runoff1, &
       RUNOFF2=runoff2, RUNOFF3=runoff3, RC=rc, PC=pc, RSMIN=rsmin, XLAI=xlai, RCS=rcs, RCT=rct, RCQ=rcq, &
       RCSOIL=rcsoil, SOILW=soilw, SOILM=soilm, Q1=q1, SMAV=smav, RDLAI2D=rdlai2d, USEMONALB=usemonalb, &
       SNOTIME1=snotime1, RIBB=ribb, SMCWLT=smcwlt, SMCDRY=smcdry, SMCREF=smcref, SMCMAX=smcmax, &
       NROOT=nroot, SFHEAD1RT=sfcheadrt, INFXS1RT=infxsrt, ETPND1=etpnd1, OPT_THCND=opt_thcnd, &
       AOASIS=aoasis, XSDA_QFX=xsda_qfx, HFX_PHY=hfx_phy, QFX_PHY=qfx_phy, XQNORM=xqnorm, fasdas=fasdas, &
       HCPCT_FASDAS=hcpct_fasdas, IRRIGATION_CHANNEL=irrigation_channel)

  call esm_dump_open_file('noah', trim(outpath))
  call esm_dump_var('dt', dt); call esm_dump_var('real_bytes', storage_size(dt)/8)
  call esm_dump_var('cmc_out', cmc); call esm_dump_var('t1_out', t1); call esm_dump_var('stc_out', stc)
  call esm_dump_var('smc_out', smc); call esm_dump_var('sh2o_out', sh2o); call esm_dump_var('snowh_out', snowh)
  call esm_dump_var('sneqv_out', sneqv); call esm_dump_var('albedo_out', albedo); call esm_dump_var('ch_out', ch)
  call esm_dump_var('emissi_out', emissi); call esm_dump_var('z0_out', z0); call esm_dump_var('z0brd_out', z0brd)
  call esm_dump_var('eta', eta); call esm_dump_var('sheat', sheat); call esm_dump_var('eta_kinematic', eta_kinematic)
  call esm_dump_var('fdown', fdown); call esm_dump_var('ec', ec); call esm_dump_var('edir', edir)
  call esm_dump_var('et', et); call esm_dump_var('ett', ett); call esm_dump_var('esnow', esnow)
  call esm_dump_var('drip', drip); call esm_dump_var('dew', dew); call esm_dump_var('beta', beta)
  call esm_dump_var('etp', etp); call esm_dump_var('ssoil', ssoil)
  call esm_dump_var('flx1', flx1); call esm_dump_var('flx2', flx2); call esm_dump_var('flx3', flx3)
  call esm_dump_var('snomlt', snomlt); call esm_dump_var('sncovr', sncovr)
  call esm_dump_var('runoff1', runoff1); call esm_dump_var('runoff2', runoff2); call esm_dump_var('runoff3', runoff3)
  call esm_dump_var('rc', rc); call esm_dump_var('pc', pc); call esm_dump_var('rsmin', rsmin)
  call esm_dump_var('xlai_out', xlai); call esm_dump_var('rcs', rcs); call esm_dump_var('rct', rct)
  call esm_dump_var('rcq', rcq); call esm_dump_var('rcsoil', rcsoil)
  call esm_dump_var('soilw', soilw); call esm_dump_var('soilm', soilm); call esm_dump_var('q1', q1)
  call esm_dump_var('smav', smav); call esm_dump_var('snotime1_out', snotime1)
  call esm_dump_var('smcwlt', smcwlt); call esm_dump_var('smcdry', smcdry); call esm_dump_var('smcref', smcref)
  call esm_dump_var('smcmax', smcmax); call esm_dump_var('nroot', nroot)
  ! SFLX/NOPAC internals published by the noah_esm_*.inc splices (see the Makefile):
  ! the soil thermal conductivity of the top layer, the surface-energy-balance closure
  ! pair yy/zz1, and the linearised Penman quantities, none of which WRF stores.
  call esm_dump_var('df1', ESM_DF1); call esm_dump_var('yy', ESM_YY)
  call esm_dump_var('zz1', ESM_ZZ1); call esm_dump_var('ssoil_pre', ESM_SSOIL_PRE)
  call esm_dump_var('rch', ESM_RCH); call esm_dump_var('rr', ESM_RR)
  call esm_dump_var('epsca', ESM_EPSCA); call esm_dump_var('t24', ESM_T24)
  call esm_dump_var('fdown_pen', ESM_FDOWN)
  ! SNOPAC's snow-surface energy-balance closure and snow-pack internals, and the
  ! frozen-soil sink per layer (HRT -> TBND/TMPAVG/SNKSRC -> FRH2O).  sno_branch is
  ! 0 when SNOPAC was not called, 1 for its sub-freezing block, 2 for its melt block;
  ! frz_path is 1 when FRH2O returned FREE = SMC, 2 for the converged Newton
  ! iteration and 3 for the Flerchinger explicit fallback (the 10-iteration cap).
  call esm_dump_var('sno_branch', ESM_SNO_BRANCH)
  call esm_dump_var('sno_t12', ESM_SNO_T12); call esm_dump_var('sno_t12a', ESM_SNO_T12A)
  call esm_dump_var('sno_t12b', ESM_SNO_T12B); call esm_dump_var('sno_denom', ESM_SNO_DENOM)
  call esm_dump_var('sno_dtot', ESM_SNO_DTOT); call esm_dump_var('sno_dsoil', ESM_SNO_DSOIL)
  call esm_dump_var('sno_df1', ESM_SNO_DF1); call esm_dump_var('sno_yy', ESM_SNO_YY)
  call esm_dump_var('sno_zz1', ESM_SNO_ZZ1); call esm_dump_var('sno_etanrg', ESM_SNO_ETANRG)
  call esm_dump_var('sno_esnow1', ESM_SNO_ESNOW1); call esm_dump_var('sno_esnow2', ESM_SNO_ESNOW2)
  call esm_dump_var('sno_ex', ESM_SNO_EX); call esm_dump_var('sno_prcp1', ESM_SNO_PRCP1)
  call esm_dump_var('sno_sndens', ESM_SNO_SNDENS); call esm_dump_var('sno_esd', ESM_SNO_ESD)
  call esm_dump_var('sno_snomlt', ESM_SNO_SNOMLT)
  call esm_dump_var('sno_flx1', ESM_SNO_FLX1); call esm_dump_var('sno_flx3', ESM_SNO_FLX3)
  call esm_dump_var('frz_tsurf', ESM_FRZ_TSURF)
  call esm_dump_var('frz_tavg', ESM_FRZ_TAVG(1:nsoil))
  call esm_dump_var('frz_tbnd', ESM_FRZ_TBND(1:nsoil))
  call esm_dump_var('frz_tsnsr', ESM_FRZ_TSNSR(1:nsoil))
  call esm_dump_var('frz_free', ESM_FRZ_FREE(1:nsoil))
  call esm_dump_var('frz_xh2o', ESM_FRZ_XH2O(1:nsoil))
  call esm_dump_var('frz_qtot', ESM_FRZ_QTOT(1:nsoil))
  call esm_dump_var('frz_sh2o_in', ESM_FRZ_SH2O_IN(1:nsoil))
  call esm_dump_var('frz_called', ESM_FRZ_CALLED(1:nsoil))
  call esm_dump_var('frz_nlog', ESM_FRZ_NLOG(1:nsoil))
  call esm_dump_var('frz_path', ESM_FRZ_PATH(1:nsoil))
  call esm_dump_close()
end program noah_driver
