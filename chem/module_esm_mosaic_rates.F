!=======================================================================
! module_esm_mosaic_rates: EarthSciML instrumentation buffers for MOSAIC's
! INSTANTANEOUS sub-process rates.
!
! Why this exists.  Each MOSAIC stage updates the column state rsub in place in
! real32, so a stage increment recovered by differencing two dumped states is
! the difference of two nearly equal single-precision numbers -- for
! coagulation only one to two significant digits survive (EqWeFiC
! FORTRAN_BUGS N87).  These buffers let each process write the rate WHERE IT
! FORMS IT, at the state the stage starts from, so the reference never depends
! on differencing (or on a promoted real64 build).
!
! Control: module_mosaic_driver sets esm_rate_on for the dumped column and
! esm_rate_k to the rsub level index before each stage; the process modules
! fill their buffer when esm_rate_on is true; the driver dumps the buffers.
! Everything here is inert when esm_rate_on is .false.
!
! Units are the ones the process itself works in, and are stated per buffer.
!=======================================================================
      module module_esm_mosaic_rates

      implicit none
      public

      logical, save :: esm_rate_on = .false.
      integer, save :: esm_rate_k  = 0

      integer, parameter :: esm_maxk    = 100   ! kmaxd of module_data_mosaic_other
      integer, parameter :: esm_maxbin  = 8     ! maxd_asize
      integer, parameter :: esm_maxcomp = 40    ! coagulation components + 3
      integer, parameter :: esm_maxgas  = 10    ! MOSAIC volatile gases

!   nucleation.  With FORTRAN_BUGS B20 fixed the scheme is an instantaneous
!   PROJECTION (the vapour relaxes to the critical concentration and the excess
!   forms particles), so what is recorded is the projection increment itself,
!   in the units of rsub (mol/mol-air for gases and aerosol species,
!   particles/mol-air for number) PER CALL, not per second.
      real, save :: esm_rate_newnuc(6, esm_maxk)
!     rows: 1 h2so4, 2 nh3, 3 so4_a(bin), 4 nh4_a(bin), 5 num_a(bin),
!           6 the bin index the new particles go into (not an increment)

!   the same increments divided by dtnuc, which is what a consumer that insists
!   on a rate would see; with B20 fixed this scales as 1/dtnuc
      real, save :: esm_rate_newnuc_fd(6, esm_maxk)

!   the composition and thresholds B20 also corrupted:
!     1 (mole aerosol nh4)/(mole aerosol so4) of the new particles
!     2 their dry density (g/cm3)
!     3 qh2so4_crit, 4 qh2so4_cur, 5 qnh3_cur (mol/mol-air)
      real, save :: esm_rate_newnuc_extra(5, esm_maxk)

!   coagulation, evaluated at the state the stage starts from:
!     d(number)/dt   in the units of num_distrib per second
!     d(volume)/dt   in the units of vol_distrib per second
      real, save :: esm_rate_coag_num(esm_maxbin, esm_maxk)
      real, save :: esm_rate_coag_vol(esm_maxbin, esm_maxcomp, esm_maxk)

!   gas-particle transfer (ASTEM), evaluated at the state the stage starts
!   from, in MOSAIC's internal units (nmol/m3 air per second):
!     non-volatiles: gas(iv)*kg(iv,ibin), the exact dt -> 0 limit of the
!                    analytical exponential decay
!     semi-volatiles: flux_s + flux_l, the fluxes ASTEM itself forms on its
!                     first sub-step
      real, save :: esm_rate_cond_nonvol(esm_maxgas, esm_maxbin, esm_maxk)
      real, save :: esm_rate_cond_semivol(esm_maxgas, esm_maxbin, esm_maxk)

!   the mass-transfer coefficients ASTEM condenses with, and the two per-bin
!   quantities they are built from, recorded right after aerosolmtc:
!     esm_mtc_kg     kg(iv,ibin), 1/s
!     esm_mtc_dpwet  dp_wet_a(ibin), the wet diameter calc_dry_n_wet_aerosol_props
!                    gives, cm
!     esm_mtc_num    num_a(ibin), particles/cm3
!     esm_mtc_tk     t_k (K) and esm_mtc_patm p_atm (atm) of the cell
      real, save :: esm_mtc_kg(esm_maxgas, esm_maxbin, esm_maxk)
      real, save :: esm_mtc_dpwet(esm_maxbin, esm_maxk)
      real, save :: esm_mtc_num(esm_maxbin, esm_maxk)
      real, save :: esm_mtc_tk(esm_maxk), esm_mtc_patm(esm_maxk)

!   the Brownian coagulation kernel and the per-bin sizes it is built from,
!   recorded inside coagsolv once the kernel matrix is filled:
!     esm_coag_cbr      cbr(i,j), cm^3/(particle*s)
!     esm_coag_dpwet    dpwet_inp(i), cm, and esm_coag_denswet denswet_inp(i), g/cm^3
!     esm_coag_tk       tkelvin (K) and esm_coag_patm patm (atm) of the cell
      real, save :: esm_coag_cbr(esm_maxbin, esm_maxbin, esm_maxk)
      real, save :: esm_coag_dpwet(esm_maxbin, esm_maxk)
!
!     The inputs and the bin-splitting tables of coagsolv's PAIR ASSEMBLY -- the
!     step that turns the pair kernel cbr into each bin's own rate.  None of
!     these is recoverable from what the driver already dumps:
!       esm_coag_dpdry    dpdry_inp(i), cm.  volume(i) = (4pi/3)(dpdry/2)^3 is
!                         what the splitting tables and the number production
!                         are written in, and only dpwet was published before.
!       esm_coag_distrib  the distribution coagsolv actually integrates,
!                         distrib(i,1) = number (#/cm3-air) and
!                         distrib(i,jaer>1) = the ncomp_coag volume/mass mixing
!                         ratios.  Reconstructing it from rsub would mean
!                         re-deriving mosaic_coag_1clm's number-conforming
!                         clamp, which rewrites rsub BEFORE coagsolv sees it.
!       esm_coag_fink     fink(i,j,k) DENSE: the volume fraction of an i+j
!                         collision that lands in bin k.  coagsolv stores it
!                         sparsely, indexed by a per-(i,k) counter jbins with
!                         the partner bin in jbinik, so the dense array is the
!                         only form a consumer can check its own transcription
!                         of the splitting rule against.
!       esm_coag_floss    floss(i,j), the fraction of bin i lost per i+j
!                         collision.
!       esm_coag_deltat   deltat = deltat_inp/nsubstep, the sub-step rrate is
!                         scaled by.
      real, save :: esm_coag_dpdry(esm_maxbin, esm_maxk)
      real, save :: esm_coag_distrib(esm_maxbin, esm_maxcomp, esm_maxk)
      real, save :: esm_coag_fink(esm_maxbin, esm_maxbin, esm_maxbin, esm_maxk)
      real, save :: esm_coag_floss(esm_maxbin, esm_maxbin, esm_maxk)
      real, save :: esm_coag_deltat = 0.0
      real, save :: esm_coag_denswet(esm_maxbin, esm_maxk)
      real, save :: esm_coag_tk(esm_maxk), esm_coag_patm(esm_maxk)

!   MOSAIC THERMODYNAMICS, compute_activities (module_mosaic_therm.F): its inputs
!   and outputs at the FIRST and the LAST call per (bin, level) of a dumped step
!   (slot 1 and slot 2), so the ZSR water, the MTEM mixing rule and the ionic
!   equilibrium it performs can be referenced call by call.  The thermodynamics
!   runs in real(kind=8), so these are real(kind=8) as well.  esm_th_* names:
!       elec(je)    electrolyte(je, jliquid, ibin), nmol/m3(air)
!       aerl/aert   aer(ja, jliquid|jtotal, ibin), nmol/m3(air)
!       ah2o_a      ah2o_a(ibin), the Kelvin-adjusted water activity; ah2o (per
!                   level) the ambient relative humidity
!       water       water_a(ibin), kg/m3(air), set by aerosol_water(jliquid)
!       mc, ma      cation/anion molalities, mol/kg(water)
!       gam, act    gam(je, ibin) and activity(je, ibin)
!       xt          calculate_xt(ibin, jliquid)
!       log_gamz    the binary MTEM log10 gamma of A in E at ah2o, per level
!       keq_ll      Keq_ll(1:3), per level
!   and the parameter tables, filled once: a_zsr, b_zsr, aw_min, b_mtem,
!   kappa_nonelectro, mw_aer_mac, dens_aer_mac.
      integer, parameter :: esm_maxel = 22, esm_maxaer = 112
      integer, save :: esm_th_ncall(esm_maxbin, esm_maxk)
      integer, save :: esm_th_slot = 0
      real(kind=8), save :: esm_th_elec(esm_maxel, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_aerl(esm_maxaer, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_aert(esm_maxaer, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_ah2o_a(esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_water(esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_xt(esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_gam_ratio(esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_mc(4, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_ma(5, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_gam(esm_maxel, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_act(esm_maxel, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_th_ah2o(esm_maxk), esm_th_tk(esm_maxk)
      real(kind=8), save :: esm_th_keq_ll(3, esm_maxk)
      real(kind=8), save :: esm_th_log_gamz(esm_maxel, esm_maxel, esm_maxk)
      real(kind=8), save :: esm_th_a_zsr(6, esm_maxel), esm_th_b_zsr(esm_maxel)
      real(kind=8), save :: esm_th_aw_min(esm_maxel)
      real(kind=8), save :: esm_th_b_mtem(6, esm_maxel, esm_maxel)
      real(kind=8), save :: esm_th_kappa(esm_maxaer), esm_th_mw_aer(esm_maxaer)
      real(kind=8), save :: esm_th_dens_aer(esm_maxaer)
      integer, save :: esm_th_iaer(7) = 0   ! iso4_a ino3_a icl_a inh4_a imsa_a ina_a ica_a


!   MOSAIC PHASE STATE, aerosol_phase_state (module_mosaic_therm.F): the Kelvin
!   iteration around MESA, at the FIRST and the LAST call per (bin, level) of a
!   dumped step (slot 1 / slot 2), real(kind=8).  esm_ps_* names:
!     entry:  rh_pc, t_k (per level), aert_in = aer(:,jtotal), elet_in =
!             electrolyte(:,jtotal) (conform_electrolytes has run), water_in,
!             hyst_in = water_a_hyst, up_in = water_a_up, num = num_a,
!             jhyst_in = jhyst_leg
!     exit:   state = jaerosolstate, phase = jphase, jhyst = jhyst_leg,
!             path = which MESA branch decided the state (esm_ps_path_*),
!             niter = Kelvin iterations, ah2o_a, kelvin, kel(iv), water,
!             water_up, mass/vol dry/wet, dpmv, sigma, growth, mdrh,
!             aers/aerl = aer(:,jsolid|jliquid), eles/elel likewise
!   and once: sigma_water, partial_molar_vol, mdrh_t, ptol_mol_astem.
      integer, parameter :: esm_maxel23 = 23, esm_ps_maxgas = 112
      integer, save :: esm_ps_ncall(esm_maxbin, esm_maxk)
      integer, save :: esm_ps_slot = 0, esm_ps_path_now = 0
      real(kind=8), save :: esm_ps_rh_pc(esm_maxk), esm_ps_tk(esm_maxk)
      real(kind=8), save :: esm_ps_aert_in(esm_maxaer, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_elet_in(esm_maxel23, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_aers(esm_maxaer, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_aerl(esm_maxaer, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_eles(esm_maxel23, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_elel(esm_maxel23, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_kel(esm_ps_maxgas, esm_maxbin, esm_maxk, 2)
!     scalars per (bin, level, slot): 1 water_in 2 hyst_in 3 up_in 4 num
!     5 jhyst_in 6 state 7 phase 8 jhyst 9 path 10 niter 11 ah2o_a 12 kelvin
!     13 water 14 water_up 15 mass_dry 16 vol_dry 17 mass_wet 18 vol_wet
!     19 dpmv 20 sigma 21 growth 22 mdrh 23 rh_pc
      integer, parameter :: esm_ps_nsc = 23
      real(kind=8), save :: esm_ps_sc(esm_ps_nsc, esm_maxbin, esm_maxk, 2)
      real(kind=8), save :: esm_ps_sigma_water = 0.0d0
      real(kind=8), save :: esm_ps_pmv(esm_ps_maxgas) = 0.0d0
      real(kind=8), save :: esm_ps_mdrh_t(64) = 0.0d0
      real(kind=8), save :: esm_ps_ptol = 0.0d0

!   shapes actually used, so a consumer does not have to guess

!   section moving (move_sections, module_mosaic_movesect.F; msectional = 20,
!   the linear-discrete / mass-number-advection algorithm), recorded around
!   the call in aerchemistry (module_mosaic_therm.F) at the routine's OWN entry
!   and exit, per level of the dumped column:
!     esm_ms_rsub_bef / _aft  the whole rsub(1:ltot2) column vector just before
!                             and just after move_sections (mol/mol-air for
!                             aerosol species, particles/mol-air for number)
!     esm_ms_grow(1:4,n)      drymass_pregrow (g/mol-air), drydens_pregrow
!                             (g/cm3), drymass_aftgrow, drydens_aftgrow
!     esm_ms_diag(1:6,n)      after the call: adrydens_sub, awetdens_sub
!                             (g/cm3), admeandry_sub, admeanwet_sub (cm),
!                             aqvoldry_sub (cm3/mol-air), aqmassdry_sub (g/mol-air)
!     esm_ms_xfer(1:7,n)      from calc_masnumadvect: iforce_movecenter,
!                             nnewsave(1:2), xferfracvol(1:2), xferfracnum(1:2)
      integer, parameter :: esm_maxl = 400
      integer, save :: esm_ms_ltot = 0
      real, save :: esm_ms_rsub_bef(esm_maxl, esm_maxk)
      real, save :: esm_ms_rsub_aft(esm_maxl, esm_maxk)
      real, save :: esm_ms_grow(4, esm_maxbin, esm_maxk)
      real, save :: esm_ms_diag(6, esm_maxbin, esm_maxk)
      real, save :: esm_ms_xfer(7, esm_maxbin, esm_maxk)
      integer, save :: esm_nbin = 0, esm_ncomp_coag = 0, esm_ngas = 0
!   MOSAIC's gas indices, dumped with the rates
      integer, save :: esm_ih2so4 = 0, esm_ihno3 = 0, esm_ihcl = 0, &
                       esm_inh3 = 0, esm_imsa = 0

      contains

      subroutine esm_rates_reset()
      esm_rate_newnuc      = 0.0
      esm_rate_newnuc_fd   = 0.0
      esm_rate_newnuc_extra = 0.0
      esm_rate_coag_num    = 0.0
      esm_rate_coag_vol    = 0.0
      esm_rate_cond_nonvol = 0.0
      esm_rate_cond_semivol = 0.0
      esm_mtc_kg = 0.0
      esm_mtc_dpwet = 0.0
      esm_mtc_num = 0.0
      esm_mtc_tk = 0.0
      esm_mtc_patm = 0.0
      esm_coag_cbr = 0.0
      esm_coag_dpwet = 0.0
      esm_coag_dpdry = 0.0
      esm_coag_distrib = 0.0
      esm_coag_fink = 0.0
      esm_coag_floss = 0.0
      esm_coag_deltat = 0.0
      esm_coag_denswet = 0.0
      esm_coag_tk = 0.0
      esm_coag_patm = 0.0
      esm_th_ncall = 0; esm_th_slot = 0
      esm_th_elec = 0.0d0; esm_th_aerl = 0.0d0; esm_th_aert = 0.0d0
      esm_th_ah2o_a = 0.0d0; esm_th_water = 0.0d0; esm_th_xt = 0.0d0
      esm_th_gam_ratio = 0.0d0; esm_th_mc = 0.0d0; esm_th_ma = 0.0d0
      esm_th_gam = 0.0d0; esm_th_act = 0.0d0
      esm_th_ah2o = 0.0d0; esm_th_tk = 0.0d0; esm_th_keq_ll = 0.0d0
      esm_th_log_gamz = 0.0d0
      esm_ps_ncall = 0; esm_ps_slot = 0; esm_ps_path_now = 0
      esm_ps_rh_pc = 0.0d0; esm_ps_tk = 0.0d0
      esm_ps_aert_in = 0.0d0; esm_ps_elet_in = 0.0d0
      esm_ps_aers = 0.0d0; esm_ps_aerl = 0.0d0
      esm_ps_eles = 0.0d0; esm_ps_elel = 0.0d0
      esm_ps_kel = 0.0d0; esm_ps_sc = 0.0d0
      esm_ms_rsub_bef = 0.0; esm_ms_rsub_aft = 0.0
      esm_ms_grow = 0.0; esm_ms_diag = 0.0; esm_ms_xfer = 0.0
      end subroutine esm_rates_reset

      end module module_esm_mosaic_rates
