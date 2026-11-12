!=======================================================================
! module_esm_mesa_tape: EarthSciML instrumentation of MOSAIC's MESA
! pseudo-transient solver, SUBROUTINE mesa_ptc (module_mosaic_therm.F), at
! EVERY iteration and at the routine's own exit (EqWeFiC N113: tape a stage at
! the routine's own exit).
!
! Inert unless esm_mp_on is .true.; only kernels/mosaic_therm_driver sets it,
! so the model's own runs are unaffected.  Records are written as text, one
! line per record, to unit esm_mp_unit:
!   TAG case call it v1 v2 ...      (values es24.16e3, real(kind=8))
! TAGs (ion species are aer(i,...) for i in iso4 ino3 icl inh4 imsa ico3 ina
! ica ioin; E = electrolyte(1:22,...)):
!   ENTRY  after mesa_ptc's initialisation (it = 0): t_k rh_pc ah2o_a
!          mass_dry_a vol_dry_a mass_dry_salt keq_sl(15) keq_ll(3)
!          jsalt_present(15) epercent(22,jtotal) E(jtotal) aer(ions,jtotal)
!   AERT   aer(1:naer,jtotal) (it = 0)
!   GAMZ   log_gamz(22,22) column-major (it = 0)
!   TOP    loop top: E(jsolid) E(jliquid) aer(ions,jsolid) aer(ions,jliquid)
!          phi_salt_old alpha_salt phi_salt water_a
!   THERM  in mesa_flux_salt right after `calcium =` (ions_to_electrolytes,
!          compute_activities and the first mesa_estimate_eleliquid done):
!          aer(ions,jliquid) aer(ions,jsolid) E(jliquid) E(jsolid) water_a
!          activity(22) eleliquid(22)
!   FLUX   after mesa_convergence_criterion: iconverge_mass iconverge_flux
!          idissolved sat_ratio phi_salt flux_sl frac_salt_solid
!          frac_salt_liq E(jsolid) aer(ions,jliquid) aer(ions,jsolid)
!          eleliquid water_a activity(22)
!   STEP   end of the iteration: phi_bar alpha_salt phi_salt_old hsalt
!          flux_sl phi_salt E(jsolid) aer(ions,jsolid) aer(ions,jliquid)
!   EXIT   outcome (1 mass, 2 flux, 3 not converged) jaerosolstate
!          jhyst_leg growth_factor mass_wet_a vol_wet_a water_a mass_dry_a
!          E(jsolid) E(jliquid) aer(ions,jsolid) aer(ions,jliquid)
!          epercent(22,jsolid)
!=======================================================================
      module module_esm_mesa_tape
      use module_data_mosaic_therm
      implicit none
      logical, save :: esm_mp_on = .false.
      integer, save :: esm_mp_unit = 0, esm_mp_case = 0, esm_mp_call = 0
      integer, save :: esm_mp_it = 0
      logical, save :: esm_at_on = .false., esm_at_file = .false.
      integer, save :: esm_at_unit = 0, esm_at_sub = 0
      character(len=*), parameter :: esm_mp_fmt = '(a,3(1x,i0),*(1x,es24.16e3))'

      contains

      subroutine esm_mp_ions(jp, ibin, v)
      integer jp, ibin
      real(kind=8) v(9)
      v = (/ aer(iso4_a,jp,ibin), aer(ino3_a,jp,ibin), aer(icl_a,jp,ibin), &
             aer(inh4_a,jp,ibin), aer(imsa_a,jp,ibin), aer(ico3_a,jp,ibin), &
             aer(ina_a,jp,ibin), aer(ica_a,jp,ibin), aer(ioin_a,jp,ibin) /)
      end subroutine esm_mp_ions

      subroutine esm_mp_entry(ibin)
      integer ibin, i
      real(kind=8) vt(9)
      esm_mp_call = esm_mp_call + 1
      esm_mp_it = 0
      call esm_mp_ions(jtotal, ibin, vt)
      write(esm_mp_unit, esm_mp_fmt) 'ENTRY', esm_mp_case, esm_mp_call, 0, &
        t_k, rh_pc, ah2o_a(ibin), mass_dry_a(ibin), vol_dry_a(ibin), &
        mass_dry_salt(ibin), keq_sl(1:nsalt), keq_ll(1:3), &
        (dble(jsalt_present(i)), i = 1, nsalt), epercent(1:nelectrolyte,jtotal,ibin), &
        electrolyte(1:nelectrolyte,jtotal,ibin), vt
      write(esm_mp_unit, esm_mp_fmt) 'AERT', esm_mp_case, esm_mp_call, 0, &
        aer(1:naer,jtotal,ibin)
      write(esm_mp_unit, esm_mp_fmt) 'GAMZ', esm_mp_case, esm_mp_call, 0, &
        log_gamz(1:nelectrolyte,1:nelectrolyte)
      end subroutine esm_mp_entry

      subroutine esm_mp_top(ibin, it)
      integer ibin, it
      real(kind=8) vs(9), vl(9)
      esm_mp_it = it
      call esm_mp_ions(jsolid, ibin, vs)
      call esm_mp_ions(jliquid, ibin, vl)
      write(esm_mp_unit, esm_mp_fmt) 'TOP', esm_mp_case, esm_mp_call, it, &
        electrolyte(1:nelectrolyte,jsolid,ibin), electrolyte(1:nelectrolyte,jliquid,ibin), &
        vs, vl, phi_salt_old(1:nsalt), alpha_salt(1:nsalt), phi_salt(1:nsalt), water_a(ibin)
      end subroutine esm_mp_top

      subroutine esm_mp_therm(ibin)
      integer ibin
      real(kind=8) vs(9), vl(9)
      call esm_mp_ions(jsolid, ibin, vs)
      call esm_mp_ions(jliquid, ibin, vl)
      write(esm_mp_unit, esm_mp_fmt) 'THERM', esm_mp_case, esm_mp_call, esm_mp_it, &
        vl, vs, electrolyte(1:nelectrolyte,jliquid,ibin), electrolyte(1:nelectrolyte,jsolid,ibin), &
        water_a(ibin), activity(1:nelectrolyte,ibin), eleliquid(1:nelectrolyte)
      end subroutine esm_mp_therm

      subroutine esm_mp_flux(ibin, it, icm, icf, idis)
      integer ibin, it, icm, icf, idis
      real(kind=8) vs(9), vl(9)
      call esm_mp_ions(jsolid, ibin, vs)
      call esm_mp_ions(jliquid, ibin, vl)
      write(esm_mp_unit, esm_mp_fmt) 'FLUX', esm_mp_case, esm_mp_call, it, &
        dble(icm), dble(icf), dble(idis), sat_ratio(1:nsalt), phi_salt(1:nsalt), &
        flux_sl(1:nsalt), frac_salt_solid(1:nsalt), frac_salt_liq(1:nsalt), &
        electrolyte(1:nelectrolyte,jsolid,ibin), vl, vs, eleliquid(1:nelectrolyte), &
        water_a(ibin), activity(1:nelectrolyte,ibin)
      end subroutine esm_mp_flux

      subroutine esm_mp_step(ibin, it)
      integer ibin, it
      real(kind=8) vs(9), vl(9)
      call esm_mp_ions(jsolid, ibin, vs)
      call esm_mp_ions(jliquid, ibin, vl)
      write(esm_mp_unit, esm_mp_fmt) 'STEP', esm_mp_case, esm_mp_call, it, &
        phi_bar(1:nsalt), alpha_salt(1:nsalt), phi_salt_old(1:nsalt), hsalt(1:nsalt), &
        flux_sl(1:nsalt), phi_salt(1:nsalt), electrolyte(1:nelectrolyte,jsolid,ibin), vs, vl
      end subroutine esm_mp_step

      subroutine esm_mp_exit(ibin, it, outcome)
      integer ibin, it, outcome
      real(kind=8) vs(9), vl(9)
      call esm_mp_ions(jsolid, ibin, vs)
      call esm_mp_ions(jliquid, ibin, vl)
      write(esm_mp_unit, esm_mp_fmt) 'EXIT', esm_mp_case, esm_mp_call, it, &
        dble(outcome), dble(jaerosolstate(ibin)), dble(jhyst_leg(ibin)), growth_factor(ibin), &
        mass_wet_a(ibin), vol_wet_a(ibin), water_a(ibin), mass_dry_a(ibin), &
        electrolyte(1:nelectrolyte,jsolid,ibin), electrolyte(1:nelectrolyte,jliquid,ibin), &
        vs, vl, epercent(1:nelectrolyte,jsolid,ibin)
      end subroutine esm_mp_exit

!=======================================================================
! ASTEM tape (2026-11-12): MOSAIC's gas-particle partitioning solver,
! SUBROUTINE ASTEM and its stages, recorded at each stage's OWN exit
! (EqWeFiC N113).  Inert unless esm_at_on; the model turns it on for a
! dumped column (esm_rate_on) only when the environment variable
! ESM_ASTEM_TAPE names a file (esm_at_open), and the therm driver may set it
! directly.  Same text layout as the MESA tape, header
!   TAG step lev sub bin v1 v2 ...
! with step = ncorecnt_aer, lev = kclm_aer, sub = ASTEM sub-step (0 before
! the semi-volatile loop), bin = ibin (0 for whole-box records).
! Ions (9) as esm_mp_ions; G = gas(1:ngas_ioa); E = electrolyte(1:22,jp).
!   AIN   after aerosolmtc: dtchem t_k p_atm rh_pc ah2o keq_gl(4) keq_ll(3)
!         kp_nh4no3 kp_nh4cl G kp_nh4no3_0 kp_nh4cl_0 keq_sg(1:2)
!   BIN   (per bin, same point) jaerosolstate jphase jhyst_leg num_a dp_wet_a
!         water_a mass_dry_a kg(1:5) kel(1:5) ions(s) ions(l) ions(t)
!         E(s) E(l) E(t) epercent(1:22,jtotal) aer(1:naer,jtotal)
!   SB    (state per bin) jaerosolstate jphase jhyst_leg water_a ions(s)
!         ions(l) ions(t) E(s) E(l) E(t)
!   SVT   sub-step top (label 10): t_old G total_species(1:5); per bin SVB = SB
!   (Only the PHYSICS is taped: the per-bin surface equilibria and driving
!   forces where the fluxes are formed; ASTEM's step updates are not, because
!   EqWeFiC does not transcribe the integrator.)
!   FXI   per bin before its flux call: route(1 dry 2 wet 3 mix 4 mixed->wet) G
!   FXO   per bin after: route jphase idry_case3a ieqblm sfc_a(1:5)
!         df_gas_s df_gas_l flux_s flux_l Heff phi_volatile_s phi_volatile_l
!         integrate(1:5,jsolid) integrate(1:5,jliquid) mc(1:4) ma(1:5)
!         gam(1:22) gam_ratio activity(1:22) water_a keq_nh4no3 keq_nh4cl
!         G ions(s) ions(l) ions(t) E(s) E(l) kel(1:5) epercent(1:22,jsolid)
!         keq_nh4no3_0 keq_nh4cl_0
!   AOUT  end of ASTEM_semi_volatiles: isteps G
!=======================================================================
      subroutine esm_at_open()
      character(len=1024) :: path
      integer :: st
      logical, save :: tried = .false.
      if (tried) return
      tried = .true.
      call get_environment_variable('ESM_ASTEM_TAPE', path, status=st)
      if (st /= 0 .or. len_trim(path) == 0) return
      esm_at_unit = 78
      open(esm_at_unit, file=trim(path), status='replace')
      esm_at_file = .true.
      end subroutine esm_at_open

      subroutine esm_at_hdr(tag, sub, ibin, v)
      character(len=*) tag
      integer sub, ibin
      real(kind=8) v(:)
      write(esm_at_unit, '(a,4(1x,i0),*(1x,es24.16e3))') tag, ncorecnt_aer, &
        kclm_aer, sub, ibin, v
      end subroutine esm_at_hdr

      subroutine esm_at_ain(dtchem)
      real(kind=8) dtchem
      integer ibin
      real(kind=8) vs(9), vl(9), vt(9)
      if (.not. esm_at_on) return
      esm_at_sub = 0
      call esm_at_hdr('AIN', 0, 0, (/ dtchem, t_k, p_atm, rh_pc, ah2o, &
        keq_gl(1:4), keq_ll(1:3), kp_nh4no3, kp_nh4cl, gas(1:ngas_ioa), &
        kp_nh4no3_0, kp_nh4cl_0, keq_sg(1:2) /))
      do ibin = 1, nbin_a
        call esm_mp_ions(jsolid, ibin, vs); call esm_mp_ions(jliquid, ibin, vl)
        call esm_mp_ions(jtotal, ibin, vt)
        call esm_at_hdr('BIN', 0, ibin, (/ dble(jaerosolstate(ibin)), &
          dble(jphase(ibin)), dble(jhyst_leg(ibin)), num_a(ibin), dp_wet_a(ibin), &
          water_a(ibin), mass_dry_a(ibin), kg(1:ngas_ioa,ibin), kel(1:ngas_ioa,ibin), &
          vs, vl, vt, electrolyte(1:nelectrolyte,jsolid,ibin), &
          electrolyte(1:nelectrolyte,jliquid,ibin), electrolyte(1:nelectrolyte,jtotal,ibin), &
          epercent(1:nelectrolyte,jtotal,ibin), aer(1:naer,jtotal,ibin) /))
      end do
      end subroutine esm_at_ain

      subroutine esm_at_sb(tag, ibin)
      character(len=*) tag
      integer ibin
      real(kind=8) vs(9), vl(9), vt(9)
      call esm_mp_ions(jsolid, ibin, vs); call esm_mp_ions(jliquid, ibin, vl)
      call esm_mp_ions(jtotal, ibin, vt)
      call esm_at_hdr(tag, esm_at_sub, ibin, (/ dble(jaerosolstate(ibin)), &
        dble(jphase(ibin)), dble(jhyst_leg(ibin)), water_a(ibin), vs, vl, vt, &
        electrolyte(1:nelectrolyte,jsolid,ibin), electrolyte(1:nelectrolyte,jliquid,ibin), &
        electrolyte(1:nelectrolyte,jtotal,ibin) /))
      end subroutine esm_at_sb



      subroutine esm_at_svt(isteps, t_old)
      integer isteps, ibin
      real(kind=8) t_old
      if (.not. esm_at_on) return
      esm_at_sub = isteps
      call esm_at_hdr('SVT', isteps, 0, (/ t_old, gas(1:ngas_ioa), &
        total_species(1:ngas_ioa) /))
      do ibin = 1, nbin_a
        call esm_at_sb('SVB', ibin)
      end do
      end subroutine esm_at_svt

      subroutine esm_at_fxi(ibin, route)
      integer ibin, route
      if (.not. esm_at_on) return
      call esm_at_hdr('FXI', esm_at_sub, ibin, (/ dble(route), gas(1:ngas_ioa) /))
      end subroutine esm_at_fxi

      subroutine esm_at_fxo(ibin, route)
      integer ibin, route, iv
      real(kind=8) vs(9), vl(9), vt(9)
      if (.not. esm_at_on) return
      call esm_mp_ions(jsolid, ibin, vs); call esm_mp_ions(jliquid, ibin, vl)
      call esm_mp_ions(jtotal, ibin, vt)
      call esm_at_hdr('FXO', esm_at_sub, ibin, (/ dble(route), dble(jphase(ibin)), &
        dble(idry_case3a(ibin)), dble(ieqblm_astem), sfc_a(1:ngas_ioa), &
        df_gas_s(1:ngas_ioa,ibin), df_gas_l(1:ngas_ioa,ibin), flux_s(1:ngas_ioa,ibin), &
        flux_l(1:ngas_ioa,ibin), heff(1:ngas_ioa,ibin), phi_volatile_s(1:ngas_ioa,ibin), &
        phi_volatile_l(1:ngas_ioa,ibin), (dble(integrate(iv,jsolid,ibin)), iv=1,ngas_ioa), &
        (dble(integrate(iv,jliquid,ibin)), iv=1,ngas_ioa), mc(1:ncation,ibin), &
        ma(1:nanion,ibin), gam(1:nelectrolyte,ibin), gam_ratio(ibin), &
        activity(1:nelectrolyte,ibin), water_a(ibin), keq_nh4no3, keq_nh4cl, &
        gas(1:ngas_ioa), vs, vl, vt, electrolyte(1:nelectrolyte,jsolid,ibin), &
        electrolyte(1:nelectrolyte,jliquid,ibin), kel(1:ngas_ioa,ibin), &
        epercent(1:nelectrolyte,jsolid,ibin), keq_nh4no3_0, keq_nh4cl_0 /))
      end subroutine esm_at_fxo




      subroutine esm_at_aout(isteps)
      integer isteps
      if (.not. esm_at_on) return
      call esm_at_hdr('AOUT', isteps, 0, (/ dble(isteps), gas(1:ngas_ioa) /))
      end subroutine esm_at_aout

      end module module_esm_mesa_tape
