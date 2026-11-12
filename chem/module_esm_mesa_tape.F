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

      end module module_esm_mesa_tape
