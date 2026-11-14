! astem_section_driver: ASTEM's per-section fluxes AT ONE STATE (EqWeFiC per-section replay).
!
!   astem_section_driver <levels.txt> <out.tape>
!
! In the model, ASTEM_semi_volatiles forms the eight sections' fluxes one after another, and a
! section's flux routine can edit the gas (absorb/degas_tiny, degas_solid_*) before the next
! section forms its flux, so the later sections' fluxes in the in-model ASTEM tape are not at
! the sub-step's state.  This driver loads the state of the FIRST sub-step exactly as the
! model's tape recorded it (SVT gas, SVB per section, BIN kg, the Kelvin factors, the dumped
! ah2o_a and aer(:,jtotal)) and, for every section, restores the gas and the section's
! aerosol to that state and runs the model's own dispatch (ASTEM_flux_dry / _wet / _mix,
! module_mosaic_therm compiled from chem/, real(kind=8)) with the ASTEM tape on.  It does it
! twice: sub-step tag 1 with the model's state edits, tag 2 with them skipped
! (esm_sec_noedit, kernels/astem_section_patch.py).  Records: SVT/SVB (tag 1, a check that the
! loaded state is the model's), then FXI/FXO per section per tag; the tape header carries
! the model step (ncorecnt_aer) and level (kclm_aer) as the model's own tape does.
!
! levels.txt: nlev, then per level (list-directed):
!   step lev
!   t_k p_atm rh_pc ah2o relhumclm  (AIN; relhumclm from the dump: the model computes the
!                                    constants with ah2o = relhumclm, and aerosol_phase_state
!                                    then resets ah2o = rh_pc*0.01, the AIN value)
!   gas(5) total_species(5)        (SVT)
!   per section 1..8:
!     jaerosolstate jphase jhyst_leg
!     water_a ah2o_a num_a dp_wet_a mass_dry_a
!     kg(5) kel(5)
!     ions_s(9) ions_l(9)          (so4 no3 cl nh4 msa co3 na ca oin)
!     aer(1:naer, jtotal)
!     E_s(22) E_l(22) E_t(22)
program astem_section_driver
  use module_data_mosaic_therm
  use module_mosaic_therm
  use module_esm_mesa_tape
  use module_esm_mosaic_rates
  implicit none
  character(len=1024) :: inpath, outpath
  integer, parameter :: nb = 8
  integer :: nlev, il, ib, iv, step, lev, vbs(1), itag, js(nb), jp(nb), jh(nb), ia(9)
  real(kind=8) :: g0(5), tot0(5), tk, pa, rh, ah, rhc, xs(9), xl(9)
  real(kind=8) :: aer0(naer,3,nb), el0(nelectrolyte,3,nb), w0(nb)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  esm_at_unit = 79
  open(esm_at_unit, file=trim(outpath), status='replace')
  esm_at_file = .true.
  open(10, file=trim(inpath), status='old')
  read(10,*) nlev
  vbs(1) = 0
  call load_mosaic_parameters
  ia = (/ iso4_a, ino3_a, icl_a, inh4_a, imsa_a, ico3_a, ina_a, ica_a, ioin_a /)
  do il = 1, nlev
    read(10,*) step, lev
    read(10,*) tk, pa, rh, ah, rhc
    read(10,*) g0, tot0
    call initialize_mosaic_variables
    nbin_a = nb
    t_k = tk; p_atm = pa; rh_pc = rh; ah2o = rhc
    call update_thermodynamic_constants(vbs)
    ah2o = ah
    gas(1:5) = g0; total_species(1:5) = tot0
    do ib = 1, nb
      read(10,*) js(ib), jp(ib), jh(ib)
      read(10,*) w0(ib), aH2O_a(ib), num_a(ib), dp_wet_a(ib), mass_dry_a(ib)
      read(10,*) kg(1:5,ib), kel(1:5,ib)
      read(10,*) xs, xl
      read(10,*) aer0(:,jtotal,ib)
      read(10,*) el0(:,jsolid,ib), el0(:,jliquid,ib), el0(:,jtotal,ib)
      aer0(:,jsolid,ib) = 0.0d0; aer0(:,jliquid,ib) = 0.0d0
      do iv = 1, 9
        aer0(ia(iv),jsolid,ib) = xs(iv); aer0(ia(iv),jliquid,ib) = xl(iv)
      end do
    end do
    ncorecnt_aer = step; kclm_aer = lev
    do itag = 1, 2
      esm_sec_noedit = (itag == 2)
      do ib = 1, nb
        call set_section(ib)
      end do
      esm_at_on = .true.
      if (itag == 1) then
        call esm_at_ain(60.0d0)
        call esm_at_svt(1, 0.0d0)
      endif
      esm_at_sub = itag
      do ib = 1, nb
        call set_section(ib)
        gas(1:5) = g0
        phi_nh4no3_s = 0.0; phi_nh4cl_s = 0.0
        ieqblm_ASTEM = mYES
        idry_case3a(ib) = mNO
        do iv = 1, ngas_ioa
          sfc_a(iv) = gas(iv)
          df_gas_s(iv,ib) = 0.0; df_gas_l(iv,ib) = 0.0
          flux_s(iv,ib) = 0.0;   flux_l(iv,ib) = 0.0
          Heff(iv,ib) = 0.0;     volatile_s(iv,ib) = 0.0
          phi_volatile_s(iv,ib) = 0.0; phi_volatile_l(iv,ib) = 0.0
          integrate(iv,jsolid,ib) = mNO; integrate(iv,jliquid,ib) = mNO
        end do
        if (jaerosolstate(ib) .eq. all_solid) then
          jphase(ib) = jsolid
          call esm_at_fxi(ib, 1); call ASTEM_flux_dry(ib); call esm_at_fxo(ib, 1)
        elseif (jaerosolstate(ib) .eq. all_liquid) then
          jphase(ib) = jliquid
          call esm_at_fxi(ib, 2); call ASTEM_flux_wet(ib); call esm_at_fxo(ib, 2)
        elseif (jaerosolstate(ib) .eq. mixed) then
          if (electrolyte(jnh4no3,jsolid,ib) .gt. 0.0 .or. &
              electrolyte(jnh4cl, jsolid,ib) .gt. 0.0) then
            call esm_at_fxi(ib, 3); call ASTEM_flux_mix(ib); call esm_at_fxo(ib, 3)
          else
            jphase(ib) = jliquid
            call esm_at_fxi(ib, 4); call ASTEM_flux_wet(ib); call esm_at_fxo(ib, 4)
          endif
        endif
      end do
      esm_at_on = .false.
    end do
  end do
  close(10)
  close(esm_at_unit)
contains
  subroutine set_section(ib)
    integer ib
    jaerosolstate(ib) = js(ib); jphase(ib) = jp(ib); jhyst_leg(ib) = jh(ib)
    water_a(ib) = w0(ib)
    aer(:,:,ib) = aer0(:,:,ib)
    electrolyte(:,:,ib) = el0(:,:,ib)
  end subroutine set_section
end program astem_section_driver
