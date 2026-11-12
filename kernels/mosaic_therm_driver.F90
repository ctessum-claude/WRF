! mosaic_therm_driver: run MOSAIC's thermodynamics (chem/module_mosaic_therm.F)
! on hand-specified bin compositions, through WRF's OWN compiled objects
! (main/libwrflib.a, built with the model's flags: default REAL is real32 and the
! thermodynamics declares real(kind=8)), so a result is what the model itself
! would compute for that bin.  It exists to reach branches the chem170 column
! never visits (EqWeFiC FORTRAN_BUGS / components/aerosol/mosaic coverage notes).
!
!   mosaic_therm_driver <cases.txt> <out.json>
!
! cases.txt, whitespace separated:
!   ncase
!   then per case:
!     mode t_k p_atm relhum num_a water_a water_a_hyst water_a_up jhyst_leg
!     nnz
!     nnz lines "iaer value"   -- aer(iaer, jtotal, 1), nmol/m3(air)
!   mode 1: compute_activities on the fully deliquesced bin: ah2o_a = relhum,
!           conform_electrolytes(jtotal), do_full_deliquescence, compute_activities
!   mode 3: as mode 1, but after conform_electrolytes a second list
!           "nextra" + "iaer value" lines is ADDED to aer(:, jtotal) -- an ion
!           excess the electrolytes do not carry, as ASTEM's flux updates leave
!           between conform calls (reaches compute_activities' a_c /= 0)
!   mode 2: the model's own phase-state step: conform_electrolytes(jtotal),
!           aerosol_phase_state (Kelvin iteration around MESA)
! relhum and t_k are read as default REAL (the model's relhumclm and rsub are
! real32) and rh_pc = 100.*relhum is formed in default REAL, as in mosaic().
! Each case is recorded through the model's own instrumentation hooks
! (esm_th_* for compute_activities, esm_ps_* for aerosol_phase_state) with the
! case number as the "level", and dumped under the same names the model uses
! (th_*_c1/_c2, ps_*_c1/_c2), with nbin = 1.
program mosaic_therm_driver
  use module_esm_dump
  use module_data_mosaic_therm
  use module_esm_mosaic_rates
  use module_mosaic_therm
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: ncase, ic, mode, nnz, i, ia, jh, vbs(1), islot
  real :: tk4, rh4, patm4
  real(kind=8) :: numa, wa, wh, wu, v, xt
  character(len=3) :: sfx

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  open(10, file=trim(inpath), status='old')
  read(10,*) ncase
  if (ncase > esm_maxk) stop 'too many cases'
  vbs(1) = 0
  nbin_a = 1
  call esm_rates_reset()
  esm_nbin = 1
  call load_mosaic_parameters
  do ic = 1, ncase
    read(10,*) mode, tk4, patm4, rh4, numa, wa, wh, wu, jh
    read(10,*) nnz
    call initialize_mosaic_variables
    t_k = tk4
    p_atm = patm4
    rh_pc = 100.*rh4
    ah2o = rh4
    call update_thermodynamic_constants(vbs)
    num_a(1) = numa; water_a(1) = wa; water_a_hyst(1) = wh; water_a_up(1) = wu
    jhyst_leg(1) = jh
    do i = 1, nnz
      read(10,*) ia, v
      aer(ia,jtotal,1) = v
    enddo
    jaerosolstate(1) = 0
    istat_mosaic_fe1 = 1
    esm_rate_on = .true.; esm_rate_k = ic
    call conform_electrolytes(jtotal,1,xt)
    if (mode == 3) then
      read(10,*) nnz
      do i = 1, nnz
        read(10,*) ia, v
        aer(ia,jtotal,1) = aer(ia,jtotal,1) + v
      enddo
    endif
    if (mode == 1 .or. mode == 3) then
      ah2o_a(1) = rh4
      do i = 1, nelectrolyte
        molality0(i) = bin_molality(i,1)
      enddo
      call do_full_deliquescence(1)
      call compute_activities(1)
    else
      call aerosol_phase_state(1)
    endif
    esm_rate_on = .false.
  enddo
  close(10)

  call esm_dump_open_file('mosaic_therm', trim(outpath))
  call esm_dump_var('ncase', ncase)
  call esm_dump_var('th_ncall', esm_th_ncall(1:1,1:ncase))
  call esm_dump_var('th_ah2o', esm_th_ah2o(1:ncase))
  call esm_dump_var('th_t_k', esm_th_tk(1:ncase))
  call esm_dump_var('th_keq_ll', esm_th_keq_ll(:,1:ncase))
  call esm_dump_var('th_log_gamz', esm_th_log_gamz(:,:,1:ncase))
  call esm_dump_var('th_iaer', esm_th_iaer)
  call esm_dump_var('ps_ncall', esm_ps_ncall(1:1,1:ncase))
  call esm_dump_var('ps_rh_pc', esm_ps_rh_pc(1:ncase))
  call esm_dump_var('ps_t_k', esm_ps_tk(1:ncase))
  call esm_dump_var('ps_sigma_water', esm_ps_sigma_water)
  call esm_dump_var('ps_partial_molar_vol', esm_ps_pmv)
  call esm_dump_var('ps_mdrh_t', esm_ps_mdrh_t)
  call esm_dump_var('ps_ptol_mol_astem', esm_ps_ptol)
  ! MESA's mutual-deliquescence tables (load_mosaic_parameters), for the step-3 path
  call esm_dump_var('jsalt_index', jsalt_index(1:nsalt))
  call esm_dump_var('jsulf_poor', jsulf_poor)
  call esm_dump_var('jsulf_rich', jsulf_rich)
  call esm_dump_var('d_mdrh', d_mdrh)
  call esm_dump_var('mw_electrolyte', mw_electrolyte(1:nelectrolyte))
  call esm_dump_var('mdrh_t_last', mdrh_t(1:63))
  do islot = 1, 2
    write(sfx,'(a,i1)') '_c', islot
    call esm_dump_var('th_elec'//sfx, esm_th_elec(:,1:1,1:ncase,islot))
    call esm_dump_var('th_aerl'//sfx, esm_th_aerl(:,1:1,1:ncase,islot))
    call esm_dump_var('th_aert'//sfx, esm_th_aert(:,1:1,1:ncase,islot))
    call esm_dump_var('th_ah2o_a'//sfx, esm_th_ah2o_a(1:1,1:ncase,islot))
    call esm_dump_var('th_water'//sfx, esm_th_water(1:1,1:ncase,islot))
    call esm_dump_var('th_xt'//sfx, esm_th_xt(1:1,1:ncase,islot))
    call esm_dump_var('th_gam_ratio'//sfx, esm_th_gam_ratio(1:1,1:ncase,islot))
    call esm_dump_var('th_mc'//sfx, esm_th_mc(:,1:1,1:ncase,islot))
    call esm_dump_var('th_ma'//sfx, esm_th_ma(:,1:1,1:ncase,islot))
    call esm_dump_var('th_gam'//sfx, esm_th_gam(:,1:1,1:ncase,islot))
    call esm_dump_var('th_act'//sfx, esm_th_act(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_aert_in'//sfx, esm_ps_aert_in(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_elet_in'//sfx, esm_ps_elet_in(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_aers'//sfx, esm_ps_aers(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_aerl'//sfx, esm_ps_aerl(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_eles'//sfx, esm_ps_eles(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_elel'//sfx, esm_ps_elel(:,1:1,1:ncase,islot))
    call esm_dump_var('ps_kel'//sfx, esm_ps_kel(1:8,1:1,1:ncase,islot))
    call esm_dump_var('ps_sc'//sfx, esm_ps_sc(:,1:1,1:ncase,islot))
  end do
  call esm_dump_close()
end program mosaic_therm_driver
