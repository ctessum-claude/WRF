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
      real, save :: esm_coag_denswet(esm_maxbin, esm_maxk)
      real, save :: esm_coag_tk(esm_maxk), esm_coag_patm(esm_maxk)

!   shapes actually used, so a consumer does not have to guess
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
      esm_coag_denswet = 0.0
      esm_coag_tk = 0.0
      esm_coag_patm = 0.0
      end subroutine esm_rates_reset

      end module module_esm_mosaic_rates
