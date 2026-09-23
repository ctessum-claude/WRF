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

!   nucleation, dt -> 0 limit, in the units of rsub (mol/mol for gases and
!   aerosol species, particles/mol-air for number) per second
      real, save :: esm_rate_newnuc(6, esm_maxk)
!     rows: 1 h2so4, 2 nh3, 3 so4_a(bin), 4 nh4_a(bin), 5 num_a(bin),
!           6 the bin index the new particles go into (not a rate)

!   nucleation as WRF actually applies it, increment/dtnuc, same rows/units
      real, save :: esm_rate_newnuc_fd(6, esm_maxk)

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

!   shapes actually used, so a consumer does not have to guess
      integer, save :: esm_nbin = 0, esm_ncomp_coag = 0, esm_ngas = 0
!   MOSAIC's gas indices, dumped with the rates
      integer, save :: esm_ih2so4 = 0, esm_ihno3 = 0, esm_ihcl = 0, &
                       esm_inh3 = 0, esm_imsa = 0

      contains

      subroutine esm_rates_reset()
      esm_rate_newnuc      = 0.0
      esm_rate_newnuc_fd   = 0.0
      esm_rate_coag_num    = 0.0
      esm_rate_coag_vol    = 0.0
      esm_rate_cond_nonvol = 0.0
      esm_rate_cond_semivol = 0.0
      end subroutine esm_rates_reset

      end module module_esm_mosaic_rates
