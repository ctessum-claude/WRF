! coagsolv_rate_driver: replay MOSAIC's coagsolv, alone, on the inputs the
! instrumented model recorded at its entry (esm dump scheme 'mosaic', fork
! f311091e: coag_dpdry, coag_dpwet, coag_denswet, coag_distrib, coag_t_k,
! coag_p_atm, coag_deltat), and write what it forms: the Brownian kernel, the
! dense bin-splitting tables and the instantaneous number and composition rates.
!   coagsolv_rate_driver <input.flat> <output.json>
! Built real64 by default (-fdefault-real-8 promotes module_mosaic_coag), PREC=
! for real32.  The real32 build reproduces the in-model rate_coag_* exactly,
! which validates the recorded inputs and this driver; the real64 build is the
! reference, because coagsolv forms the high share of a split product as
! 1 - finklow, which in real32 leaves only a few digits of the small shares a
! tiny partner contributes (EqWeFiC FORTRAN_BUGS N117).
!
! This calls coagsolv directly, not mosaic_coag_1clm: the number-conforming
! clamp and the diameter/density preparation upstream of it are already applied
! in the recorded inputs, which also keeps this driver clear of the unresolved
! real64 behaviour of kernels/mosaic_subproc_driver (FORTRAN_BUGS N87).  fink is
! written flattened to rank 3, (i, j, k + nbin*(lev-1)), because module_esm_dump
! has no real64 rank-4 writer.
program coagsolv_rate_driver
  use esm_flat_input
  use module_esm_dump
  use module_esm_mosaic_rates
  use module_mosaic_coag, only: coagsolv
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: nbin, ncomp, nk, k, iok, j
  real :: deltat
  real, allocatable :: dpdry(:,:), dpwet(:,:), denswet(:,:), distrib(:,:,:), tk(:), patm(:)
  real, allocatable :: num(:), vol(:,:), fink3(:,:,:)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  nbin = flat_i0('rate_nbin'); ncomp = flat_i0('rate_ncomp_coag')
  nk = flat_i0('kclm_calcend') - flat_i0('kclm_calcbgn') + 1
  deltat = flat_r0('coag_deltat')
  allocate(dpdry(nbin,nk), dpwet(nbin,nk), denswet(nbin,nk), distrib(nbin,ncomp+1,nk), tk(nk), patm(nk))
  allocate(num(nbin), vol(nbin,ncomp), fink3(nbin,nbin,nbin*nk))
  call flat_r2('coag_dpdry', dpdry); call flat_r2('coag_dpwet', dpwet)
  call flat_r2('coag_denswet', denswet); call flat_r3('coag_distrib', distrib)
  call flat_r1('coag_t_k', tk); call flat_r1('coag_p_atm', patm)

  call esm_rates_reset()
  esm_rate_on = .true.
  do k = 1, nk
    if (maxval(distrib(:,1,k)) <= 0.0) cycle     ! mosaic_coag_1clm skips no level; a zero column has zero rates
    esm_rate_k = k
    num = distrib(:,1,k)
    vol = distrib(:,2:ncomp+1,k)
    call coagsolv(nbin, nbin, ncomp, ncomp, tk(k), patm(k), deltat, 1, 6, -1, iok, &
                  dpdry(:,k), dpwet(:,k), denswet(:,k), num, vol)
    if (iok < 0) stop 'coagsolv failed'
  end do
  do k = 1, nk
    do j = 1, nbin
      fink3(:,:,j + nbin*(k-1)) = esm_coag_fink(1:nbin,1:nbin,j,k)
    end do
  end do

  call esm_dump_open_file('coagsolv_rate', trim(outpath))
  call esm_dump_var('rate_coag_num', esm_rate_coag_num(1:nbin,1:nk))
  call esm_dump_var('rate_coag_vol', esm_rate_coag_vol(1:nbin,1:ncomp,1:nk))
  call esm_dump_var('coag_cbr', esm_coag_cbr(1:nbin,1:nbin,1:nk))
  call esm_dump_var('coag_floss', esm_coag_floss(1:nbin,1:nbin,1:nk))
  call esm_dump_var('coag_fink_flat', fink3)
  call esm_dump_close()
end program coagsolv_rate_driver
