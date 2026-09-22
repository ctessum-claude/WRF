! mosaic_drydep_driver: replay one dumped column (esm dump scheme
! 'mosaic_drydep', written by chem/module_mosaic_drydep.F) through
! mosaic_drydep_1clm / aerosol_depvel_2, extracted verbatim from that file at
! build time, and write the per-bin deposition velocities as an esm dump.
!   mosaic_drydep_driver <input.flat> <output.json>
! Built real64 by default (-fdefault-real-8 promotes the WRF-native modules),
! PREC= for a real32 build, so the two can be compared.
! The routines read their section geometry and composition pointers from
! module_data_mosaic_asect and the column state from module_data_mosaic_other;
! the dump carries both, so nothing of WRF's initialisation is needed.
program mosaic_drydep_driver
  use esm_flat_input
  use module_esm_dump
  use module_data_mosaic_asect
  use module_data_mosaic_other
  use mosaic_drydep_kernel, only: mosaic_drydep_1clm
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: n, l, nsize, ncomp, ll
  real :: ustar_in, depresist_a_in, cair
  real :: vdep(maxd_asize,maxd_atype,maxd_aphase)
  real, allocatable :: rsub_in(:), dcen(:), dlo(:), dhi(:), vlo(:), vhi(:), dens(:), mw(:)
  integer, allocatable :: mptr(:), nptr(:), wptr(:)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  ltot2 = flat_i0('ltot2'); ktemp = flat_i0('ktemp')
  nsize = flat_i0('nsize_aer'); ncomp = flat_i0('ncomp_aer')
  allocate(rsub_in(ltot2), dcen(nsize), dlo(nsize), dhi(nsize), vlo(nsize), vhi(nsize))
  allocate(dens(ncomp), mw(ncomp), mptr(ncomp*nsize), nptr(nsize), wptr(nsize))
  call flat_r1('rsub', rsub_in)
  call flat_r1('dcen_sect', dcen); call flat_r1('dlo_sect', dlo); call flat_r1('dhi_sect', dhi)
  call flat_r1('volumlo_sect', vlo); call flat_r1('volumhi_sect', vhi)
  call flat_r1('dens_aer', dens); call flat_r1('mw_aer', mw)
  call flat_i1('massptr_aer', mptr); call flat_i1('numptr_aer', nptr); call flat_i1('waterptr_aer', wptr)
  ustar_in = flat_r0('ust'); depresist_a_in = flat_r0('aer_res'); cair = flat_r0('cairclm')
  ! mw_water_aer and dens_water_aer are parameters in module_data_mosaic_asect; the dump carries them only as a check

  ! module state the kernel reads
  ntype_aer = 1; nphase_aer = 1; ai_phase = 1
  nsize_aer(1) = nsize; ncomp_aer(1) = ncomp; ncomp_plustracer_aer(1) = ncomp
  do n = 1, nsize
    dcen_sect(n,1) = dcen(n); dlo_sect(n,1) = dlo(n); dhi_sect(n,1) = dhi(n)
    volumlo_sect(n,1) = vlo(n); volumhi_sect(n,1) = vhi(n)
    numptr_aer(n,1,ai_phase) = nptr(n); waterptr_aer(n,1) = wptr(n)
    do ll = 1, ncomp
      massptr_aer(ll,n,1,ai_phase) = mptr(ll + (n-1)*ncomp)
    end do
  end do
  do ll = 1, ncomp
    dens_aer(ll,1) = dens(ll); mw_aer(ll,1) = mw(ll)
  end do
  rsub = 0.0
  do l = 1, ltot2
    rsub(l,1,1) = rsub_in(l)
  end do
  cairclm(1) = cair

  vdep = 0.0
  call mosaic_drydep_1clm(0, 1, 1, ustar_in, depresist_a_in, vdep)

  call esm_dump_open_file('mosaic_drydep', trim(outpath))
  call esm_dump_var('real_bytes', storage_size(ustar_in)/8)
  call esm_dump_var('nsize_aer', nsize); call esm_dump_var('ncomp_aer', ncomp)
  call esm_dump_var('vdep_aer', vdep(1:nsize,1,ai_phase))
  call esm_dump_close()
end program mosaic_drydep_driver
