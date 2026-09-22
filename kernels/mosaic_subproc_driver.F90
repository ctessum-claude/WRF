! mosaic_subproc_driver: replay the MOSAIC nucleation and coagulation stages of
! one dumped column (esm dump scheme 'mosaic', written by
! chem/module_mosaic_driver.F) through WRF's own mosaic_newnuc_1clm and
! mosaic_coag_1clm, and write the resulting column state as an esm dump.
!   mosaic_subproc_driver <input.flat> <output.json>   env ESM_DT overrides dtchem
! Built real64 by default (-fdefault-real-8 promotes the WRF-native modules),
! PREC= for a real32 build.  NOTE (2026-09-22): the real32 build reproduces WRF
! bit for bit, which validates both the dump and this driver, but the real64
! build returns zero nucleation and a different coagulation state even though
! every module variable is verified correct at the call site; see
! data/eqwefic/notes/chem170_scm_reference_run.md.  Use the real32 build until
! that is resolved.  The two stages are operator-split increments over
! dtchem, and in the in-model real32 dumps the coagulation increment is
! cancellation-limited (EqWeFiC FORTRAN_BUGS N87), which is what this driver is
! for.
!
! The dump carries the section geometry, the composition pointers and the
! post-growth column diagnostics (aqvoldry_sub, aqmassdry_sub, adrydens_sub)
! that the two stages read from the data modules, so none of WRF's
! initialisation is needed.  Stage order matches the model: the input state is
! rsub AFTER aerchemistry (rsub_therm), nucleation runs first, then coagulation.
program mosaic_subproc_driver
  use esm_flat_input
  use module_esm_dump
  use module_data_mosaic_asect
  use module_data_mosaic_other
  use module_mosaic_newnuc, only: mosaic_newnuc_1clm
  use module_mosaic_coag, only: mosaic_coag_1clm
  implicit none
  character(len=1024) :: inpath, outpath, dtstr
  integer :: n, l, k, ll, nsize, ncomp, kbgn, kend, nk, stat, istat, ktau, ktauc
  real :: dtchem
  integer, allocatable :: mptr(:), nptr(:), wptr(:), hptr(:), lso4(:), lnh4(:)
  real, allocatable :: r1(:), r2(:,:), dens(:), mw(:)
  real, allocatable :: rsub_therm(:,:), rsub0_in(:,:), rsub_newnuc(:,:), rsub_coag(:,:)
  real, allocatable :: rsub0_buf(:,:,:)   ! full module extents; too big for the stack

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  ltot2 = flat_i0('ltot2'); kbgn = flat_i0('kclm_calcbgn'); kend = flat_i0('kclm_calcend')
  ktemp = flat_i0('ktemp'); kh2o = flat_i0('kh2o')
  ! the gas rows nucleation reads; dumped since the hook was extended, and
  ! overridable by ESM_KH2SO4 / ESM_KNH3 / ESM_KSO2 for older dumps
  kh2so4 = idx_of('kh2so4', 'ESM_KH2SO4')
  knh3   = idx_of('knh3',   'ESM_KNH3')
  kso2   = idx_of('kso2',   'ESM_KSO2')
  nk = kend - kbgn + 1
  nsize = flat_i0('nsize_aer'); ncomp = flat_i0('ncomp_aer')
  ktau = flat_i0('ktau'); ktauc = flat_i0('ktauc')
  dtchem = flat_r0('dtchem')
  call get_environment_variable('ESM_DT', dtstr, status=stat)
  if (stat == 0 .and. len_trim(dtstr) > 0) read(dtstr, *) dtchem

  ntype_aer = 1; nphase_aer = flat_i0('nphase_aer'); ai_phase = flat_i0('ai_phase')
  nsize_aer(1) = nsize; ncomp_aer(1) = ncomp; ncomp_plustracer_aer(1) = flat_i0('ncomp_plustracer_aer')

  allocate(mptr(ncomp*nsize), nptr(nsize), wptr(nsize), hptr(nsize), lso4(nsize), lnh4(nsize))
  allocate(dens(ncomp), mw(ncomp), r1(nsize), r2(nsize,nk))
  allocate(rsub_therm(ltot2,nk), rsub0_in(ltot2,nk), rsub_newnuc(ltot2,nk), rsub_coag(ltot2,nk))
  call flat_i1('massptr_aer', mptr); call flat_i1('numptr_aer', nptr)
  call flat_i1('waterptr_aer', wptr); call flat_i1('hyswptr_aer', hptr)
  call flat_i1('lptr_so4_aer', lso4); call flat_i1('lptr_nh4_aer', lnh4)
  call flat_r1('dens_aer', dens); call flat_r1('mw_aer', mw)
  do ll = 1, ncomp
    dens_aer(ll,1) = dens(ll); mw_aer(ll,1) = mw(ll)
  end do
  do n = 1, nsize
    numptr_aer(n,1,ai_phase) = nptr(n); waterptr_aer(n,1) = wptr(n); hyswptr_aer(n,1) = hptr(n)
    lptr_so4_aer(n,1,ai_phase) = lso4(n); lptr_nh4_aer(n,1,ai_phase) = lnh4(n)
    do ll = 1, ncomp
      massptr_aer(ll,n,1,ai_phase) = mptr(ll + (n-1)*ncomp)
    end do
  end do
  call flat_r1('dcen_sect', r1);     dcen_sect(1:nsize,1) = r1
  call flat_r1('dlo_sect', r1);      dlo_sect(1:nsize,1) = r1
  call flat_r1('dhi_sect', r1);      dhi_sect(1:nsize,1) = r1
  call flat_r1('volumcen_sect', r1); volumcen_sect(1:nsize,1) = r1
  call flat_r1('volumlo_sect', r1);  volumlo_sect(1:nsize,1) = r1
  call flat_r1('volumhi_sect', r1);  volumhi_sect(1:nsize,1) = r1
  call flat_r2('aqvoldry_sub', r2);  aqvoldry_sub(1:nsize,1,kbgn:kend,1) = r2
  call flat_r2('aqmassdry_sub', r2); aqmassdry_sub(1:nsize,1,kbgn:kend,1) = r2
  call flat_r2('adrydens_sub', r2);  adrydens_sub(1:nsize,1,kbgn:kend,1) = r2

  call flat_r2('rsub_therm', rsub_therm); call flat_r2('rsub0', rsub0_in)
  rsub = 0.0
  do k = 1, nk
    do l = 1, ltot2
      rsub(l,kbgn+k-1,1) = rsub_therm(l,k)
    end do
  end do
  cairclm = 0.0; ptotclm = 0.0; relhumclm = 0.0; afracsubarea = 0.0
  block
    real, allocatable :: c(:)
    allocate(c(nk))
    call flat_r1('cairclm', c);   cairclm(kbgn:kend) = c
    call flat_r1('ptotclm', c);   ptotclm(kbgn:kend) = c
    call flat_r1('relhumclm', c); relhumclm(kbgn:kend) = c
    call flat_r1('afracsubarea', c); afracsubarea(kbgn:kend,1) = c
  end block
  nsubareas = 1; lunerr = -1; lunout = -1; ncorecnt = ktau - 1

  allocate(rsub0_buf(l2maxd, kmaxd, nsubareamaxd))
  rsub0_buf = 0.0
  do k = 1, nk
    do l = 1, ltot2
      rsub0_buf(l, kbgn+k-1, 1) = rsub0_in(l, k)
    end do
  end do

  istat = 0
  call mosaic_newnuc_1clm(istat, 1, 1, kbgn, kend, flat_i0('idiagbb_dum'), dtchem, dtchem, &
       rsub0_buf, 1, ktau, ktauc, 1, 1, 1, 1, 1, nk)
  do k = 1, nk
    rsub_newnuc(:,k) = rsub(1:ltot2,kbgn+k-1,1)
  end do
  call mosaic_coag_1clm(istat, 1, 1, kbgn, kend, flat_i0('idiagbb_dum'), dtchem, dtchem, &
       1, ktau, ktauc, 1, 1, 1, 1, 1, nk)
  do k = 1, nk
    rsub_coag(:,k) = rsub(1:ltot2,kbgn+k-1,1)
  end do

  call esm_dump_open_file('mosaic', trim(outpath))
  call esm_dump_var('real_bytes', storage_size(dtchem)/8)
  call esm_dump_var('dtchem', dtchem); call esm_dump_var('istat', istat)
  call esm_dump_var('rsub_newnuc_replay', rsub_newnuc)
  call esm_dump_var('rsub_coag_replay', rsub_coag)
  call esm_dump_close()

contains

  ! the rsub row index of a gas nucleation reads: from the dump when the hook
  ! carried it, else from an environment variable (older dumps)
  integer function idx_of(name, envvar) result(k)
    character(len=*), intent(in) :: name, envvar
    character(len=64) :: val
    integer :: st
    if (flat_has(name)) then
      k = flat_i0(name); return
    end if
    call get_environment_variable(envvar, val, status=st)
    if (st /= 0 .or. len_trim(val) == 0) then
      write(0,'(A)') 'mosaic_subproc_driver: need '//trim(name)//' in the dump or '//trim(envvar)
      error stop 2
    end if
    read(val,*) k
  end function idx_of

end program mosaic_subproc_driver
