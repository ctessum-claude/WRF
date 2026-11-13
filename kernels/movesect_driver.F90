! movesect_driver: replay MOSAIC's move_sections (chem/module_mosaic_movesect.F,
! msectional = 20, the linear-discrete / mass-number-advection algorithm) on
! the state recorded at its OWN entry for every level of one dumped column
! (esm dump scheme 'mosaic', fields ms_* written by module_mosaic_therm.F's
! hook), and write the state at its exit.
!   movesect_driver <input.flat> <output.json>
! env ESM_MS_DELTA (default 1) scales the growth the call sees: the entry state
! becomes  pre + delta*(aft - pre)  for every rsub row (pre = rsub_in, the
! aerchemistry entry state of that level, aft = ms_rsub_bef) and for the dry
! mass; the dry VOLUME is interpolated the same way and the dry density is
! their ratio.  delta = 1 is the in-model call (real32 cross-check); delta -> 0
! gives the instantaneous transfer rate  d(transfer)/d(delta) / dtchem,  which
! is what the .esm operator is tested against.  ESM_MS_FREEZE=1 holds the rsub
! rows at their pre-growth values (see below).
! Built real64 (-fdefault-real-8 promotes the WRF-native modules).
program movesect_driver
  use esm_flat_input
  use module_esm_dump
  use module_data_mosaic_asect
  use module_data_mosaic_other
  use module_mosaic_movesect, only: move_sections
  use module_esm_mosaic_rates, only: esm_rate_on, esm_rate_k, esm_ms_xfer
  implicit none
  character(len=1024) :: inpath, outpath, dstr
  integer :: n, l, k, ll, nsize, ncomp, kbgn, kend, nk, stat, lms
  logical :: freeze
  real :: delta, vpre, vaft, vd, md
  integer, allocatable :: mptr(:), nptr(:), wptr(:), hptr(:)
  real, allocatable :: r1(:), dens(:), mw(:)
  real, allocatable :: rin(:,:), rbef(:,:), grow(:,:,:)
  real, allocatable :: o_bef(:,:), o_aft(:,:), o_grow(:,:,:), o_diag(:,:,:)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  delta = 1.0
  call get_environment_variable('ESM_MS_DELTA', dstr, status=stat)
  if (stat == 0 .and. len_trim(dstr) > 0) read(dstr, *) delta
  ! ESM_MS_FREEZE=1 keeps every rsub row at its PRE-growth value and scales
  ! only the growth (dry mass and volume): the transfer is then m0 x(delta)
  ! exactly, so its delta-derivative m0 x'(0) is free of the m(delta) x(delta)
  ! cross term that dominates for a component that grows by orders of
  ! magnitude in one step (e.g. nitrate at the first step).
  call get_environment_variable('ESM_MS_FREEZE', dstr, status=stat)
  freeze = (stat == 0 .and. trim(dstr) == '1')

  lms = flat_i0('ms_ltot'); ltot2 = flat_i0('ltot2')
  kbgn = flat_i0('kclm_calcbgn'); kend = flat_i0('kclm_calcend'); nk = kend - kbgn + 1
  nsize = flat_i0('nsize_aer'); ncomp = flat_i0('ncomp_aer')
  msectional = flat_i0('msectional')
  ntype_aer = 1; nphase_aer = flat_i0('nphase_aer'); ai_phase = flat_i0('ai_phase')
  cw_phase = 2
  nsize_aer(1) = nsize; ncomp_aer(1) = ncomp; ncomp_plustracer_aer(1) = flat_i0('ncomp_plustracer_aer')
  nsubareas = 1; lunerr = -1; lunout = -1

  allocate(mptr(ncomp*nsize), nptr(nsize), wptr(nsize), hptr(nsize), dens(ncomp), mw(ncomp), r1(nsize))
  call flat_i1('massptr_aer', mptr); call flat_i1('numptr_aer', nptr)
  call flat_i1('waterptr_aer', wptr); call flat_i1('hyswptr_aer', hptr)
  call flat_r1('dens_aer', dens); call flat_r1('mw_aer', mw)
  do ll = 1, ncomp
    dens_aer(ll,1) = dens(ll); mw_aer(ll,1) = mw(ll)
  end do
  do n = 1, nsize
    numptr_aer(n,1,ai_phase) = nptr(n); waterptr_aer(n,1) = wptr(n); hyswptr_aer(n,1) = hptr(n)
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

  allocate(rin(ltot2,nk), rbef(lms,nk), grow(4,nsize,nk))
  allocate(o_bef(lms,nk), o_aft(lms,nk), o_grow(4,nsize,nk), o_diag(6,nsize,nk))
  call flat_r2('rsub_in', rin); call flat_r2('ms_rsub_bef', rbef); call flat_r3('ms_grow', grow)

  esm_ms_xfer = 0.0
  do k = 1, nk
    rsub(:,kbgn+k-1,1) = 0.0
    do l = 1, lms
      if (freeze) then
        rsub(l,kbgn+k-1,1) = rin(l,k)
      else
        rsub(l,kbgn+k-1,1) = rin(l,k) + delta*(rbef(l,k) - rin(l,k))
      end if
    end do
    do n = 1, nsize
      drymass_pregrow(n,1) = grow(1,n,k); drydens_pregrow(n,1) = grow(2,n,k)
      md = grow(1,n,k) + delta*(grow(3,n,k) - grow(1,n,k))
      if (grow(2,n,k) > 0.0 .and. grow(4,n,k) > 0.0) then
        vpre = grow(1,n,k)/grow(2,n,k); vaft = grow(3,n,k)/grow(4,n,k)
        vd = vpre + delta*(vaft - vpre)
        drydens_aftgrow(n,1) = md/vd
      else
        drydens_aftgrow(n,1) = grow(4,n,k)
      end if
      drymass_aftgrow(n,1) = md
      o_grow(1,n,k) = drymass_pregrow(n,1); o_grow(2,n,k) = drydens_pregrow(n,1)
      o_grow(3,n,k) = drymass_aftgrow(n,1); o_grow(4,n,k) = drydens_aftgrow(n,1)
    end do
    o_bef(:,k) = rsub(1:lms,kbgn+k-1,1)
    esm_rate_on = .true.; esm_rate_k = kbgn+k-1
    call move_sections(1, 1, 1, kbgn+k-1, 1)
    esm_rate_on = .false.
    o_aft(:,k) = rsub(1:lms,kbgn+k-1,1)
    do n = 1, nsize
      o_diag(1,n,k) = adrydens_sub(n,1,kbgn+k-1,1); o_diag(2,n,k) = awetdens_sub(n,1,kbgn+k-1,1)
      o_diag(3,n,k) = admeandry_sub(n,1,kbgn+k-1,1); o_diag(4,n,k) = admeanwet_sub(n,1,kbgn+k-1,1)
      o_diag(5,n,k) = aqvoldry_sub(n,1,kbgn+k-1,1); o_diag(6,n,k) = aqmassdry_sub(n,1,kbgn+k-1,1)
    end do
  end do

  call esm_dump_open_file('mosaic', trim(outpath))
  call esm_dump_var('real_bytes', storage_size(delta)/8)
  call esm_dump_var('delta', delta)
  call esm_dump_var('freeze', merge(1, 0, freeze))
  call esm_dump_var('ms_rsub_bef_replay', o_bef)
  call esm_dump_var('ms_rsub_aft_replay', o_aft)
  call esm_dump_var('ms_grow_replay', o_grow)
  call esm_dump_var('ms_diag_replay', o_diag)
  call esm_dump_var('ms_xfer_replay', esm_ms_xfer(:,1:nsize,kbgn:kend))
  call esm_dump_close()
end program movesect_driver
