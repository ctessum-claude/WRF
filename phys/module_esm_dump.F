!=======================================================================
! module_esm_dump: structured (JSON) dumps of physics-kernel inputs and
! outputs for the EarthSciML EqWeFiC project (branch earthsciml-instrumented).
!
! Usage in a physics wrapper, around a kernel call:
!   use module_esm_dump
!   ...
!   dumping = esm_dump_want('ysu')          ! counts calls per scheme
!   if (dumping) call esm_dump_open('ysu')
!   if (dumping) call esm_dump_var('tx', t3d_hv)   ! any rank 0-3 real/int/logical
!   call bl_ysu_run(...)
!   if (dumping) call esm_dump_var('ttnp', rthblten_hv)
!   if (dumping) call esm_dump_close()
!
! Sub-process diagnostics inside a kernel (guarded so the kernel stays clean
! when built without the dump module, e.g. #ifdef ESM_DUMP in physics_mmm):
!   if (esm_dump_inner) call esm_dump_var('praut_raw', praut)
! The wrapper sets esm_dump_inner = dumping (or dumping .and. i == its for
! per-column kernels) around the call; esm_dump_active() tells whether a
! dump file is currently open.
!
! Control (environment variables):
!   ESM_DUMP       comma-separated scheme names to dump, or "all" (unset = off)
!   ESM_DUMP_CALLS comma-separated 1-based call indices per scheme to dump
!                  (default "1"); every call increments the scheme's counter.
!   ESM_DUMP_DIR   output directory (default ".")
!
! Location mode, for 3-D multi-rank, multi-domain runs (all optional; each
! filter applies only when its variable is set):
!   ESM_DUMP_DOMAIN   dump only while WRF is solving this grid id
!                     (esm_dump_domain, set by solve_em)
!   ESM_DUMP_STEPS    comma-separated step numbers or ranges "a-b" of that
!                     domain's itimestep (esm_dump_step, set by solve_em);
!                     when set, ESM_DUMP_CALLS is ignored
!   ESM_DUMP_COLUMNS  comma-separated "i:j" global mass-point indices; a hook
!                     that passes its row (want(scheme, j=, its=, ite=)) dumps
!                     only when one of the columns lies in that tile row
!   ESM_DUMP_FIRE_BOX "i0:i1:j0:j1" global fire-mesh indices; a hook that
!                     passes fbox = [ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme]
!                     dumps only when its tile meets the box, and every 2-D
!                     array whose shape is the memory or the tile extent is
!                     written cut down to (box intersect tile); the cut is
!                     recorded as box_i0/box_i1/box_j0/box_j1
! With COLUMNS or FIRE_BOX set, a hook that passes no location never dumps.
! In location mode the file is
!   <dir>/esm_dump_<scheme>_d<domain>_s<step>_<loc>_<call>.json
! where <loc> is j<j>_i<its> (row hooks) or f<ifts>_<jfts> (fire tiles), so
! files from different MPI ranks never collide.
!
! One file per dumped call: <dir>/esm_dump_<scheme>_<call>.json with
!   {"scheme": "...", "call": n, "vars": {"name": {"kind": ..., "shape": [..],
!    "data": [...]}, ...}}
! Arrays are written in Fortran (column-major) order, flattened; "shape"
! gives the extents so an extractor can reshape. Reals use ES24.16E3.
!=======================================================================
module module_esm_dump
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  private
  public :: esm_dump_want, esm_dump_open, esm_dump_open_file, esm_dump_close, esm_dump_var
  public :: esm_dump_active, esm_dump_inner
  public :: esm_dump_domain, esm_dump_step

  ! Set by the host model (WRF solve_em) before any physics of a step.
  integer, save :: esm_dump_domain = 0
  integer, save :: esm_dump_step = 0

  ! Location-mode configuration, parsed once from the environment.
  logical, save :: cfg_done = .false.
  integer, save :: cfg_domain = -1
  integer, parameter :: max_items = 64
  integer, save :: n_steps = -1, step_lo(max_items), step_hi(max_items)
  integer, save :: n_cols = -1, col_i(max_items), col_j(max_items)
  logical, save :: have_fbox = .false.
  integer, save :: fbox_cfg(4)
  ! Current file's location tag and 2-D cut window.
  character(len=64), save :: cur_loc = ''
  logical, save :: win_on = .false.
  integer, save :: win_mem(4), win_tile(4), win_cut(4)

  ! esm_dump_inner: set by a wrapper (or driver) to let dump calls placed
  ! INSIDE a kernel (sub-process diagnostics) append to the open file.  The
  ! wrapper sets it .true. only for the call/column it wants, so kernels that
  ! are called once per column do not write duplicate keys.
  logical, save :: esm_dump_inner = .false.

  integer, parameter :: max_schemes = 32
  character(len=32) :: scheme_names(max_schemes) = ''
  integer :: scheme_calls(max_schemes) = 0
  integer :: nschemes = 0
  integer :: unit = -1
  logical :: first_var = .true.
  character(len=32) :: cur_scheme = ''
  integer :: cur_call = 0

  interface esm_dump_var
    module procedure dump_r8_0, dump_r8_1, dump_r8_2, dump_r8_3
    module procedure dump_r4_0, dump_r4_1, dump_r4_2, dump_r4_3
    module procedure dump_i_0, dump_i_1, dump_i_2
    module procedure dump_l_0, dump_l_1
    module procedure dump_str
  end interface

contains

  logical function env_list_has(varname, item, default) result(has)
    character(len=*), intent(in) :: varname, item, default
    character(len=1024) :: val
    integer :: stat, p, q
    call get_environment_variable(varname, val, status=stat)
    if (stat /= 0) val = default
    has = .false.
    if (len_trim(val) == 0) return
    if (trim(val) == 'all') then
      has = .true.; return
    end if
    p = 1
    do
      q = index(val(p:), ',')
      if (q == 0) then
        if (trim(adjustl(val(p:))) == trim(item)) has = .true.
        exit
      end if
      if (trim(adjustl(val(p:p+q-2))) == trim(item)) has = .true.
      p = p + q
      if (p > len_trim(val)) exit
    end do
  end function env_list_has

  ! Split a comma-separated environment variable into at most max_items items.
  subroutine env_items(varname, items, n)
    character(len=*), intent(in) :: varname
    character(len=64), intent(out) :: items(max_items)
    integer, intent(out) :: n
    character(len=4096) :: val
    integer :: stat, p, q
    n = -1
    call get_environment_variable(varname, val, status=stat)
    if (stat /= 0 .or. len_trim(val) == 0) return
    n = 0
    p = 1
    do while (p <= len_trim(val) .and. n < max_items)
      q = index(val(p:), ',')
      n = n + 1
      if (q == 0) then
        items(n) = adjustl(val(p:len_trim(val)))
        exit
      end if
      items(n) = adjustl(val(p:p+q-2))
      p = p + q
    end do
  end subroutine env_items

  ! Replace every separator character in s by a blank, so a list-directed
  ! read can take the integers out of "a-b", "i:j" or "i0:i1:j0:j1".
  function blanked(s, sep) result(r)
    character(len=*), intent(in) :: s, sep
    character(len=len(s)) :: r
    integer :: k
    r = s
    do k = 1, len(r)
      if (index(sep, r(k:k)) > 0) r(k:k) = ' '
    end do
  end function blanked

  subroutine read_config()
    character(len=64) :: items(max_items)
    character(len=64) :: s, b
    integer :: n, k, stat
    cfg_done = .true.
    call get_environment_variable('ESM_DUMP_DOMAIN', s, status=stat)
    if (stat == 0 .and. len_trim(s) > 0) read(s, *) cfg_domain
    call env_items('ESM_DUMP_STEPS', items, n)
    n_steps = n
    do k = 1, max(n, 0)
      if (index(items(k), '-') > 0) then
        b = blanked(items(k), '-')
        read(b, *) step_lo(k), step_hi(k)
      else
        read(items(k), *) step_lo(k)
        step_hi(k) = step_lo(k)
      end if
    end do
    call env_items('ESM_DUMP_COLUMNS', items, n)
    n_cols = n
    do k = 1, max(n, 0)
      b = blanked(items(k), ':')
      read(b, *) col_i(k), col_j(k)
    end do
    call get_environment_variable('ESM_DUMP_FIRE_BOX', s, status=stat)
    have_fbox = stat == 0 .and. len_trim(s) > 0
    if (have_fbox) then
      b = blanked(s, ':')
      read(b, *) fbox_cfg
    end if
  end subroutine read_config

  ! j, its, ite: the tile row of an atmosphere-grid hook.
  ! fbox: [ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme] of a fire-mesh hook.
  logical function esm_dump_want(scheme, j, its, ite, fbox) result(want)
    character(len=*), intent(in) :: scheme
    integer, intent(in), optional :: j, its, ite
    integer, intent(in), optional :: fbox(8)
    integer :: is, k
    character(len=16) :: callstr
    logical :: hit
    want = .false.
    if (.not. env_list_has('ESM_DUMP', scheme, '')) return
    is = 0
    do is = 1, nschemes
      if (trim(scheme_names(is)) == trim(scheme)) exit
    end do
    if (is > nschemes) then
      nschemes = nschemes + 1
      is = nschemes
      scheme_names(is) = scheme
    end if
    scheme_calls(is) = scheme_calls(is) + 1
    if (.not. cfg_done) call read_config()
    if (cfg_domain >= 0 .and. cfg_domain /= esm_dump_domain) return
    if (n_steps >= 0) then
      hit = .false.
      do k = 1, n_steps
        if (esm_dump_step >= step_lo(k) .and. esm_dump_step <= step_hi(k)) hit = .true.
      end do
      if (.not. hit) return
    else
      write(callstr, '(I0)') scheme_calls(is)
      if (.not. env_list_has('ESM_DUMP_CALLS', trim(callstr), '1')) return
    end if
    cur_loc = ''
    win_on = .false.
    if (n_cols >= 0 .or. have_fbox) then
      if (present(j) .and. present(its) .and. present(ite) .and. n_cols >= 0) then
        hit = .false.
        do k = 1, n_cols
          if (col_j(k) == j .and. col_i(k) >= its .and. col_i(k) <= ite) hit = .true.
        end do
        if (.not. hit) return
        write(cur_loc, '("j",I0,"_i",I0)') j, its
      else if (present(fbox) .and. have_fbox) then
        win_cut(1) = max(fbox_cfg(1), fbox(1)); win_cut(2) = min(fbox_cfg(2), fbox(2))
        win_cut(3) = max(fbox_cfg(3), fbox(3)); win_cut(4) = min(fbox_cfg(4), fbox(4))
        if (win_cut(1) > win_cut(2) .or. win_cut(3) > win_cut(4)) return
        win_tile = fbox(1:4)
        win_mem = fbox(5:8)
        win_on = .true.
        write(cur_loc, '("f",I0,"_",I0)') fbox(1), fbox(3)
      else
        return
      end if
    end if
    want = .true.
    cur_scheme = scheme
    cur_call = scheme_calls(is)
  end function esm_dump_want

  subroutine esm_dump_open(scheme)
    character(len=*), intent(in) :: scheme
    character(len=1024) :: dir, fname
    integer :: stat
    logical :: cut
    call get_environment_variable('ESM_DUMP_DIR', dir, status=stat)
    if (stat /= 0 .or. len_trim(dir) == 0) dir = '.'
    if (len_trim(cur_loc) > 0) then
      write(fname, '(A,"/esm_dump_",A,"_d",I0,"_s",I0,"_",A,"_",I0,".json")') trim(dir), trim(scheme), &
          esm_dump_domain, esm_dump_step, trim(cur_loc), cur_call
    else
      write(fname, '(A,"/esm_dump_",A,"_",I0,".json")') trim(dir), trim(scheme), cur_call
    end if
    open(newunit=unit, file=trim(fname), status='replace', action='write')
    write(unit, '(A,A,A,I0,A)') '{"scheme": "', trim(scheme), '", "call": ', cur_call, ', "vars": {'
    first_var = .true.
    if (len_trim(cur_loc) > 0) then
      ! The window variables are written uncut.
      cut = win_on
      win_on = .false.
      call dump_i_0('esm_domain', esm_dump_domain)
      call dump_i_0('esm_itimestep', esm_dump_step)
      if (cut) then
        call dump_i_0('box_i0', win_cut(1)); call dump_i_0('box_i1', win_cut(2))
        call dump_i_0('box_j0', win_cut(3)); call dump_i_0('box_j1', win_cut(4))
      end if
      win_on = cut
    end if
  end subroutine esm_dump_open

  ! Lower bounds (as 1-based offsets into an assumed-shape array) of the cut
  ! window, or .false. when the array is not memory- or tile-shaped.
  logical function cut_offsets(shp, o1, o2) result(cut)
    integer, intent(in) :: shp(2)
    integer, intent(out) :: o1, o2
    cut = .false.
    o1 = 0; o2 = 0
    if (.not. win_on) return
    if (shp(1) == win_mem(2) - win_mem(1) + 1 .and. shp(2) == win_mem(4) - win_mem(3) + 1) then
      o1 = win_cut(1) - win_mem(1) + 1; o2 = win_cut(3) - win_mem(3) + 1
      cut = .true.
    else if (shp(1) == win_tile(2) - win_tile(1) + 1 .and. shp(2) == win_tile(4) - win_tile(3) + 1) then
      o1 = win_cut(1) - win_tile(1) + 1; o2 = win_cut(3) - win_tile(3) + 1
      cut = .true.
    end if
  end function cut_offsets

  ! Open a dump at an explicit path (used by the standalone kernel drivers).
  subroutine esm_dump_open_file(scheme, path)
    character(len=*), intent(in) :: scheme, path
    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit, '(A,A,A,I0,A)') '{"scheme": "', trim(scheme), '", "call": ', 0, ', "vars": {'
    first_var = .true.
  end subroutine esm_dump_open_file

  subroutine esm_dump_close()
    write(unit, '(A)') '}}'
    close(unit)
    unit = -1
    esm_dump_inner = .false.
    win_on = .false.
    cur_loc = ''
  end subroutine esm_dump_close

  logical function esm_dump_active() result(active)
    active = unit /= -1
  end function esm_dump_active

  subroutine header(name, kind, shape)
    character(len=*), intent(in) :: name, kind
    integer, intent(in) :: shape(:)
    integer :: d
    if (.not. first_var) write(unit, '(A)') ','
    first_var = .false.
    write(unit, '(A,A,A,A,A)', advance='no') '"', trim(name), '": {"kind": "', kind, '", "shape": ['
    do d = 1, size(shape)
      if (d > 1) write(unit, '(A)', advance='no') ', '
      write(unit, '(I0)', advance='no') shape(d)
    end do
    write(unit, '(A)', advance='no') '], "data": ['
  end subroutine header

  subroutine footer()
    write(unit, '(A)', advance='no') ']}'
  end subroutine footer

  subroutine write_r8(vals)
    real(real64), intent(in) :: vals(:)
    integer :: n
    do n = 1, size(vals)
      if (n > 1) write(unit, '(A)', advance='no') ', '
      if (mod(n, 4) == 1 .and. n > 1) write(unit, '(A)') ''
      write(unit, '(ES24.16E3)', advance='no') vals(n)
    end do
  end subroutine write_r8

  subroutine dump_r8_0(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v
    call header(name, 'real64', [integer ::])
    call write_r8([v])
    call footer()
  end subroutine dump_r8_0
  subroutine dump_r8_1(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v(:)
    call header(name, 'real64', shape(v))
    call write_r8(v)
    call footer()
  end subroutine dump_r8_1
  subroutine dump_r8_2(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v(:,:)
    integer :: o1, o2, n1, n2
    if (cut_offsets(shape(v), o1, o2)) then
      n1 = win_cut(2) - win_cut(1); n2 = win_cut(4) - win_cut(3)
      call header(name, 'real64', [n1 + 1, n2 + 1])
      call write_r8(reshape(v(o1:o1+n1, o2:o2+n2), [(n1 + 1) * (n2 + 1)]))
    else
      call header(name, 'real64', shape(v))
      call write_r8(reshape(v, [size(v)]))
    end if
    call footer()
  end subroutine dump_r8_2
  subroutine dump_r8_3(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v(:,:,:)
    call header(name, 'real64', shape(v))
    call write_r8(reshape(v, [size(v)]))
    call footer()
  end subroutine dump_r8_3

  subroutine dump_r4_0(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v
    call header(name, 'real32', [integer ::])
    call write_r8([real(v, real64)])
    call footer()
  end subroutine dump_r4_0
  subroutine dump_r4_1(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v(:)
    call header(name, 'real32', shape(v))
    call write_r8(real(v, real64))
    call footer()
  end subroutine dump_r4_1
  subroutine dump_r4_2(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v(:,:)
    integer :: o1, o2, n1, n2
    if (cut_offsets(shape(v), o1, o2)) then
      n1 = win_cut(2) - win_cut(1); n2 = win_cut(4) - win_cut(3)
      call header(name, 'real32', [n1 + 1, n2 + 1])
      call write_r8(real(reshape(v(o1:o1+n1, o2:o2+n2), [(n1 + 1) * (n2 + 1)]), real64))
    else
      call header(name, 'real32', shape(v))
      call write_r8(real(reshape(v, [size(v)]), real64))
    end if
    call footer()
  end subroutine dump_r4_2
  subroutine dump_r4_3(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v(:,:,:)
    call header(name, 'real32', shape(v))
    call write_r8(real(reshape(v, [size(v)]), real64))
    call footer()
  end subroutine dump_r4_3

  subroutine write_i(vals)
    integer, intent(in) :: vals(:)
    integer :: n
    do n = 1, size(vals)
      if (n > 1) write(unit, '(A)', advance='no') ', '
      write(unit, '(I0)', advance='no') vals(n)
    end do
  end subroutine write_i
  subroutine dump_i_0(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v
    call header(name, 'int', [integer ::])
    call write_i([v])
    call footer()
  end subroutine dump_i_0
  subroutine dump_i_1(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v(:)
    call header(name, 'int', shape(v))
    call write_i(v)
    call footer()
  end subroutine dump_i_1
  subroutine dump_i_2(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v(:,:)
    integer :: o1, o2, n1, n2
    if (cut_offsets(shape(v), o1, o2)) then
      n1 = win_cut(2) - win_cut(1); n2 = win_cut(4) - win_cut(3)
      call header(name, 'int', [n1 + 1, n2 + 1])
      call write_i(reshape(v(o1:o1+n1, o2:o2+n2), [(n1 + 1) * (n2 + 1)]))
    else
      call header(name, 'int', shape(v))
      call write_i(reshape(v, [size(v)]))
    end if
    call footer()
  end subroutine dump_i_2

  subroutine write_l(vals)
    logical, intent(in) :: vals(:)
    integer :: n
    do n = 1, size(vals)
      if (n > 1) write(unit, '(A)', advance='no') ', '
      if (vals(n)) then
        write(unit, '(A)', advance='no') 'true'
      else
        write(unit, '(A)', advance='no') 'false'
      end if
    end do
  end subroutine write_l
  subroutine dump_l_0(name, v)
    character(len=*), intent(in) :: name
    logical, intent(in) :: v
    call header(name, 'bool', [integer ::])
    call write_l([v])
    call footer()
  end subroutine dump_l_0
  subroutine dump_l_1(name, v)
    character(len=*), intent(in) :: name
    logical, intent(in) :: v(:)
    call header(name, 'bool', shape(v))
    call write_l(v)
    call footer()
  end subroutine dump_l_1

  subroutine dump_str(name, v)
    character(len=*), intent(in) :: name, v
    if (.not. first_var) write(unit, '(A)') ','
    first_var = .false.
    write(unit, '(A,A,A,A,A)', advance='no') '"', trim(name), '": {"kind": "str", "shape": [], "data": ["', trim(v), '"]}'
  end subroutine dump_str

end module module_esm_dump
