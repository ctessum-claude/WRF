! esm_flat_input: reader for the flat text produced by EqWeFiC tools/esm_dump.py --flat.
! Records: "<name> <kind> <rank> <extents...>" then one line of values (Fortran order).
module esm_flat_input
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: flat_read, flat_r0, flat_r1, flat_r2, flat_r3, flat_i0, flat_i1, flat_l0, flat_has

  type :: rec_t
    character(len=64) :: name = ''
    character(len=8)  :: kind = ''
    integer :: rank = 0
    integer :: dims(3) = 1
    real(real64), allocatable :: r(:)
    integer, allocatable :: i(:)
    logical, allocatable :: l(:)
  end type
  type(rec_t), allocatable :: recs(:)
  integer :: nrec = 0

contains

  subroutine flat_read(path)
    character(len=*), intent(in) :: path
    integer :: u, ios, n, k, rank
    character(len=64) :: name
    character(len=8) :: kind
    character(len=32) :: tok(3)
    allocate(recs(512))
    open(newunit=u, file=path, status='old', action='read')
    do
      read(u, *, iostat=ios) name, kind, rank
      if (ios /= 0) exit
      backspace(u)
      nrec = nrec + 1
      recs(nrec)%name = name; recs(nrec)%kind = kind; recs(nrec)%rank = rank
      recs(nrec)%dims = 1
      select case (rank)
      case (0); read(u, *) name, kind, rank
      case (1); read(u, *) name, kind, rank, recs(nrec)%dims(1)
      case (2); read(u, *) name, kind, rank, recs(nrec)%dims(1:2)
      case (3); read(u, *) name, kind, rank, recs(nrec)%dims(1:3)
      end select
      n = product(recs(nrec)%dims)
      select case (trim(kind))
      case ('real'); allocate(recs(nrec)%r(n)); read(u, *) recs(nrec)%r
      case ('int');  allocate(recs(nrec)%i(n)); read(u, *) recs(nrec)%i
      case ('bool'); allocate(recs(nrec)%l(n)); read(u, *) recs(nrec)%l
      case ('str');  read(u, '(A)') name   ! strings are ignored by the drivers
      end select
    end do
    close(u)
  end subroutine flat_read

  integer function find(name) result(k)
    character(len=*), intent(in) :: name
    do k = 1, nrec
      if (trim(recs(k)%name) == trim(name)) return
    end do
    write(*, '(A,A)') 'esm_flat_input: missing variable ', trim(name)
    stop 2
  end function find

  logical function flat_has(name)
    character(len=*), intent(in) :: name
    integer :: k
    flat_has = .false.
    do k = 1, nrec
      if (trim(recs(k)%name) == trim(name)) flat_has = .true.
    end do
  end function flat_has

  real(real64) function flat_r0(name)
    character(len=*), intent(in) :: name
    flat_r0 = recs(find(name))%r(1)
  end function flat_r0
  subroutine flat_r1(name, a)
    character(len=*), intent(in) :: name
    real(real64), intent(out) :: a(:)
    a = recs(find(name))%r(1:size(a))
  end subroutine flat_r1
  subroutine flat_r2(name, a)
    character(len=*), intent(in) :: name
    real(real64), intent(out) :: a(:,:)
    a = reshape(recs(find(name))%r(1:size(a)), shape(a))
  end subroutine flat_r2
  subroutine flat_r3(name, a)
    character(len=*), intent(in) :: name
    real(real64), intent(out) :: a(:,:,:)
    a = reshape(recs(find(name))%r(1:size(a)), shape(a))
  end subroutine flat_r3
  integer function flat_i0(name)
    character(len=*), intent(in) :: name
    flat_i0 = recs(find(name))%i(1)
  end function flat_i0
  subroutine flat_i1(name, a)
    character(len=*), intent(in) :: name
    integer, intent(out) :: a(:)
    a = recs(find(name))%i(1:size(a))
  end subroutine flat_i1
  logical function flat_l0(name)
    character(len=*), intent(in) :: name
    flat_l0 = recs(find(name))%l(1)
  end function flat_l0
end module esm_flat_input
