! Link stubs for kernels/mosaic_therm_driver: three symbols libwrflib.a's start_domain
! and io_int reference but that the thermodynamics never calls (they live in the
! em core and in frame/ outside the library).  Calling any of them is an error.
subroutine start_domain_em()
  stop 'start_domain_em stub called'
end subroutine
subroutine wrf_sizeof_integer(retval)
  integer retval
  retval = 4
end subroutine
subroutine wrf_sizeof_real(retval)
  integer retval
  retval = 4
end subroutine
