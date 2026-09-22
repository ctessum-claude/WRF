! Minimal stand-in for frame/module_wrf_error, so WRF-native physics modules
! that only report fatal errors (e.g. module_sf_noahlsm) build standalone.
module module_wrf_error
  implicit none
contains
  subroutine wrf_error_fatal(msg)
    character(len=*), intent(in) :: msg
    write(0, '(A)') 'FATAL: '//trim(msg)
    error stop 1
  end subroutine wrf_error_fatal
  subroutine wrf_message(msg)
    character(len=*), intent(in) :: msg
    write(0, '(A)') trim(msg)
  end subroutine wrf_message
  subroutine wrf_debug(level, msg)
    integer, intent(in) :: level
    character(len=*), intent(in) :: msg
  end subroutine wrf_debug
end module module_wrf_error
