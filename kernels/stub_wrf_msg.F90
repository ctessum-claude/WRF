! External wrf_message / wrf_debug / wrf_error_fatal entry points, for WRF
! modules (e.g. chem/module_peg_util.F) that call them as plain externals.
subroutine wrf_message(msg)
  character(len=*), intent(in) :: msg
  write(0,'(A)') trim(msg)
end subroutine wrf_message
subroutine wrf_debug(level, msg)
  integer, intent(in) :: level
  character(len=*), intent(in) :: msg
end subroutine wrf_debug
subroutine wrf_error_fatal(msg)
  character(len=*), intent(in) :: msg
  write(0,'(A)') 'FATAL: '//trim(msg)
  error stop 1
end subroutine wrf_error_fatal
