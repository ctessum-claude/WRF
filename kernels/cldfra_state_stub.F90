! Stub of the four module_state_description parameters that cal_cldfra1 imports.
! frame/module_state_description.F is generated from the Registry and pulls in
! the whole package/namelist description; cal_cldfra1 reads only these, so the
! standalone kernel build declares them with the generated file's values
! (frame/module_state_description.F:14,15,161,168 of WRF 4.8).
module module_state_description
  implicit none
  integer, parameter :: fer_mp_hires = 5
  integer, parameter :: fer_mp_hires_advect = 15
  integer, parameter :: kfetascheme = 1
  integer, parameter :: kfcupscheme = 10
end module module_state_description
