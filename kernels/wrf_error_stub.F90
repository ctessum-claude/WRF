! wrf_error_stub: empty stand-in for WRF's module_wrf_error, which radconst
! USEs but never calls; lets solar_driver compile radconst outside WRF.
module module_wrf_error
end module module_wrf_error
