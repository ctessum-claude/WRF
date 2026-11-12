! Minimal stand-in for WRF's generated module_configure, carrying only the
! grid_config_rec_type fields the WRF-Chem emission kernels read
! (module_cbmz_addemiss, module_mosaic_addemiss, module_bioemi_simple,
! add_biogenics), so they build standalone in the kernel drivers.
module module_configure
  implicit none
  type grid_config_rec_type
    integer :: chem_opt = 0
    integer :: emiss_opt = 0
    integer :: emiss_inpt_opt = 0
    integer :: bio_emiss_opt = 0
    integer :: biomass_burn_opt = 0
    integer :: kemit = 1
    integer :: num_soil_layers = 4
    integer :: start_month = 1
  end type grid_config_rec_type
end module module_configure
