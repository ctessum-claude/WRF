! emissions_driver: replay one dumped column (esm dump scheme 'emissions',
! written by the hook around the emissions_driver call in chem/chem_driver.F)
! through the routines chem/emissions_driver.F calls for chem_opt = 170 with
! emiss_opt = 4, bio_emiss_opt = 1 and no dust, sea salt, DMS, volcanic or
! fire emissions, in emissions_driver's own order:
!   bio_emissions        (module_bioemi_simple, GUNTHER1: e_bio)
!   cbmz_addemiss_anthro (module_cbmz_addemiss: anthropogenic gases)
!   add_biogenics        (module_emissions_anthropogenics: e_bio into chem)
!   mosaic_addemiss      (module_mosaic_addemiss: anthropogenic aerosol)
! and write the emitted column, e_bio and each stage's increment as an esm dump.
!   emissions_driver <input.flat> <output.json>
! Built real64 (-fdefault-real-8).  The flat must carry the run's species and
! emission pointers (p_*), which WRF sets at run time from the registry
! packages; every pointer the flat omits keeps its default of 1, as it does in
! the model for a species the package lacks.  The column is i = j = 1.
program emissions_driver_replay
  use esm_flat_input
  use module_esm_dump
  use module_configure, only: grid_config_rec_type
  use module_state_description
  use module_data_mosaic_asect
  use module_bioemi_simple, only: bio_emissions
  use module_cbmz_addemiss, only: cbmz_addemiss_anthro
  use add_biogenics_kernel, only: add_biogenics
  use module_mosaic_addemiss, only: mosaic_addemiss
  implicit none
  character(len=1024) :: inpath, outpath
  type(grid_config_rec_type) :: cf
  integer :: kts, kte, kme, kemit, ne_area, ktau, julday, nsize, n
  real :: dt, gmt, dx, g
  integer :: ivgtyp(1,1), isltyp(1,1)
  real :: gsw(1,1), vegfra(1,1), z2(1,1)
  real, allocatable :: chem(:,:,:,:), moist(:,:,:,:), e_bio(:,:,:), emis_ant(:,:,:,:)
  real, allocatable :: t_phy(:,:,:), p_phy(:,:,:), rho(:,:,:), alt(:,:,:), dz8w(:,:,:), z_at_w(:,:,:), zero3(:,:,:)
  real, allocatable :: chem_in(:,:), col(:), r2(:,:), r1(:), c0(:,:), c1(:,:), c2(:,:)
  real, allocatable :: ebu(:,:,:,:), smois(:,:,:), erod(:,:,:), emis_dust(:,:,:,:)
  integer, allocatable :: i1(:)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  call flat_read(trim(inpath))
  kts = flat_i0('kts'); kte = flat_i0('kte'); kemit = flat_i0('kemit')
  num_chem = flat_i0('num_chem'); num_emis_ant = flat_i0('num_emis_ant'); ne_area = flat_i0('ne_area')
  cf%chem_opt = flat_i0('chem_opt'); cf%emiss_opt = flat_i0('emiss_opt')
  cf%emiss_inpt_opt = flat_i0('emiss_inpt_opt'); cf%bio_emiss_opt = flat_i0('bio_emiss_opt')
  cf%biomass_burn_opt = 0; cf%kemit = kemit
  ktau = flat_i0('ktau'); julday = flat_i0('julday')
  dt = flat_r0('dt'); gmt = flat_r0('gmt'); dx = 20000.; g = 9.81
  ivgtyp(1,1) = flat_i0('ivgtyp'); isltyp = 0
  gsw(1,1) = flat_r0('gsw'); vegfra(1,1) = flat_r0('vegfra'); z2 = 0.
#include "emissions_pointers.inc"
  num_moist = 4; p_qv = 2; p_qc = 3; p_qr = 4
  num_ebu = 1

  ! MOSAIC section data that init_data_mosaic_asect sets at run time
  nsize = flat_i0('nsize_aer')
  ntype_aer = 1; nphase_aer = 1; ai_phase = 1; nsize_aer(1) = nsize
  allocate(r1(nsize), i1(nsize))
  call flat_r1('volumcen_sect', r1); volumcen_sect(1:nsize,1) = r1
  lptr_so4_aer = -999888777; lptr_no3_aer = -999888777; lptr_cl_aer = -999888777
  lptr_msa_aer = -999888777; lptr_co3_aer = -999888777; lptr_nh4_aer = -999888777
  lptr_na_aer = -999888777; lptr_ca_aer = -999888777; lptr_oin_aer = -999888777
  lptr_oc_aer = -999888777; lptr_bc_aer = -999888777; numptr_aer = -999888777
  call flat_i1('lptr_so4_aer', i1); lptr_so4_aer(1:nsize,1,1) = i1
  call flat_i1('lptr_no3_aer', i1); lptr_no3_aer(1:nsize,1,1) = i1
  call flat_i1('lptr_oc_aer', i1);  lptr_oc_aer(1:nsize,1,1) = i1
  call flat_i1('lptr_bc_aer', i1);  lptr_bc_aer(1:nsize,1,1) = i1
  call flat_i1('lptr_oin_aer', i1); lptr_oin_aer(1:nsize,1,1) = i1
  call flat_i1('numptr_aer', i1);   numptr_aer(1:nsize,1,1) = i1

  kme = kte + 1   ! WRF's memory extent: kte = kme - 1
  allocate(chem(1,kts:kme,1,num_chem), moist(1,kts:kme,1,num_moist), e_bio(1,1,ne_area))
  allocate(emis_ant(1,1:kemit,1,num_emis_ant), r2(kemit,num_emis_ant))
  allocate(t_phy(1,kts:kme,1), p_phy(1,kts:kme,1), rho(1,kts:kme,1), alt(1,kts:kme,1))
  allocate(dz8w(1,kts:kme,1), z_at_w(1,kts:kme,1), zero3(1,kts:kme,1))
  allocate(chem_in(kts:kte,num_chem), col(kts:kte), c0(kts:kte,num_chem), c1(kts:kte,num_chem), c2(kts:kte,num_chem))
  allocate(ebu(1,kts:kme,1,1), smois(1,cf%num_soil_layers,1), erod(1,1,3), emis_dust(1,1,1,1))
  ebu = 0.; smois = 0.; erod = 0.; emis_dust = 0.; zero3 = 0.
  chem = 0.; t_phy = 300.; p_phy = 1.e5; rho = 1.; alt = 1.; dz8w = 1.
  moist = 0.
  call flat_r1('qv', col); moist(1,kts:kte,1,p_qv) = col
  call flat_r1('qc', col); moist(1,kts:kte,1,p_qc) = col
  call flat_r1('qr', col); moist(1,kts:kte,1,p_qr) = col
  call flat_r1('t_phy', col); t_phy(1,kts:kte,1) = col
  call flat_r1('p_phy', col); p_phy(1,kts:kte,1) = col
  call flat_r1('rho', col);   rho(1,kts:kte,1) = col
  call flat_r1('alt', col);   alt(1,kts:kte,1) = col
  call flat_r1('dz8w', col);  dz8w(1,kts:kte,1) = col
  block
    real :: zw(kts:kte+1)
    call flat_r1('z_at_w', zw); z_at_w(1,:,1) = zw
  end block
  call flat_r2('emis_ant', r2); emis_ant(1,:,1,:) = r2
  emis_ant(1,:,1,1) = 0.   ! the dummy slot every absent p_e_* points at
  call flat_r2('chem_in', chem_in); chem(1,kts:kte,1,:) = chem_in
  e_bio = 0.

  call bio_emissions(1, ktau, dt, dx, cf, gmt, julday, t_phy, moist, zero3, zero3, &
       e_bio, p_phy, chem, rho, dz8w, ne_area, ivgtyp, gsw, vegfra, z2, z2, z2, z2, z2, z_at_w, &
       1,1, 1,1, kts,kme,  1,1, 1,1, kts,kme,  1,1, 1,1, kts,kte)
  c0 = chem(1,kts:kte,1,:)
  call cbmz_addemiss_anthro(1, dt, dz8w, cf, rho, chem, emis_ant, alt, &
       1,1, 1,1, kts,kme,  1,1, 1,1, kts,kme,  1,1, 1,1, kts,kte)
  c1 = chem(1,kts:kte,1,:)
  call add_biogenics(1, dt, dz8w, cf, rho, chem, e_bio, ne_area, &
       z2,z2,z2,z2,z2, z2,z2,z2,z2,z2, z2,z2,z2,z2,z2,z2, z2,z2, &
       1,1, 1,1, kts,kme,  1,1, 1,1, kts,kme,  1,1, 1,1, kts,kte)
  c2 = chem(1,kts:kte,1,:)
  call mosaic_addemiss(1, dt, z2, z2, alt, dz8w, z2, cf, chem, z2, z2, smois, ivgtyp, isltyp, &
       emis_ant, ebu, 0, 0, ktau, zero3, zero3, rho, g, dx, erod, 0, 0, z2, emis_dust, &
       1,1, 1,1, kts,kme,  1,1, 1,1, kts,kme,  1,1, 1,1, kts,kte)

  call esm_dump_open_file('emissions', trim(outpath))
  call esm_dump_var('real_bytes', storage_size(dt)/8)
  call esm_dump_var('kts', kts); call esm_dump_var('kte', kte); call esm_dump_var('num_chem', num_chem)
  call esm_dump_var('e_bio', e_bio(1,1,:))
  call esm_dump_var('d_anthro_gas', c1 - c0)
  call esm_dump_var('d_bio', c2 - c1)
  call esm_dump_var('d_aer', chem(1,kts:kte,1,:) - c2)
  call esm_dump_var('chem_out', chem(1,kts:kte,1,:))
  call esm_dump_close()
end program emissions_driver_replay
