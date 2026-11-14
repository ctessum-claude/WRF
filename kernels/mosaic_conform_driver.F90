! mosaic_conform_driver: run MOSAIC's conform_electrolytes(jtotal) (chem/module_mosaic_therm.F)
! on given bin compositions, through WRF's OWN compiled objects (main/libwrflib.a, the
! chem170 build; the routine declares real(kind=8)), and print its inputs and outputs
! at full precision.  EqWeFiC components/aerosol/mosaic/conform_electrolytes.esm.
!
!   mosaic_conform_driver <cases.txt> <out.txt>
!
! cases.txt: ncase, then per case one line of 14 reals:
!   so4 no3 cl nh4 msa na ca co3 oin oc bc  (aer(:, jtotal, 1), nmol/m3 (oin, oc, bc ng/m3))
!   hno3_g hcl_g nh3_g
! A NEGATIVE hno3_g carries all three gases over from the previous case's output, so a
! sequence of cases replays mosaic_dynamic_solver's bin loop, where the gases one bin's
! conform_electrolytes degasses are what the next bin sees.  Each case first runs
! check_aerosol_mass, as that loop does; a bin it finds below mass_cutoff is NOT conformed
! (called = 0, outputs are the inputs), as in the model.
! out.txt, per case one line:
!   case called xt  aer(so4 no3 cl nh4 msa na ca co3, jtotal)  electrolyte(1:22, jtotal)
!   gas(hno3 hcl nh3)  total_species(ihcl_g)  tot_cl_in
! (total_species(ihcl_g) and tot_cl_in start at 0, so they print form_nacl's
! chlorine creation when the gas goes negative.)
program mosaic_conform_driver
  use module_data_mosaic_therm
  use module_mosaic_therm
  implicit none
  character(len=1024) :: inpath, outpath
  integer :: ncase, ic, j, idx(11), called
  real(kind=8) :: v(14), xt, g3(3)

  call get_command_argument(1, inpath); call get_command_argument(2, outpath)
  nbin_a = 1
  call load_mosaic_parameters
  idx = (/ iso4_a, ino3_a, icl_a, inh4_a, imsa_a, ina_a, ica_a, ico3_a, ioin_a, ioc_a, ibc_a /)
  open(10, file=trim(inpath), status='old')
  open(11, file=trim(outpath), status='replace')
  read(10,*) ncase
  do ic = 1, ncase
    read(10,*) v
    g3 = (/ gas(ihno3_g), gas(ihcl_g), gas(inh3_g) /)
    call initialize_mosaic_variables
    do j = 1, 11
      aer(idx(j), jtotal, 1) = v(j)
    enddo
    if (v(12) >= 0.0d0) g3 = v(12:14)
    gas(ihno3_g) = g3(1); gas(ihcl_g) = g3(2); gas(inh3_g) = g3(3)
    total_species(ihcl_g) = 0.0d0
    tot_cl_in = 0.0d0
    electrolyte(:, :, 1) = 0.0d0
    call check_aerosol_mass(1)   ! initialize_mosaic_variables set jaerosolstate = -1
    called = 0; xt = 0.0d0
    if (jaerosolstate(1) /= no_aerosol) then
      call conform_electrolytes(jtotal, 1, xt)
      called = 1
    endif
    write(11,'(2i6,40es26.17e3)') ic, called, xt, (aer(idx(j), jtotal, 1), j=1,8), &
      (electrolyte(j, jtotal, 1), j=1,nelectrolyte), gas(ihno3_g), gas(ihcl_g), gas(inh3_g), &
      total_species(ihcl_g), tot_cl_in
  enddo
  close(10); close(11)
end program mosaic_conform_driver
