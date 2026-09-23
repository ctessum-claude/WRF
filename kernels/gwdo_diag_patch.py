#!/usr/bin/env python3
"""Generate build/bl_gwdo_diag.F90: a copy of physics_mmm/bl_gwdo.F90 that also
exposes the branch decisions of each column, so the same row can be run at
real32 and real64 and the two branch records compared.

Nothing in the physics is altered: only module-level diagnostic arrays are
added, filled where the branch is taken, and the module/subroutine are renamed
so the diagnostic copy can be linked beside the original.

  gwdo_diag_patch.py <bl_gwdo.F90> <out.F90>
"""
import sys

src, out = sys.argv[1], sys.argv[2]
s = open(src).read()

def sub(old, new, n=1):
    global s
    assert s.count(old) == n, (s.count(old), old[:70])
    s = s.replace(old, new)

sub("\n module bl_gwdo\n use ccpp_kind_types", "\n module bl_gwdo_diag\n use ccpp_kind_types")
sub("\n end module bl_gwdo\n", "\n end module bl_gwdo_diag\n")
sub("   subroutine bl_gwdo_run(", "   subroutine bl_gwdo_diag_run(")
sub("   end subroutine bl_gwdo_run\n", "   end subroutine bl_gwdo_diag_run\n")
sub(""" public:: bl_gwdo_run,     &
          bl_gwdo_init,    &
          bl_gwdo_finalize""", """ public:: bl_gwdo_diag_run,     &
          bl_gwdo_init,    &
          bl_gwdo_finalize
! EarthSciML diagnostics: the branch record of the last call.
 integer, allocatable, public:: d_kref(:), d_komax(:), d_kbomax(:), d_kblk(:)
 logical, allocatable, public:: d_ldrag(:), d_icrilv(:)
 integer, allocatable, public:: d_sat(:,:)   ! per level: 0 no test, 1 saturation, 2 no wavebreaking
 real(kind=kind_phys), allocatable, public:: d_bnv2low(:), d_ulow(:), d_fr(:), d_taub(:), d_zref(:)
 real(kind=kind_phys), allocatable, public:: d_dtfac(:,:), d_taup(:,:), d_velco(:,:), d_usqj(:,:)
! EarthSciML stage 2: the intermediate fields each .esm sub-component produces.
 integer, allocatable, public:: d_nwd(:), d_icrilvk(:,:)
 real(kind=kind_phys), allocatable, public:: d_ubar(:), d_vbar(:), d_rhobar(:), d_delks(:), d_delks1(:), &
          d_usqjlow(:), d_oa(:), d_ol(:), d_olp(:), d_oc(:), d_od(:), d_dxy(:), d_dxyp(:), d_dxeff(:), &
          d_area(:), d_coefm(:), d_efac(:), d_gfac(:), d_xlinv(:), d_nhd(:), d_xn(:), d_yn(:), &
          d_zblk(:), d_fbdcd(:), d_dusfc(:), d_dvsfc(:)
 real(kind=kind_phys), allocatable, public:: d_rho(:,:), d_del(:,:), d_u1(:,:), d_v1(:,:), d_vtj(:,:), &
          d_vtk(:,:), d_usqjr(:,:), d_bnv2r(:,:), d_taupw(:,:), d_taufb(:,:), d_taudr(:,:), &
          d_dtaux2d(:,:), d_dtauy2d(:,:), d_satfac(:,:), d_rim(:,:)""")

# allocate and clear at the top of the routine, after the constants block
sub("""   lcap   = kte
   lcapp1 = lcap + 1""", """   if (allocated(d_kref)) deallocate(d_kref, d_komax, d_kbomax, d_kblk, d_ldrag, d_icrilv, d_sat, &
                                     d_bnv2low, d_ulow, d_fr, d_taub, d_zref, d_dtfac, d_taup, d_velco, d_usqj)
   allocate(d_kref(its:ite), d_komax(its:ite), d_kbomax(its:ite), d_kblk(its:ite), &
            d_ldrag(its:ite), d_icrilv(its:ite), d_sat(its:ite,kts:kte), &
            d_bnv2low(its:ite), d_ulow(its:ite), d_fr(its:ite), d_taub(its:ite), d_zref(its:ite), &
            d_dtfac(its:ite,kts:kte), d_taup(its:ite,kts:kte+1), d_velco(its:ite,kts:kte-1), &
            d_usqj(its:ite,kts:kte))
   d_kblk = 0 ; d_sat = 0 ; d_fr = 0. ; d_taub = 0.
   if (allocated(d_nwd)) deallocate(d_nwd, d_icrilvk, d_ubar, d_vbar, d_rhobar, d_delks, d_delks1, &
            d_usqjlow, d_oa, d_ol, d_olp, d_oc, d_od, d_dxy, d_dxyp, d_dxeff, d_area, d_coefm, &
            d_efac, d_gfac, d_xlinv, d_nhd, d_xn, d_yn, d_zblk, d_fbdcd, d_dusfc, d_dvsfc, &
            d_rho, d_del, d_u1, d_v1, d_vtj, d_vtk, d_usqjr, d_bnv2r, d_taupw, d_taufb, d_taudr, &
            d_dtaux2d, d_dtauy2d, d_satfac, d_rim)
   allocate(d_nwd(its:ite), d_icrilvk(its:ite,kts:kte), d_ubar(its:ite), d_vbar(its:ite), &
            d_rhobar(its:ite), d_delks(its:ite), d_delks1(its:ite), d_usqjlow(its:ite), &
            d_oa(its:ite), d_ol(its:ite), d_olp(its:ite), d_oc(its:ite), d_od(its:ite), &
            d_dxy(its:ite), d_dxyp(its:ite), d_dxeff(its:ite), d_area(its:ite), d_coefm(its:ite), &
            d_efac(its:ite), d_gfac(its:ite), d_xlinv(its:ite), d_nhd(its:ite), d_xn(its:ite), &
            d_yn(its:ite), d_zblk(its:ite), d_fbdcd(its:ite), d_dusfc(its:ite), d_dvsfc(its:ite), &
            d_rho(its:ite,kts:kte), d_del(its:ite,kts:kte), d_u1(its:ite,kts:kte), &
            d_v1(its:ite,kts:kte), d_vtj(its:ite,kts:kte), d_vtk(its:ite,kts:kte), &
            d_usqjr(its:ite,kts:kte), d_bnv2r(its:ite,kts:kte), d_taupw(its:ite,kts:kte+1), &
            d_taufb(its:ite,kts:kte+1), d_taudr(its:ite,kts:kte), d_dtaux2d(its:ite,kts:kte), &
            d_dtauy2d(its:ite,kts:kte), d_satfac(its:ite,kts:kte), d_rim(its:ite,kts:kte))
   d_icrilvk = 0 ; d_coefm = 0. ; d_efac = 0. ; d_gfac = 0. ; d_xlinv = 0. ; d_nhd = 0.
   d_zblk = 0. ; d_fbdcd = 0. ; d_satfac = 0. ; d_rim = 0.

   lcap   = kte
   lcapp1 = lcap + 1""")

# the RAW Richardson number and N^2 profiles, before the low-level weighted
# average overwrites level 1 and floods it down through kref-1
sub("""!
! compute the low-level wind speed
!""", """   d_usqjr = usqj ; d_bnv2r = bnv2
!
! compute the low-level wind speed
!""")

# the wind-rose sector, a scalar inside the orography loop
sub("""     nwd   = nwdir(idir)
""", """     nwd   = nwdir(idir)
     d_nwd(i) = nwd
""")

# the per-level critical-level flag (icrilv is a running OR over k >= kref)
sub("""         icrilv(i) = icrilv(i) .or. ( usqj(i,k) < ric)                         &
                               .or. (velco(i,k) <= 0.)
""", """         icrilv(i) = icrilv(i) .or. ( usqj(i,k) < ric)                         &
                               .or. (velco(i,k) <= 0.)
         if (icrilv(i)) d_icrilvk(i,k) = 1
""")

# the wave-launch scalars, inside the .not.ldrag branch where they are defined
# (fr and coefm are never initialized, so they must be read there, not at the end)
sub("""       bnv(i) = sqrt( bnv2(i,1) )
       fr(i) = bnv(i)  * rulow(i) * var(i) * od(i)
       fr(i) = min(fr(i),frmax)""", """       bnv(i) = sqrt( bnv2(i,1) )
       fr(i) = bnv(i)  * rulow(i) * var(i) * od(i)
       fr(i) = min(fr(i),frmax)
       d_fr(i) = fr(i)""")
sub("""       taub(i)  = rhobar(i) * efac * xlinv * gfac * ulow(i)*ulow(i)*ulow(i)    &
                  / bnv(i)
""", """       taub(i)  = rhobar(i) * efac * xlinv * gfac * ulow(i)*ulow(i)*ulow(i)    &
                  / bnv(i)
       d_efac(i) = efac ; d_gfac(i) = gfac ; d_xlinv(i) = xlinv ; d_coefm(i) = coefm(i)
       d_nhd(i) = 0.
""")
sub("""         taub(i)  = taub(i) * (1.+nhd_effect)
""", """         taub(i)  = taub(i) * (1.+nhd_effect)
         d_nhd(i) = nhd_effect
""")

# the minimum Richardson number of Shutts (1985) and the saturation factor
# 2*sqrt(temc)-temc, which vanishes at Ri = ric = 0.25 (FORTRAN_BUGS N85)
sub("""           rim  = usqj(i,k) * (1.-fro) / (tem * tem)
""", """           rim  = usqj(i,k) * (1.-fro) / (tem * tem)
           d_rim(i,k) = rim
""")
sub("""             temc = 2.0 + 1.0 / tem2
""", """             temc = 2.0 + 1.0 / tem2
             d_satfac(i,k) = 2.*sqrt(temc)-temc
""")

# the wave stress profile BEFORE the flow-blocking stress is added to it
sub("""   do i = its,ite
     if (.not.ldrag(i)) then
!
! determine the height of flow-blocking layer""", """   d_taupw = taup
   do i = its,ite
     if (.not.ldrag(i)) then
!
! determine the height of flow-blocking layer""")

# the flow-blocking drag coefficient and layer depth
sub("""         fbdcd = max(2.0-1.0/od(i),0.)
""", """         fbdcd = max(2.0-1.0/od(i),0.)
         d_fbdcd(i) = fbdcd ; d_zblk(i) = zblk
""")

# the stress divergence BEFORE the wind-reversal / mesosphere limiter
sub("""! apply damping factor to prevent wind reversal""",
    """   d_taudr = taud
! apply damping factor to prevent wind reversal""")

# the flow-blocking layer index, which is a scalar inside the i loop
sub("""       if (kblk /= 0) then
!
! compute flow-blocking stress""", """       d_kblk(i) = kblk
       if (kblk /= 0) then
!
! compute flow-blocking stress""")

# the saturation-hypothesis branch
sub("""           if (rim <= ric) then  ! saturation hypothesis!""",
    """           d_sat(i,k) = 2
           if (rim <= ric) then  ! saturation hypothesis!
             d_sat(i,k) = 1""")

# record the rest just before the tendencies are rotated back to the grid
sub("""! rotate tendencies from zonal/meridional back to model grid""",
    """   d_kref = kref ; d_komax = komax ; d_kbomax = kbomax
   d_ldrag = ldrag ; d_icrilv = icrilv
   d_ulow = ulow ; d_taub = taub ; d_zref = zref
   d_bnv2low(its:ite) = bnv2(its:ite,1)
   d_dtfac = dtfac ; d_taup = taup ; d_velco = velco ; d_usqj = usqj
   d_ubar = ubar ; d_vbar = vbar ; d_rhobar = rhobar ; d_delks = delks ; d_delks1 = delks1
   d_usqjlow(its:ite) = usqj(its:ite,1)
   d_oa = oa ; d_ol = ol ; d_olp = olp ; d_oc = oc ; d_od = od
   d_dxy = dxy ; d_dxyp = dxyp ; d_dxeff = dx_eff ; d_area = area
   d_xn = xn ; d_yn = yn ; d_dusfc = dusfc ; d_dvsfc = dvsfc
   d_rho = rho ; d_del = del ; d_u1 = u1 ; d_v1 = v1 ; d_vtj = vtj ; d_vtk = vtk
   d_taufb = taufb ; d_dtaux2d = dtaux2d ; d_dtauy2d = dtauy2d
!
! rotate tendencies from zonal/meridional back to model grid""")

open(out, "w").write(s)
print(f"wrote {out}")
