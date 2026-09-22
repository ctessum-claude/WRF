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
 real(kind=kind_phys), allocatable, public:: d_dtfac(:,:), d_taup(:,:), d_velco(:,:), d_usqj(:,:)""")

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

   lcap   = kte
   lcapp1 = lcap + 1""")

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
   d_fr(its:ite) = fr(its:ite)
   d_dtfac = dtfac ; d_taup = taup ; d_velco = velco ; d_usqj = usqj
!
! rotate tendencies from zonal/meridional back to model grid""")

open(out, "w").write(s)
print(f"wrote {out}")
