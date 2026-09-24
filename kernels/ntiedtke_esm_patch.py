#!/usr/bin/env python3
"""Generate build/cu_ntiedtke_esm.F90: a copy of physics_mmm/cu_ntiedtke.F90 with
the EarthSciML instrumentation spliced in.  The submodule source stays untouched.

Two things are added, and NOTHING in the physics changes:

1. The existing inner dump after the `cumastrn` call (ntiedtke_esm_inner.inc),
   which publishes the convection type, cloud base and top and cumastrn's own
   tendencies and mass fluxes.  It used to be spliced by a `sed` expression in
   the Makefile; it is done here instead so that one script owns every splice.

2. A TAPE of every `cuadjtqn` call.  cuadjtqn is the scheme's moist-adjustment
   kernel and is called from four places (the environmental saturation in
   cuinin, the test parcel in cutypen, the updraught in cuascn and cubasmcn,
   and the downdraught in cudlfsn/cuddrafn), once per level per call site, so a
   single dumped step yields several hundred independent (p, T, q) adjustments
   across all three of its kcall modes.  Nothing else in the scheme publishes
   them: they are in-place updates of arrays the caller owns, and by the time
   cumastrn returns they have been overwritten many times over.  The tape
   records, per call and per column, the input state, the state after the FIRST
   Newton step, and the final state -- the two-step structure the routine has,
   split so that a test can pin each half on its own.

  ntiedtke_esm_patch.py <cu_ntiedtke.F90> <out.F90>
"""
import sys

src, out = sys.argv[1], sys.argv[2]
s = open(src).read()

def sub(old, new, n=1):
    global s
    assert s.count(old) == n, (s.count(old), old[:80])
    s = s.replace(old, new)

# ---------------------------------------------------------------- module use
sub("""
 module cu_ntiedtke
 use ccpp_kind_types,only: kind_phys
 use cu_ntiedtke_common
""", """
 module cu_ntiedtke
 use ccpp_kind_types,only: kind_phys
 use cu_ntiedtke_common
 use module_esm_dump
""")

# ------------------------------------------------------- the cuadjtqn tape
sub(""" public:: cu_ntiedtke_run,     &
          cu_ntiedtke_init,    &
          cu_ntiedtke_finalize
""", """ public:: cu_ntiedtke_run,     &
          cu_ntiedtke_init,    &
          cu_ntiedtke_finalize
! EarthSciML: a tape of every cuadjtqn call of the last cu_ntiedtke_run.
! Rows are (call index, column); adj_n is how many calls were recorded.
 integer,parameter,public:: adj_max = 4000
 integer,public:: adj_n = 0
 integer,allocatable,public:: adj_kcall(:), adj_kk(:)
 integer,allocatable,public:: adj_flag(:,:)
 real(kind=kind_phys),allocatable,public:: adj_psp(:,:), adj_t0(:,:), adj_q0(:,:), &
                                           adj_tm(:,:), adj_qm(:,:), adj_t1(:,:), adj_q1(:,:)
""")

# allocate and reset at the top of cu_ntiedtke_run, before anything is computed
sub("""!-----------------------------------------------------------------------------------------------------------------
!
      ztmst=dt
""", """!-----------------------------------------------------------------------------------------------------------------
!
      if (allocated(adj_kcall)) deallocate(adj_kcall, adj_kk, adj_flag, adj_psp, &
                                           adj_t0, adj_q0, adj_tm, adj_qm, adj_t1, adj_q1)
      allocate(adj_kcall(adj_max), adj_kk(adj_max), adj_flag(adj_max,lq), &
               adj_psp(adj_max,lq), adj_t0(adj_max,lq), adj_q0(adj_max,lq), &
               adj_tm(adj_max,lq), adj_qm(adj_max,lq), adj_t1(adj_max,lq), adj_q1(adj_max,lq))
      adj_kcall = -1 ; adj_kk = 0 ; adj_flag = 0
      adj_psp = 0. ; adj_t0 = 0. ; adj_q0 = 0.
      adj_tm = 0. ; adj_qm = 0. ; adj_t1 = 0. ; adj_q1 = 0.
      adj_n = 0
!
      ztmst=dt
""")

# record the entry state of every cuadjtqn call
sub("""      zqmax=0.5

!     2.           calculate condensation and adjust t and q accordingly
!                  -----------------------------------------------------
""", """      zqmax=0.5

      if (adj_n < adj_max) then
        adj_n = adj_n + 1
        adj_kcall(adj_n) = kcall
        adj_kk(adj_n)    = kk
        do jl = 1,klon
          adj_psp(adj_n,jl) = psp(jl)
          adj_t0(adj_n,jl)  = pt(jl,kk)
          adj_q0(adj_n,jl)  = pq(jl,kk)
          adj_tm(adj_n,jl)  = pt(jl,kk)
          adj_qm(adj_n,jl)  = pq(jl,kk)
          if (kcall == 0) then
            adj_flag(adj_n,jl) = 1
          else if (ldflag(jl)) then
            adj_flag(adj_n,jl) = 1
          end if
        end do
      end if

!     2.           calculate condensation and adjust t and q accordingly
!                  -----------------------------------------------------
""")

# the mid-state, once per kcall branch, right after the FIRST Newton step
sub("""            pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk) - zcond
            zl = 1./(pt(jl,kk)-c4les)
""", """            pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk) - zcond
            if (adj_n > 0) then
              adj_tm(adj_n,jl) = pt(jl,kk) ; adj_qm(adj_n,jl) = pq(jl,kk)
            end if
            zl = 1./(pt(jl,kk)-c4les)
""")
sub("""          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          zqsat = foeewm(pt(jl,kk))*zqp
""", """          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          if (adj_n > 0) then
            adj_tm(adj_n,jl) = pt(jl,kk) ; adj_qm(adj_n,jl) = pq(jl,kk)
          end if
          zqsat = foeewm(pt(jl,kk))*zqp
""")
sub("""        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
        zqsat = foeewm(pt(jl,kk))*zqp
""", """        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
        if (adj_n > 0) then
          adj_tm(adj_n,jl) = pt(jl,kk) ; adj_qm(adj_n,jl) = pq(jl,kk)
        end if
        zqsat = foeewm(pt(jl,kk))*zqp
""")

# the exit state
sub("""      end if

      return
      end subroutine cuadjtqn""", """      end if

      if (adj_n > 0) then
        do jl = 1,klon
          adj_t1(adj_n,jl) = pt(jl,kk) ; adj_q1(adj_n,jl) = pq(jl,kk)
        end do
      end if

      return
      end subroutine cuadjtqn""")

# --------------------------------------------------- the cuinin stage boundary
# cuinin's outputs are cumastrn locals that cuascn and the rest overwrite many
# times, so they have to be read where they are produced.
sub("""     &     plude,    ilab)

!----------------------------------
!*    3.0   cloud base calculations
!----------------------------------""", """     &     plude,    ilab)

      if (esm_dump_inner) then
        call esm_dump_var('ci_pten', pten); call esm_dump_var('ci_pqen', pqen)
        call esm_dump_var('ci_pqsen', pqsen(1:klon*klev))
        call esm_dump_var('ci_puen', puen); call esm_dump_var('ci_pven', pven)
        call esm_dump_var('ci_pverv', pverv); call esm_dump_var('ci_pgeo', pgeo)
        call esm_dump_var('ci_paph', paph); call esm_dump_var('ci_pgeoh', zgeoh)
        call esm_dump_var('ci_ztenh', ztenh); call esm_dump_var('ci_zqenh', zqenh)
        call esm_dump_var('ci_zqsenh', zqsenh); call esm_dump_var('ci_ilwmin', ilwmin)
        call esm_dump_var('ci_ptu', ptu); call esm_dump_var('ci_pqu', pqu)
        call esm_dump_var('ci_ztd', ztd); call esm_dump_var('ci_zqd', zqd)
        call esm_dump_var('ci_plu', plu)
        call esm_dump_var('ci_zuu', zuu); call esm_dump_var('ci_zvu', zvu)
        call esm_dump_var('ci_zud', zud); call esm_dump_var('ci_zvd', zvd)
        call esm_dump_var('ci_ilab', ilab)
      endif

!----------------------------------
!*    3.0   cloud base calculations
!----------------------------------""")

# ------------------------------------------- the existing inner cumastrn dump
sub("""     &     scale_fac, scale_fac2)
!
""", """     &     scale_fac, scale_fac2)
#include "ntiedtke_esm_inner.inc"
!
""")

open(out, "w").write(s)
assert 'ntiedtke_esm_inner.inc' in s
print(f"wrote {out}")
