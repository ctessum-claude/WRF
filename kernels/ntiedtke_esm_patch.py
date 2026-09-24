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

# ------------------------------- the shallow test-parcel ascent inside cutypen
# The shallow phase's arrays are reused and overwritten by the deep phase that
# follows it in the same routine, so they have to be read where they are left.
sub("""!-----------------------------------------------------------
! next, let's check the deep convection""", """      if (esm_dump_inner) then
        call esm_dump_var('sh_dh', dh); call esm_dump_var('sh_dhen', dhen)
        call esm_dump_var('sh_ptu', ptu); call esm_dump_var('sh_pqu', pqu)
        call esm_dump_var('sh_plu', plu); call esm_dump_var('sh_kup', kup)
        call esm_dump_var('sh_vptu', vptu); call esm_dump_var('sh_vten', vten)
        call esm_dump_var('sh_zbuo', zbuo); call esm_dump_var('sh_abuoy', abuoy)
        call esm_dump_var('sh_klab', klab); call esm_dump_var('sh_eta', eta)
        call esm_dump_var('sh_dz', dz); call esm_dump_var('sh_coef', coef)
        call esm_dump_var('sh_kcbot', kcbot); call esm_dump_var('sh_kctop', kctop)
        call esm_dump_var('sh_ktype', ktype); call esm_dump_var('sh_wbase', wbase)
        call esm_dump_var('sh_ptenh', ptenh); call esm_dump_var('sh_pqenh', pqenh)
        call esm_dump_var('sh_pgeoh', pgeoh); call esm_dump_var('sh_paph', paph)
        call esm_dump_var('sh_pgeo', pgeo); call esm_dump_var('sh_pten', pten)
        call esm_dump_var('sh_pqen', pqen); call esm_dump_var('sh_pap', pap)
        call esm_dump_var('sh_hfx', hfx); call esm_dump_var('sh_qfx', qfx)
      endif

!-----------------------------------------------------------
! next, let's check the deep convection""")

# ------------------- inside cuflxn: the state the precipitation sweeps start from
# Blocks 1 and 1a of cuflxn zero pdmfup/pdmfdp outside the cloud and below the
# downdraught, so the sweeps do NOT see the values cumastrn passed in.
sub("""!*    2.            calculate rain/snow fall rates                             ""","""      if (esm_dump_inner) then
        call esm_dump_var('pr_pdmfup', pdmfup); call esm_dump_var('pr_pdmfdp', pdmfdp)
        call esm_dump_var('pr_plglac', plglac); call esm_dump_var('pr_pqsen', pqsen)
        call esm_dump_var('pr_pmflxr', pmflxr); call esm_dump_var('pr_pmflxs', pmflxs)
        call esm_dump_var('pr_pdpmel', pdpmel); call esm_dump_var('pr_rhevap', rhevap)
        call esm_dump_var('pr_ktopm2', ktopm2)
      endif
!*    2.            calculate rain/snow fall rates                             """)

# ---- immediately after the cuflxn CALL, before cumastrn's downdraught-rescaling
# post-correction adds zmfuub back into the rain flux and recomputes zdmfup.
sub("""     &  ,  prain,    pmfdde_rate, pmflxr, pmflxs )

! some adjustments needed""", """     &  ,  prain,    pmfdde_rate, pmflxr, pmflxs )
      if (esm_dump_inner) then
        call esm_dump_var('cx_pmflxr', pmflxr); call esm_dump_var('cx_pmflxs', pmflxs)
        call esm_dump_var('cx_pdmfup', zdmfup); call esm_dump_var('cx_pdmfdp', zdmfdp)
        call esm_dump_var('cx_pdpmel', zdpmel); call esm_dump_var('cx_prain', prain)
        call esm_dump_var('cx_pqsen', pqsen(1:klon*klev))
      endif

! some adjustments needed""")

# ------------------------------------------------- the cuflxn stage boundary
sub("""!*    8.0      update tendencies for t and q in subroutine cudtdq""",
    """      if (esm_dump_inner) then
        call esm_dump_var('fx_pmflxr', pmflxr); call esm_dump_var('fx_pmflxs', pmflxs)
        call esm_dump_var('fx_pdmfup', zdmfup); call esm_dump_var('fx_pdmfdp', zdmfdp)
        call esm_dump_var('fx_pdpmel', zdpmel); call esm_dump_var('fx_plglac', zlglac)
        call esm_dump_var('fx_prain', prain); call esm_dump_var('fx_pmfu', pmfu)
        call esm_dump_var('fx_pmfd', pmfd); call esm_dump_var('fx_pmfus', zmfus)
        call esm_dump_var('fx_pmfds', zmfds); call esm_dump_var('fx_pmfuq', zmfuq)
        call esm_dump_var('fx_pmfdq', zmfdq); call esm_dump_var('fx_pmful', zmful)
        call esm_dump_var('fx_plude', plude); call esm_dump_var('fx_pqsen', pqsen(1:klon*klev))
        call esm_dump_var('fx_kcbot', kcbot); call esm_dump_var('fx_kctop', kctop)
        call esm_dump_var('fx_ktopm2', itopm2); call esm_dump_var('fx_pten', pten)
        call esm_dump_var('fx_pqen', pqen); call esm_dump_var('fx_paph', paph)
        call esm_dump_var('fx_pap', pap); call esm_dump_var('fx_lndj', lndj)
      endif

!*    8.0      update tendencies for t and q in subroutine cudtdq""")

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
