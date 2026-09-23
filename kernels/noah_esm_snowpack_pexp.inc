! The polynomial series SUM_{k=1..IPOL} x^k/(k+1)! with x = BFAC*ESDCX, taken
! BEFORE WRF adds the 1 to it.  (e**x - 1)/x differs from 1 by O(x), and x is
! proportional to DTSEC, so at any reference step small enough to be a dt -> 0
! limit the stored PEXP is 1 + O(1e-6) and the compaction signal is gone from
! it.  Publishing the sum itself keeps that signal at full precision: it is the
! whole of the compaction, and it is what the .esm component can be pinned to.
! Same argument as FORTRAN_BUGS N97/N98 one routine over.
      ESM_SNP_PEXP_M1 = PEXP
