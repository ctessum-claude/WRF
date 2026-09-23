#!/usr/bin/env python3
"""Splice the COLD-SEASON EqWeFiC instrumentation includes into the build-time
copy of module_sf_noahlsm.F that kernels/Makefile already makes for the three
warm-season splices.  Each anchor must occur exactly once or this fails, so a
change in the WRF source cannot silently drop a hook (the same guarantee the
Makefile's grep guard gives the sed splices).  The WRF source is untouched:
only the copy under kernels/build is modified.
    noah_cold_splice.py <in.F> <out.F>
"""
import sys

# (anchor line, include file, "after" | "before")
SPLICES = [
    ('#include "noah_esm_state.inc"',
     "noah_esm_cold_state.inc", "after"),
    ("  END SUBROUTINE FRH2O",
     "noah_esm_frh2o.inc", "before"),
    ("      XH2O = SH2O + QTOT * DT / (DH2O * HLICE * DZ)",
     "noah_esm_snksrc_raw.inc", "after"),
    ("      SH2O = XH2O",
     "noah_esm_snksrc.inc", "after"),
    ("            EX = FLX3*0.001/ LSUBF",
     "noah_esm_snopac_exraw.inc", "after"),
    ("      SICE = SMC (1) - SH2O (1)",
     "noah_esm_hrt_sice1.inc", "after"),
    ("         SICE = SMC (K) - SH2O (K)",
     "noah_esm_hrt_sicek.inc", "after"),
    ("         TSURF = (YY + (ZZ1-1) * STC (1)) / ZZ1",
     "noah_esm_hrt_tsurf.inc", "after"),
    ("         CALL TBND (STC (1),STC (2),ZSOIL,ZBOT,1,NSOIL,TBK)",
     "noah_esm_hrt_tbnd1.inc", "after"),
    ("               CALL TBND (STC (K),STC (K +1),ZSOIL,ZBOT,K,NSOIL,TBK1)",
     "noah_esm_hrt_tbndk.inc", "after"),
    ("      T12 = (SFCTMP + T12A + T12B) / DENOM",
     "noah_esm_snopac_t12.inc", "after"),
    ("         T1 = T12",
     "noah_esm_snopac_cold.inc", "after"),
    ("         T1 = TFREEZ * max(0.01,SNCOVR ** SNOEXP) + T12 * "
     "(1.0- max(0.01,SNCOVR ** SNOEXP))",
     "noah_esm_snopac_melt.inc", "after"),
    ("  END SUBROUTINE SNOPAC",
     "noah_esm_snopac_end.inc", "before"),
    ("      TSOILC = TSOIL -273.15",
     "noah_esm_snowpack_in.inc", "after"),
    ("      PEXP = PEXP + 1.",
     "noah_esm_snowpack_pexp.inc", "before"),
    ("      IF (DSX < 0.05) DSX = 0.05",
     "noah_esm_snowpack_dsx.inc", "after"),
    ("         SNDENS = SNDENS * (1. - DW) + DW",
     "noah_esm_snowpack_wet.inc", "after"),
    ("  END SUBROUTINE SNOWPACK",
     "noah_esm_snowpack_out.inc", "before"),
]

def main(src, dst):
    lines = open(src).read().split("\n")
    for anchor, inc, where in SPLICES:
        hits = [i for i, l in enumerate(lines) if l.rstrip() == anchor]
        if len(hits) != 1:
            sys.exit("noah_cold_splice: anchor for %s matched %d lines, expected 1:\n  %r"
                     % (inc, len(hits), anchor))
        i = hits[0]
        lines.insert(i + 1 if where == "after" else i, '#include "%s"' % inc)
    open(dst, "w").write("\n".join(lines))
    print("noah_cold_splice: %d cold-path includes spliced" % len(SPLICES))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
