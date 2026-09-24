#!/usr/bin/env python3
"""Is New Tiedtke's trigger a FROZEN STATE in the N107 sense?  No -- and this
measures the dt dependence that IS there, so the absence is evidence, not a claim.

cu_ntiedtke_run computes qsat, calls cumastrn ONCE, and applies the returned
tendency as x + tendency*dt.  cumastrn itself runs each stage exactly once in a
straight line (cuinin -> cutypen -> cuascn -> cudlfsn -> cuddrafn -> closure ->
cuflxn -> cudtdqn -> cududvn); the closure RESCALES the already-computed mass
fluxes rather than re-running the ascent.  So there is no internal integration
during which the trigger could go stale relative to a state moving underneath it.

What dt DOES do is enter the scheme as a PARAMETER: ztmst floors the CAPE
adjustment time scale (ztauc = max(ztmst, ztauc)) and scales the nonequil
closure.  That makes the returned tendency dt-dependent, which this probe
measures, and it is the same situation as GWDO's deltim in the wind-reversal
limiter -- visible in the argument list, not hidden inside an integration.
"""
import os, subprocess, sys, numpy as np
sys.path.insert(0,'/projects/illinois/eng/cee/ctessum/ctessum/code/EqWeFiC/tools')
import esm_dump
K='wrf-chem170/kernels/build/ntiedtke_driver'
OUT='chem170_scm/nt_dtprobe'; os.makedirs(OUT,exist_ok=True)
DTS=[60.0, 30.0, 6.0, 0.6, 0.06]
for st in (1, 601, 2101):
    fl=f'chem170_scm/nt_flat/{st}.flat'
    ref=None
    print(f'--- step {st}')
    for dt in DTS:
        e=dict(os.environ, ESM_DT=repr(dt))
        o=f'{OUT}/{st}_dt{dt}.json'
        subprocess.run([K,fl,o],check=True,capture_output=True,env=e)
        b=esm_dump.load(o)['vars']
        r=np.asarray(b['rthcuten'],float)
        kt=np.asarray(b['ktype']).reshape(-1); cb=np.asarray(b['kcbot']).reshape(-1)
        if ref is None: ref=r; amp=np.abs(r).max(); kt0=kt.copy(); cb0=cb.copy()
        print(f'   dt={dt:7.2f}  ktype={kt}  kcbot={cb}  max|rthcuten|={np.abs(r).max():.4e}  '
              f'|r-r(60s)|/amp={np.abs(r-ref).max()/max(amp,1e-300):.3e}  '
              f'trigger same as 60 s: {bool((kt==kt0).all() and (cb==cb0).all())}')
