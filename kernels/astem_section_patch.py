#!/usr/bin/env python3
"""Build-time patch for kernels/astem_section_driver (EqWeFiC per-section ASTEM replay).

    astem_section_patch.py module_mosaic_therm.f90

Applied to the preprocessed COPY in kernels/build_astem_sec only; chem/module_mosaic_therm.F
and the model are untouched.  ASTEM's flux routines make state edits inside a section's own
flux call -- the absorb/degas_tiny seeds, degas_solid_nh4no3/nh4cl, equilibrate_acids -- which
move gas and aerosol at once.  The patch renames each to <name>_esmorig and adds a wrapper of
the same name that skips it when the module flag esm_sec_noedit is set.  The driver forms
every section's flux twice, with the edits (what the model does) and without; where the two
fluxes are bit-identical the flux was formed at the section's state and did not depend on the
edit.  The no-edit run is a DETECTOR only: nothing it computes is used as a reference unless
it equals the model's own result.
"""
import re
import sys

EDITS = ['absorb_tiny_nh4no3', 'absorb_tiny_nh4cl', 'absorb_tiny_hno3', 'absorb_tiny_hcl',
         'degas_tiny_nh4no3', 'degas_tiny_nh4cl', 'degas_solid_nh4no3', 'degas_solid_nh4cl',
         'equilibrate_acids']
p = sys.argv[1]
s = open(p).read()
for e in EDITS:
    s, n1 = re.subn(rf'(?im)^(\s*subroutine\s+){e}(\s*\(\s*ibin\s*\))', rf'\g<1>{e}_esmorig\2', s)
    s, n2 = re.subn(rf'(?im)^(\s*end\s+subroutine\s+){e}\b', rf'\g<1>{e}_esmorig', s)
    assert n1 == 1 and n2 == 1, (e, n1, n2)
# module flag before the module's first CONTAINS
m = re.search(r'(?im)^\s*contains\s*$', s)
s = s[:m.start()] + '\n      logical :: esm_sec_noedit = .false.\n' + s[m.start():]
wr = ''.join(f"""
      subroutine {e}(ibin)
      integer ibin
      if (esm_sec_noedit) return
      call {e}_esmorig(ibin)
      end subroutine {e}
""" for e in EDITS)
m = list(re.finditer(r'(?im)^\s*end\s+module\b', s))[-1]
s = s[:m.start()] + wr + s[m.start():]
open(p, 'w').write(s)
print('patched', len(EDITS), 'state edits')
