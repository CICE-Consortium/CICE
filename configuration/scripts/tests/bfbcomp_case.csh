#!/bin/csh -f

# Redo the bfbcomp comparison for one case, after the suite has drained.
# Run from inside the case directory; see bfbcomp_redo.csh.

if (! -e ./cice.settings) then
  echo "${0}: ERROR, no cice.settings here"
  exit 1
endif

source ./cice.settings || exit 2

if (${ICE_BFBCOMP} == ${ICE_SPVAL}) exit 0

# Only redo the comparison for a case whose own run completed.  If the run
# failed or is still pending, the existing status is the one worth keeping.
if (! -e ${ICE_CASEDIR}/test_output) exit 0
set ok = `grep -c "^PASS ${ICE_TESTNAME} run" ${ICE_CASEDIR}/test_output`
if (${ok} == 0) exit 0

# Tells bfbcomp.csh not to wait on the comparison target's queue job
set ICE_BFBCOMP_REDO = 1

source ${ICE_CASEDIR}/casescripts/bfbcomp.csh

exit 0
