#!/usr/bin/env -S csh -f

# Redo every bfbcomp comparison in this test suite.
#
# Each test runs its own bfbcomp comparison at the end of its batch job, and
# can only wait a bounded time for the test it compares against to finish.
# When that wait expires the comparison reads a log that is still being
# written and reports different-data, or one that has not yet produced any
# diagnostics and reports usage-error.  Neither verdict says anything about
# the model.  Run this once every job has completed (poll_queue.csh) and the
# comparisons are redone against finished output; results.csh calls it.
#
# Safe to run repeatedly: each comparison strips its previous bfbcomp line
# from test_output before writing the new one.
#
# Also works on a suite created before this script existed -- the per-case
# scripts it needs are copied in from the sandbox that built each case.

set logfile = bfbcomp_redo.log
rm -f ${logfile}

set cnt = 0
foreach casedir (*)
  if (! -d ${casedir}) continue
  if (! -e ${casedir}/cice.settings) continue

  # Older suites have no bfbcomp scripts in casescripts; install them
  if (! -e ${casedir}/casescripts/bfbcomp_case.csh) then
    set sbox = `grep "setenv ICE_SANDBOX" ${casedir}/cice.settings | awk '{print $3}'`
    if ("${sbox}" != "") then
      if (-e ${sbox}/configuration/scripts/tests/bfbcomp_case.csh) then
        cp -f -p ${sbox}/configuration/scripts/tests/bfbcomp.csh      ${casedir}/casescripts/
        cp -f -p ${sbox}/configuration/scripts/tests/bfbcomp_case.csh ${casedir}/casescripts/
      endif
    endif
  endif
  if (! -e ${casedir}/casescripts/bfbcomp_case.csh) continue

  echo "--- ${casedir}" >> ${logfile}
  ( cd ${casedir} ; ./casescripts/bfbcomp_case.csh ) >>& ${logfile}
  @ cnt = $cnt + 1
end

echo "bfbcomp comparisons redone for ${cnt} cases, see ${logfile}"

exit 0
