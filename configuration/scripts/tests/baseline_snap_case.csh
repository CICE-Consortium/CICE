#!/bin/csh -f

# Copy one finished case into a baseline, replacing whatever is there.
# Run from inside the case directory; see baseline_snap.csh.

if ($#argv != 1) then
  echo "${0}: ERROR, usage: ${0} <baseline_name>"
  exit 1
endif

set basename = $argv[1]

if (! -e ./cice.settings) then
  echo "${0}: ERROR, no cice.settings here"
  exit 1
endif

source ./cice.settings || exit 2

if (! -e ${ICE_CASEDIR}/test_output) exit 0

# Report through test_output, the same way --bgen does, so results.csh sees a
# snapshot exactly as it sees a generating run and counts a bad one under
# failgen.  Replace any previous generate line rather than appending, so this
# stays re-runnable: see the same pattern in bfbcomp.csh.
mv -f ${ICE_CASEDIR}/test_output ${ICE_CASEDIR}/test_output.prev
cat ${ICE_CASEDIR}/test_output.prev | grep -iv "${ICE_TESTNAME} generate " >! ${ICE_CASEDIR}/test_output
rm -f ${ICE_CASEDIR}/test_output.prev

# Only snapshot a case whose own run completed.  A failed or pending run
# leaves a partial directory, and a baseline built from one is worse than no
# baseline at all: every later comparison against it is quietly wrong and
# reads as a code regression.
set ok = `grep -c "^PASS ${ICE_TESTNAME} run" ${ICE_CASEDIR}/test_output`
if (${ok} == 0) then
  echo "FAIL ${ICE_TESTNAME} generate ${basename} run-did-not-pass" >> ${ICE_CASEDIR}/test_output
  echo "SKIP ${ICE_TESTNAME} (run did not pass)"
  exit 0
endif

set baseline_dir = ${ICE_BASELINE}/${basename}/${ICE_TESTNAME}

# Replace: the point of this script is that it can be re-run.  cice.setup's
# --bgen refuses an existing directory, deliberately, which is why it cannot.
rm -rf ${baseline_dir}
mkdir -p ${baseline_dir}
cp -p -r ${ICE_RUNDIR}/* ${baseline_dir}/
set cpstatus = $status

if (${cpstatus} != 0) then
  echo "FAIL ${ICE_TESTNAME} generate ${basename} copy-failed" >> ${ICE_CASEDIR}/test_output
  echo "FAIL ${ICE_TESTNAME} (copy failed)"
  exit 1
endif

echo "PASS ${ICE_TESTNAME} generate ${basename}" >> ${ICE_CASEDIR}/test_output
echo "SNAP ${ICE_TESTNAME}"

exit 0
