#!/usr/bin/env -S csh -f

# Stamp a finished test suite as a baseline, replacing whatever is there.
#
#   cd <suite_dir>
#   ./baseline_snap.csh <baseline_name>
#
# Baseline generation is only `cp -p -r` of each run directory, but it lives
# inside the test job, so cice.setup makes you choose --bgen before the suite
# runs and then refuses to overwrite what it wrote.  That forces a full rerun
# whenever a baseline needs refreshing, even when the runs you already have
# are the ones you want.
#
# This does the same copy afterwards, over a suite that has already finished,
# and replaces the whole baseline rather than failing on what exists.  So:
# run the suite once, snapshot it, and snapshot it again later if the code
# moves on -- no --bgen decision up front, and no rerun to change your mind.
#
# Cases whose own run did not pass are skipped, not copied, and recorded as
# FAIL ... generate in their test_output -- so results.csh counts them under
# failgen exactly as it would a bad --bgen run, rather than the shortfall
# being visible only here.
#
# Also works on a suite created before this script existed: the per-case
# helper is installed from the sandbox that built each case.

if ($#argv != 1) then
  echo "${0}: ERROR, usage: ${0} <baseline_name>"
  echo "  run from a test suite directory"
  exit 1
endif

set basename = $argv[1]
set logfile = baseline_snap.log
rm -f ${logfile}

foreach casedir (*)
  if (! -d ${casedir}) continue
  if (! -e ${casedir}/cice.settings) continue

  if (! -e ${casedir}/casescripts/baseline_snap_case.csh) then
    set sbox = `grep "setenv ICE_SANDBOX" ${casedir}/cice.settings | awk '{print $3}'`
    if ("${sbox}" != "") then
      if (-e ${sbox}/configuration/scripts/tests/baseline_snap_case.csh) then
        cp -f -p ${sbox}/configuration/scripts/tests/baseline_snap_case.csh ${casedir}/casescripts/
      endif
    endif
  endif
  if (! -e ${casedir}/casescripts/baseline_snap_case.csh) continue

  ( cd ${casedir} ; ./casescripts/baseline_snap_case.csh ${basename} ) >>& ${logfile}
end

set nsnap = `grep -c "^SNAP " ${logfile}`
set nskip = `grep -c "^SKIP " ${logfile}`
set nfail = `grep -c "^FAIL " ${logfile}`

echo ""
echo "baseline '${basename}': ${nsnap} snapshotted, ${nskip} skipped, ${nfail} failed"
echo "per-case list in ${logfile}; each case also records the result in its"
echo "own test_output, so results.csh reports it as it would a --bgen run"

exit 0
