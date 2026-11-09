#!/usr/bin/env -S csh -f

# Run the bfbcomp comparison for one case versus another case
# Run from inside the case directory

echo ""
echo "Case bfb Compare Mode with bfbcomp_update.csh:"

if (! -e ./cice.settings) then
  echo "Skipping, not a case directory"
  echo ""
  exit 1
endif

source ./cice.settings || exit 2
source ./env.${ICE_MACHCOMP} -nomodules || exit 2

# Check status of the case

echo "CASE: ${ICE_CASENAME}"

if (! -e ${ICE_CASEDIR}/test_output) then
  echo "Skipping, test_output missing, not a test"
  echo ""
  exit 0
endif

if (${ICE_BFBCOMP} == ${ICE_SPVAL}) then
  echo "Skipping, not a bfbcomp test"
  echo ""
  exit 0
endif

set ok = `grep -c "^PASS ${ICE_TESTNAME} run" ${ICE_CASEDIR}/test_output`
if (${ok} == 0) then
  echo "Skipping, test run incomplete"
  echo ""
  exit 0
endif

set ok = `grep -c "^PASS ${ICE_TESTNAME} bfbcomp" ${ICE_CASEDIR}/test_output`
if (${ok} > 0) then
  echo "Skipping, bfbcomp already passed"
  echo ""
  exit 0
endif

# BFB Compare runs
# Checks if baseline run is still runnings
# Waits ICE_BFBWAIT number of 60s polls for baseline run
# Compares log, rest, and qc output as required
# Updates bfbcomp output in test_output

mv -f ${ICE_CASEDIR}/test_output ${ICE_CASEDIR}/test_output.prev
cat ${ICE_CASEDIR}/test_output.prev | grep -iv "${ICE_TESTNAME} bfbcomp" >! ${ICE_CASEDIR}/test_output
rm -f ${ICE_CASEDIR}/test_output.prev
echo "PEND ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP}" >> ${ICE_CASEDIR}/test_output

if (${ICE_BFBCOMP} != ${ICE_TESTNAME} ) then

  # Check if the baseline job is complete, wait if requested
  set maxwait = 5
  if ($?ICE_BFBWAIT) set maxwait = ${ICE_BFBWAIT}
  set job = `grep " ${ICE_BFBCOMP} " ../suite.jobs | sed  's|^[^0-9]*\([0-9]*\).*$|\1|g' | head -1`
  set qstatjob = 1
  set cnt = 0
  if (${job} =~ [0-9]*) then
    echo "Checking on job $job"
    while ($qstatjob)
      # historical avoids completed jobs on PBS (-x) and extra $job avoids superfluous header lines
      set qstatus = `${ICE_MACHINE_QSTAT} $job | grep -iv " historical " | grep $job | wc -l`
      if ($qstatus == 0) then
        echo "Job $job completed"
        set qstatjob = 0
      else
        @ cnt = $cnt + 1
        if ($cnt > $maxwait) then
          echo "No longer waiting for $job to complete ${cnt}"
          set qstatjob = 0 # Abandon check after cnt sleep 60 checks
        else
          echo "Waiting for $job to complete ${cnt}/${maxwait}"
          sleep 60   # Sleep for 1 minute, so as not to overwhelm the queue manager
        endif
      endif
    end
  endif
endif

if (${ICE_BFBTYPE} == "log") then
  set test_file = `ls -1t ${ICE_RUNDIR}/cice.runlog* | head -1`
  set base_file = `ls -1t ${ICE_RUNDIR}/../${ICE_BFBCOMP}.${ICE_TESTID}/cice.runlog* | head -1`

  echo "bfb Log Compare Mode:"
  echo "base_file: ${base_file}"
  echo "test_file: ${test_file}"

  if (${ICE_TARGET} == "opticep") then
     ${ICE_CASEDIR}/casescripts/comparelog.csh ${base_file} ${test_file} modcicefile
     set bfbstatus = $status
  else
     ${ICE_CASEDIR}/casescripts/comparelog.csh ${base_file} ${test_file}
     set bfbstatus = $status
  endif

else if (${ICE_BFBTYPE} == "logrest") then
  set test_file = `ls -1t ${ICE_RUNDIR}/cice.runlog* | head -1`
  set base_file = `ls -1t ${ICE_RUNDIR}/../${ICE_BFBCOMP}.${ICE_TESTID}/cice.runlog* | head -1`

  echo "bfb Log Compare Mode:"
  echo "base_file: ${base_file}"
  echo "test_file: ${test_file}"

  if (${ICE_TARGET} == "opticep") then
     ${ICE_CASEDIR}/casescripts/comparelog.csh ${base_file} ${test_file} modcicefile
     set bfbstatus = $status
  else
     ${ICE_CASEDIR}/casescripts/comparelog.csh ${base_file} ${test_file}
     set bfbstatusl = $status
  endif

  set test_dir = ${ICE_RUNDIR}/restart
  set base_dir = ${ICE_RUNDIR}/../${ICE_BFBCOMP}.${ICE_TESTID}/restart

  echo "bfb Restart Compare Mode:"
  echo "base_dir: ${base_dir}"
  echo "test_dir: ${test_dir}"

  ${ICE_CASEDIR}/casescripts/comparerest.csh ${base_dir} ${test_dir}
  set bfbstatusr = $status

  if ({$bfbstatusl} == ${bfbstatusr}) then
     set bfbstatus = ${bfbstatusl}
  else if (${bfbstatusl} == 1 || ${bfbstatusr} == 1) then
     set bfbstatus = 1
  else if ({$bfbstatusl} > ${bfbstatusr}) then
     set bfbstatus = ${bfbstatusl}
  else
     set bfbstatus = ${bfbstatusr}
  endif

  echo "bfb log, rest, combined status = ${bfbstatusl},${bfbstatusr},${bfbstatus}"

else if (${ICE_BFBTYPE} =~ qcchk*) then
  set test_dir = ${ICE_RUNDIR}
  set base_dir = ${ICE_RUNDIR}/../${ICE_BFBCOMP}.${ICE_TESTID}
  echo "qcchk Compare Mode:"
  echo "base_dir: ${base_dir}"
  echo "test_dir: ${test_dir}"
  ${ICE_SANDBOX}/configuration/scripts/tests/QC/cice.t-test.py ${base_dir} ${test_dir}
  set bfbstatus = $status
  # expecting failure, so switch value
  if (${ICE_BFBTYPE} == "qcchkf") then
     @ bfbstatus = 1 - $bfbstatus
  endif

else
  set test_dir = ${ICE_RUNDIR}/restart
  set base_dir = ${ICE_RUNDIR}/../${ICE_BFBCOMP}.${ICE_TESTID}/restart

  echo "bfb Restart Compare Mode:"
  echo "base_dir: ${base_dir}"
  echo "test_dir: ${test_dir}"

  ${ICE_CASEDIR}/casescripts/comparerest.csh ${base_dir} ${test_dir}
  set bfbstatus = $status
endif

mv -f ${ICE_CASEDIR}/test_output ${ICE_CASEDIR}/test_output.prev
cat ${ICE_CASEDIR}/test_output.prev | grep -iv "${ICE_TESTNAME} bfbcomp" >! ${ICE_CASEDIR}/test_output
rm -f ${ICE_CASEDIR}/test_output.prev
if (${bfbstatus} == 0) then
  echo "PASS ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP}" >> ${ICE_CASEDIR}/test_output
  echo "bfbcomp baseline and test dataset PASS"
else if (${bfbstatus} == 1) then
  echo "FAIL ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP} different-data" >> ${ICE_CASEDIR}/test_output
  echo "bfbcomp baseline and test dataset FAIL"
else if (${bfbstatus} == 2) then
  echo "MISS ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP} missing-data" >> ${ICE_CASEDIR}/test_output
  echo "bfbcomp baseline and test dataset missing, MISS"
else
  echo "FAIL ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP} usage-error" >> ${ICE_CASEDIR}/test_output
  echo "bfbcomp and test dataset usage ERROR"
endif

echo "bfbcomp_update.csh: done"
echo ""
exit 0
