
# BFB Compare runs
#
# Sourced from baseline.script at the end of a test job, and again from
# bfbcomp_redo.csh after the whole suite has drained.  The in-job pass can
# only wait a bounded time for its comparison target to finish (a compute
# node sitting idle in a sleep loop costs allocation), so when the target is
# still running the in-job comparison reads a partial log and reports
# different-data or usage-error.  The redo pass runs on the login node once
# every job is done and overwrites that verdict, which is why the bfbcomp
# line in test_output is stripped and rewritten rather than appended to.
#
# ICE_BFBCOMP_REDO    set by the redo pass; skips the wait, everything is done
# ICE_BFBCOMP_MAXWAIT number of 60s polls to wait in-job (default 5)

if (${ICE_BFBCOMP} != ${ICE_SPVAL}) then

  echo "PEND ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP}" >> ${ICE_CASEDIR}/test_output
  if (${ICE_BFBCOMP} != ${ICE_TESTNAME} && ! $?ICE_BFBCOMP_REDO) then
    set maxwait = 5
    if ($?ICE_BFBCOMP_MAXWAIT) set maxwait = ${ICE_BFBCOMP_MAXWAIT}
    # Check if the baseline job is complete
    set job = `grep " ${ICE_BFBCOMP} " ../suite.jobs | sed  's|^[^0-9]*\([0-9]*\).*$|\1|g' | head -1`
    echo "checking on Job $job"
    set qstatjob = 1
    set cnt = 0
    if (${job} =~ [0-9]*) then
      while ($qstatjob)
        # historical avoids completed jobs on PBS (-x) and extra $job avoids superfluous header lines
        set qstatus = `${ICE_MACHINE_QSTAT} $job | grep -iv " historical " | grep $job | wc -l`
#        ${ICE_MACHINE_QSTAT} $job
#        echo $job $qstatus
        if ($qstatus == 0) then
          echo "Job $job completed"
          set qstatjob = 0
        else
          @ cnt = $cnt + 1
          echo "Waiting for $job to complete $cnt"
          sleep 60   # Sleep for 1 minute, so as not to overwhelm the queue manager
          if ($cnt > $maxwait) then
            # The comparison below will read a partial log and is not
            # meaningful; bfbcomp_redo.csh redoes it after the suite drains.
            echo "No longer waiting for $job to complete"
            set qstatjob = 0 # Abandon check after cnt sleep 60 checks
          endif
        endif
#        echo $qstatjob
      end
    endif
  endif

  if (${ICE_BFBTYPE} == "log") then
    set test_file = `ls -1t ${ICE_RUNDIR}/cice.runlog* | head -1`
    set base_file = `ls -1t ${ICE_RUNDIR}/../${ICE_BFBCOMP}.${ICE_TESTID}/cice.runlog* | head -1`

    echo ""
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

    echo ""
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

    echo ""
    echo "bfb Restart Compare Mode:"
    echo "base_dir: ${base_dir}"
    echo "test_dir: ${test_dir}"

    ${ICE_CASEDIR}/casescripts/comparebfb.csh ${base_dir} ${test_dir}
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
    echo ""
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

    echo ""
    echo "bfb Restart Compare Mode:"
    echo "base_dir: ${base_dir}"
    echo "test_dir: ${test_dir}"

    ${ICE_CASEDIR}/casescripts/comparebfb.csh ${base_dir} ${test_dir}
    set bfbstatus = $status
  endif

  mv -f ${ICE_CASEDIR}/test_output ${ICE_CASEDIR}/test_output.prev
  cat ${ICE_CASEDIR}/test_output.prev | grep -iv "${ICE_TESTNAME} bfbcomp " >! ${ICE_CASEDIR}/test_output
  rm -f ${ICE_CASEDIR}/test_output.prev
  if (${bfbstatus} == 0) then
    echo "PASS ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP}" >> ${ICE_CASEDIR}/test_output
    echo "bfbcomp baseline and test dataset pass"
  else if (${bfbstatus} == 1) then
    echo "FAIL ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP} different-data" >> ${ICE_CASEDIR}/test_output
    echo "bfbcomp baseline and test dataset fail"
  else if (${bfbstatus} == 2) then
    echo "MISS ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP} missing-data" >> ${ICE_CASEDIR}/test_output
    echo "Missing data"
  else
    echo "FAIL ${ICE_TESTNAME} bfbcomp ${ICE_BFBCOMP} usage-error" >> ${ICE_CASEDIR}/test_output
    echo "bfbcomp and test dataset usage error"
  endif
endif
