#! /bin/csh -f

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.test
set subfile = cice.submit

set nthrds = ${CICE_NTHRDS}

#==========================================

# Print information about this test to stdout
echo ""
echo "Test name    : ${CICE_TESTNAME}"
echo "Test case dir: ${CICE_CASEDIR}"
echo "Test         : ${CICE_TEST}"
echo "BaseGen      : ${CICE_BASEGEN}"
echo "BaseCom      : ${CICE_BASECOM}"

# Create test script that runs cice.run, and validates
#==========================================

# Write the batch code into the job file
${CICE_SCRIPTS}/cice.batch.csh ${jobfile}
if ($? != 0) then
  exit -1
endif

cat >> ${jobfile} << EOF2

cd ${CICE_CASEDIR}
source ./cice.settings || exit 2
source ./env.\${CICE_MACHINE} || exit 2

# Check to see if executable exists in CICE_RUNDIR
if ( ! -f ${CICE_RUNDIR}/cice ) then
  echo "cice executable does not exist in ${CICE_RUNDIR}  "
  echo "Please run cice.build before this test."
  exit 99
endif

EOF2

if ( -f ${CICE_SCRIPTS}/tests/test_${CICE_TEST}.script) then
  echo "${0:t} using test_${CICE_TEST}.script"
  cat >> ${jobfile} < ${CICE_SCRIPTS}/tests/test_${CICE_TEST}.script
else
  echo "${0:t} ERROR: ${CICE_SCRIPTS}tests/test_${CICE_TEST}.script not found"
  exit -1
endif
cat >> ${jobfile} < ${CICE_SCRIPTS}/tests/baseline.script

chmod +x ${jobfile}

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
exit 0
