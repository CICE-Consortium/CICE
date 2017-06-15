#! /bin/csh -f

echo ${0}

# $1 = ${CICE_SCRDIR}
# $2 = $baseline
# $3 = baseline (compare) directory

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.test
set subfile = cice.submit

set nthrds = ${CICE_NTHRDS}

#==========================================

# Print information about this test to stdout
echo ""
echo "Test case directory: ${CICE_CASEDIR}"
if ($2 == 1) then
  # This is a baseline generating run
  echo "This is a baseline-generating test."
  echo "Baseline datasets will be stored in: ${CICE_RUNDIR}"
else
  # This is not a baseline-generating run
  echo "This is not a baseline-generating run."
  echo "Test data will be compared to data in $3"
endif

# Create test script that runs cice.run, and validates
#==========================================

# Write the batch code into the job file
$1/cice.batch.csh ${jobfile}
if ($? == -1) then
  exit -1
endif

cat >> ${jobfile} << EOF2
cd ${CICE_CASEDIR}
source ./cice.settings || exit 2
source ./env.\${CICE_MACHINE} || exit 2

# Check to see if executable exists in CICE_RUNDIR
if ( ! -f \${CICE_RUNDIR}/cice ) then
  echo "cice executable does not exist in \${CICE_RUNDIR}.  "
  echo "Please run cice.build before this test."
  exit 99
endif

# Run the CICE model
./cice.run
if ( \$? != 0 ) then
  # Run failed
  echo "FAIL \${CICE_CASENAME} run" >> \${CICE_CASEDIR}/test_output
  exit 99
else
  # Run succeeded
  echo "PASS \${CICE_CASENAME} run" >> \${CICE_CASEDIR}/test_output
endif

EOF2

if ($2 == 0) then
  # Not a baseline-generating test
cat >> ${jobfile} << EOF3
  # Get the final output filename
  foreach file (\${CICE_RUNDIR}/restart/*)
    set test_data = \$file
  end
  
  set baseline_data = $2/\$test_data:t

  echo "Performing binary comparison between files:"
  echo "baseline: \$baseline_data"
  echo "test:     \$test_data"
  if ( { cmp -s \$test_data \$baseline_data } ) then
    echo "PASS \${CICE_CASENAME} test" >> ${CICE_CASEDIR}/test_output
  else
    echo "FAIL \${CICE_CASENAME} test" >> ${CICE_CASEDIR}/test_output
  endif
EOF3
else
  # Baseline generating run
cat >> ${jobfile} << EOF3
  # Get the final output filename
  foreach file (\${CICE_RUNDIR}/restart/*)
    set test_data = \$file
  end
  
  if ( \$test_data != "" ) then
    echo "PASS \${CICE_CASENAME} generate" >> ${CICE_CASEDIR}/test_output
  else
    echo "FAIL \${CICE_CASENAME} generate" >> ${CICE_CASEDIR}/test_output
  endif
EOF3

endif

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
