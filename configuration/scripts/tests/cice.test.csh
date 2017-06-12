#! /bin/csh -f

echo ${0}

# Define this as a test case in the cice.settings file
cat >> ./cice.settings << EOF

# Add variable defining this case_dir as a test
setenv CICE_TEST true
EOF

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.test
set subfile = cice.submit

set nthrds = ${CICE_NTHRDS}

#==========================================

# Create test script that runs cice.run, and validates
#==========================================

# Write the batch code into the job file
$1/cice.batch.csh ${jobfile}

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

EOF2

if ($2 != "") then
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
