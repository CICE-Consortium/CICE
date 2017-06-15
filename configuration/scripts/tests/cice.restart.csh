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

# Run the CICE model for the 2-day simulation
./cice.run
if ( \$? != 0 ) then
  echo "FAIL \${CICE_CASENAME} 2-day-run" >> \${CICE_CASEDIR}/test_output
  exit 99
else
  echo "PASS \${CICE_CASENAME} 2-day-run" >> \${CICE_CASEDIR}/test_output
endif

# Prepend 'baseline_' to the final restart file to save for comparison
foreach file (\${CICE_RUNDIR}/restart/*)
  set test_data = \$file
end
set baseline_data = \${CICE_RUNDIR}/restart/baseline_\$test_data:t
mv \$test_data \$baseline_data

# Modify the contents of the pointer file for restart
perl -i -pe's/(\d{4})-(\d{2})-(\d{2})/sprintf("%04d-%02d-%02d",\$1,\$2,\$3-1)/e' \${CICE_RUNDIR}/ice.restart_file

# Modify the namelist for the restart simulation
cp ice_in ice_in.2-day
\${CICE_CASEDIR}/casescripts/parse_namelist.sh ice_in \${CICE_CASEDIR}/casescripts/test_nml.restart2

# Run the CICE model for the restart simulation
./cice.run
cp ice_in ice_in.restart
cp ice_in.2-day ice_in
if ( \$? != 0 ) then
  echo "FAIL \${CICE_CASENAME} restart-run" >> \${CICE_CASEDIR}/test_output
  exit 99
else
  echo "PASS \${CICE_CASENAME} restart-run" >> \${CICE_CASEDIR}/test_output
endif

EOF2

cat >> ${jobfile} << EOF3
echo "Performing binary comparison between files:"
echo "baseline: \$baseline_data"
echo "test:     \$test_data"
if ( { cmp -s \$test_data \$baseline_data } ) then
  echo "PASS \${CICE_CASENAME} compare" >> ${CICE_CASEDIR}/test_output
else
  echo "FAIL \${CICE_CASENAME} compare" >> ${CICE_CASEDIR}/test_output
endif
EOF3

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
