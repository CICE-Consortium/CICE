#! /bin/csh -f

echo ${0}

# $1 = ${CICE_SCRDIR}
# $2 = $baseline
# $3 = baseline (compare) directory
# $4 = $test

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.test
set subfile = cice.submit

set nthrds = ${CICE_NTHRDS}

#==========================================

# Print information about this test to stdout
echo ""
echo "Test name: $4"
echo "Test case directory: ${CICE_CASEDIR}"
if ($4 != 'restart') then
  if ($2 == 2) then
    echo "This is a regression run"
    echo "Baseline datasets will be stored in: ${CICE_RUNDIR}"
    echo "Data will be compared to data in $3"
  else if ($2 == 1) then
    # This is a baseline generating run
    echo "This is a baseline-generating test."
    echo "Baseline datasets will be stored in: ${CICE_RUNDIR}"
  else
    # This is not a baseline-generating run
    echo "This is not a baseline-generating run."
    echo "Test data will be compared to data in $3"
  endif
endif

# Create test script that runs cice.run, and validates
#==========================================

# Write the batch code into the job file
$1/cice.batch.csh ${jobfile}
if ($? == -1) then
  exit -1
endif

if ($4 == 'restart') then
  goto restart
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
  
  set baseline_data = $3/restart/\$test_data:t

  echo "Performing binary comparison between files:"
  echo "baseline: \$baseline_data"
  echo "test:     \$test_data"
  if ( { cmp -s \$test_data \$baseline_data } ) then
    echo "PASS \${CICE_CASENAME} compare" >> ${CICE_CASEDIR}/test_output
  else
    echo "FAIL \${CICE_CASENAME} compare" >> ${CICE_CASEDIR}/test_output
  endif
EOF3
else if ($2 == 1) then
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
else
  # Regression run
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

  set baseline_data = $3/restart/\$test_data:t

  echo "Performing binary comparison between files:"
  echo "baseline: \$baseline_data"
  echo "test:     \$test_data"
  if ( { cmp -s \$test_data \$baseline_data } ) then
    echo "PASS \${CICE_CASENAME} compare" >> ${CICE_CASEDIR}/test_output
  else
    echo "FAIL \${CICE_CASENAME} compare" >> ${CICE_CASEDIR}/test_output
  endif
EOF3
endif

goto chmod

restart:
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
chmod:
  chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
