#!/bin/csh -f

# Read in the case directories for the 4 different QC cases
set QC_DIR = "./qc_logs"
set base_case_dir = `sed -n 1p $QC_DIR/qc_dirs.txt`
set bfb_case_dir = `sed -n 2p $QC_DIR/qc_dirs.txt`
set nonbfb_case_dir = `sed -n 3p $QC_DIR/qc_dirs.txt`
set fail_case_dir = `sed -n 4p $QC_DIR/qc_dirs.txt`

# Get the history file location for each case
set rundir = `grep 'setenv ICE_RUNDIR' $base_case_dir/cice.settings | awk '{print $NF}'`
set base_histdir = "$rundir/history"

set rundir = `grep 'setenv ICE_RUNDIR' $bfb_case_dir/cice.settings | awk '{print $NF}'`
set bfb_histdir = "$rundir/history"

set rundir = `grep 'setenv ICE_RUNDIR' $nonbfb_case_dir/cice.settings | awk '{print $NF}'`
set nonbfb_histdir = "$rundir/history"

set rundir = `grep 'setenv ICE_RUNDIR' $fail_case_dir/cice.settings | awk '{print $NF}'`
set fail_histdir = "$rundir/history"

set QC_SUCCESS=0
set QC_FAIL=1

# Perform the comparisons
echo "===== Running QC tests and writing output to $QC_DIR/validate_qc.log ====="
echo "Running QC test on base and bfb directories."
echo "Expected result: PASSED"
./configuration/scripts/tests/QC/cice.t-test.py $base_histdir $bfb_histdir >& $QC_DIR/validate_qc.log
set case1="$?"
if ($case1 == $QC_SUCCESS) then
    echo "Result: PASSED"
else
    echo "Result: FAILED"
endif

echo "-----------------------------------------------"
echo "Running QC test on base and non-bfb directories."
echo "Expected result: PASSED"
./configuration/scripts/tests/QC/cice.t-test.py $base_histdir $nonbfb_histdir >>& $QC_DIR/validate_qc.log
set case2="$?"
if ($case2 == $QC_SUCCESS) then
    echo "Result: PASSED"
else
    echo "Result: FAILED"
endif

echo "-----------------------------------------------"
echo "Running QC test on base and climate-changing directories."
echo "Expected result: FAILED"
./configuration/scripts/tests/QC/cice.t-test.py $base_histdir $fail_histdir >>& $QC_DIR/validate_qc.log
set case3="$?"
if ($case3 == $QC_SUCCESS) then
    echo "Result: PASSED"
else
    echo "Result: FAILED"
endif

echo ""
echo ""

if ($case1 == $QC_SUCCESS && $case2 == $QC_SUCCESS && $case3 == $QC_FAIL ) then
    echo "QC Test has validated"
else
    echo "QC Test has failed validation"
endif
