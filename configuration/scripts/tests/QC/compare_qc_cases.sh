# Read in the case directories for the 4 different QC cases
export base_case_dir=`sed '1!d' qc_dirs.txt`
export bfb_case_dir=`sed '2!d' qc_dirs.txt`
export nonbfb_case_dir=`sed '3!d' qc_dirs.txt`
export fail_case_dir=`sed '4!d' qc_dirs.txt`

# Get the history file location for each case
export rundir=`grep 'setenv ICE_RUNDIR' $base_case_dir/cice.settings | awk '{print $NF}'`
export base_histdir="$rundir/history"

export rundir=`grep 'setenv ICE_RUNDIR' $bfb_case_dir/cice.settings | awk '{print $NF}'`
export bfb_histdir="$rundir/history"

export rundir=`grep 'setenv ICE_RUNDIR' $nonbfb_case_dir/cice.settings | awk '{print $NF}'`
export nonbfb_histdir="$rundir/history"

export rundir=`grep 'setenv ICE_RUNDIR' $fail_case_dir/cice.settings | awk '{print $NF}'`
export fail_histdir="$rundir/history"

export QC_SUCCESS=0
export QC_FAIL=1

# Perform the comparisons
echo "===== Running QC tests and writing output to validate_qc.log ====="
echo "Running QC test on base and bfb directories."
echo "Expected result: PASSED"
./configuration/scripts/tests/QC/cice.t-test.py $base_histdir $bfb_histdir &> validate_qc.log
export case1="$?"
if [ $case1 == $QC_SUCCESS ]; then
    echo "Result: PASSED"
else
    echo "Result: FAILED"
fi

echo "-----------------------------------------------"
echo "Running QC test on base and non-bfb directories."
echo "Expected result: PASSED"
./configuration/scripts/tests/QC/cice.t-test.py $base_histdir $nonbfb_histdir >> validate_qc.log 2>&1
export case2="$?"
if [ $case2 == $QC_SUCCESS ]; then
    echo "Result: PASSED"
else
    echo "Result: FAILED"
fi

echo "-----------------------------------------------"
echo "Running QC test on base and climate-changing directories."
echo "Expected result: FAILED"
./configuration/scripts/tests/QC/cice.t-test.py $base_histdir $fail_histdir >> validate_qc.log 2>&1
export case3="$?"
if [ $case3 == $QC_SUCCESS ]; then
    echo "Result: PASSED"
else
    echo "Result: FAILED"
fi

echo ""
echo ""

if [ $case1 == $QC_SUCCESS ] && [ $case2 == $QC_SUCCESS ] && [ $case3 == $QC_FAIL ]; then
    echo "QC Test has validated"
else
    echo "QC Test has failed validation"
fi
