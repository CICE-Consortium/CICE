#!/bin/csh -f

# This script is passed the test name, and whether or not bfbcomp is to be performed.
# It then builds the CTestTestfile.cmake file per the inputs

# $1 = $tsdir (the directory where CTestTestfile.cmake is to be put)
# $2 = $testname_noid (the name of the test)
# $3 = $bfbcomp 
# $4 = $spval (used to test for bfbcomp)

# Every test needs to have the "build" phase
echo "add_test($2_build grep "\""PASS.*$2 .*build"\"" results.log)" >> $1/CTestTestfile.cmake

# If this is a restart test, add the 'run-initial' and 'run-restart' tests
if ( "$2" =~ *"_restart_"* ) then
  echo "add_test($2_run_initial grep "\""PASS.*$2 .*run-initial"\"" results.log)" >> $1/CTestTestfile.cmake
  echo "add_test($2_run_restart grep "\""PASS.*$2 .*run-restart"\"" results.log)" >> $1/CTestTestfile.cmake
  echo "add_test($2_exact_restart grep "\""PASS.*$2 .*exact-restart"\"" results.log)" >> $1/CTestTestfile.cmake
else
  echo "add_test($2_run grep "\""PASS.*$2 .*run"\"" results.log)" >> $1/CTestTestfile.cmake
  echo "add_test($2_compare grep "\""PASS.*$2 .*compare"\"" results.log)" >> $1/CTestTestfile.cmake
  # Check for bfbcomp in sets
  if ( $3 != $4 ) then
    echo "add_test($2_bfbcomp grep "\""PASS.*$2 .*bfbcomp"\"" results.log)" >> $1/CTestTestfile.cmake
  endif
endif

