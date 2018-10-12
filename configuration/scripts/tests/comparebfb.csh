#!/bin/csh -f

# Compare the restart files or log files between the 2 cases
#-----------------------------------------------------------

# usage: comparebfb.script base_dir test_dir
#
#    base_dir: directory of either restart files or log files for the base 
#              simulation.
#    test_dir: directory of either restart files or log files for the test
#              simulation.  
#
#    To run a `restart` test, only 1 directory is passed

# Check to see if this is a restart case (i.e., only 1 directory is passed)
if ( $#argv == 1 ) then
  set restart = 1
  set compare = 0
  set base_dir = $argv[1]
else if ( $#argv == 2 ) then
  if ( -f $argv[1] ) then
    # Files were passed
    set compare = 1
    set restart = 0
    set base_file = $argv[1]
    set test_file = $argv[2]
  else
    # Directories were passed
    set compare = 0
    set restart = 0
    set base_dir = $argv[1]
    set test_dir = $argv[2]
  endif
else
  echo "Error in ${0}"
  echo "Usage: ${0} <base_dir> <test_dir>"
  echo "<test_dir> only included for non-restart tests"
endif

if ( $compare == 0 ) then
  # Check to see if the base directory includes runlogs, or restart files
  set numfiles = `find $base_dir -maxdepth 1 -name 'cice.runlog*' | wc -l`
  if ( $numfiles > 0 ) then
    # Compare log files
    set logcmp = 1
  else
    # Compare restart files
    set numfiles = `find $base_dir -maxdepth 1 -name '*.nc' | wc -l`
    if ( $numfiles > 0 ) then
      # Compare netcdf files
      set binary = 0
    else
      # Compare binary files
      set binary = 1
    endif
    set logcmp = 0
  endif
else
  set logcmp = 0
  set binary = 0
endif

if ( $logcmp == 1 ) then
  # Compare the diagnostic output in the log files
  # ---------------------
  echo "Performing comparison between log files:"
  echo ""
  echo ""

  if ( $restart == 1 ) then
    # This is a restart test.  Grab the base and test log files from the same directory
    set base_log = `ls -t1 $base_dir/cice.runlog* | head -2 | tail -1`
    set test_log = `ls -t1 $base_dir/cice.runlog* | head -1`
  else
    # This is a test comparing 2 separate directories
    set base_log = `ls -t1 $base_dir/cice.runlog* | head -1`
    set test_log = `ls -t1 $test_dir/cice.runlog* | head -1`
  endif
  echo "base: $base_log"
  echo "test: $test_log"

  set base_out = `tac $base_log | awk 'BEGIN{found=1;} /istep1:/ {found=0} {if (found) print}' | tac | grep '= ' | grep -v 'min, max, sum' | tr '\n' ','`
  set test_out = `tac $test_log | awk 'BEGIN{found=1;} /istep1:/ {found=0} {if (found) print}' | tac | grep '= ' | grep -v 'min, max, sum'`

  # Ensure that there is diagnostic output
  if ( ${#base_out} < 10 || ${#test_out} < 10 ) then
    echo "No diagnostic output for comparison"
    exit 1
  endif

  # Replace all asterisks (*) with a period (!) as workaround for errors
  #   encountered looping through words with asterisks in csh
  set base_out = `echo "$base_out" | sed 's/\*/./g'`
  set test_out = `echo "$test_out" | sed 's/\*/./g'`

  set failure = 0
  # Loop through each line of diagnostic output and check for differences
  foreach line ( "`echo '$base_out' | tr ',' '\n'`" )
    foreach word ( $line )
      if ( "$word" != "$test_out[1]" ) then
        # Print the difference to the log
        echo "Difference in:"
        echo "$line"
        echo "Base value: $word"
        echo "Test value: $test_out[1]"
        set failure = 1
      endif
      shift test_out
    end
  end

  if ( $failure == 0 ) then
    exit 0
    #echo "PASS ${ICE_TESTNAME} log " >> ${ICE_CASEDIR}/test_output
  else
    exit 1
    #echo "FAIL ${ICE_TESTNAME} log " >> ${ICE_CASEDIR}/test_output
  endif
else if ( $compare == 0 ) then
  echo "Exact Restart Comparison Mode:"
  if ( $binary == 1 ) then
    if ( $restart == 1 ) then
      # Restart case.  Compare restart files (iced.*) to base files (base_iced.*)
      set end_date = `ls -t1 $base_dir | head -1 | awk -F'.' '{print $NF}'`
      set failure = 0
      foreach file (${base_dir}/base_*${end_date})
        echo "Performing binary comparison between files:"
        set base_data = `echo $file | awk -F'/' '{print $NF}'`
        set test_data = `echo $file | awk -F'/' '{print $NF}' | cut -c6-`
        echo "base: $base_data"
        echo "test: $test_data"
        cmp -s $base_data $test_data
        if ( $? == 1 ) then
          set failure = 1
          echo "Failure in data comparison"
          break
        endif
      end
    else
      # bfbcmp case.  Compare restart files (iced.*) between the 2 directories
      set end_date = `ls -t1 $base_dir | head -1 | awk -F'.' '{print $NF}'`
      set failure = 0
      foreach file (${base_dir}/*${end_date})
        echo "Performing binary comparison between files:"
        set base_data = $base_dir/`echo $file | awk -F'/' '{print $NF}'`
        set test_data = $test_dir/`echo $file | awk -F'/' '{print $NF}'`
        echo "base: $base_data"
        echo "test: $test_data"
        cmp -s $base_data $test_data
        if ( $? == 1 ) then
          set failure = 1
          echo "Failure in data comparison"
          break
        endif
      end
    endif
    if ( $failure == 0 ) then
      exit 0
      #echo "PASS ${ICE_TESTNAME} test " >> ${ICE_CASEDIR}/test_output
    else
      exit 1
      #echo "FAIL ${ICE_TESTNAME} test " >> ${ICE_CASEDIR}/test_output
    endif
  else
    echo "Performing binary comparison between files:"
    if ( $restart == 1 ) then
      # This is a restart test.  Grab the restart files from the same directory
      set base_file = $base_dir/`ls -t1 $base_dir | head -2 | tail -1`
      set test_file = $base_dir/`ls -t1 $base_dir | head -1`
    else
      # This is a test comparing 2 separate directories
      set base_file = $base_dir/`ls -t1 $base_dir | head -1`
      set test_file = $test_dir/`ls -t1 $test_dir | head -1`
    endif
    echo "base: $base_file"
    echo "test: $test_file"
    if ( { cmp -s $test_file $base_file } ) then
      exit 0
      #echo "PASS ${ICE_TESTNAME} test " >> ${ICE_CASEDIR}/test_output
    else
      exit 1
      #echo "FAIL ${ICE_TESTNAME} test " >> ${ICE_CASEDIR}/test_output
    endif
  endif
else if ( $compare == 1 ) then
  # Compare restart files for differing cases (bfbcomp)
  echo ""
  echo "BFB Compare Mode:"
  if ( "$base_file" =~ *.nc && "$test_file" =~ *.nc ) then
    echo "Comparing netcdf files"
  else if ( "$base_file" !=~ *.nc && "$test_file" !=~ *.nc ) then
    echo "Comparing binary files"
  else
    echo "${0}: A comparison cannot be performed between netcdf and binary files."
    exit 1
  endif
  echo "base: $base_file"
  echo "test: $test_file"
  if ( { cmp -s $test_file $base_file } ) then
    exit 0
  else
    exit 1
  endif
else
  echo "${0}: script failure"
endif  # if logcmp
