#!/bin/csh -f

# Compare the binary files 
#-----------------------------------------------------------

# usage: comparebfb.csh base_file test_file
#     does binary diff of two files
# usage: comparebfb.csh base_dir test_dir
#     looks for base_iced and iced binary files for comparison
# usage: comparebfb.csh base_dir
#     looks for iced binary files in both directories for comparison
#
# Return Codes (depends on quality of error checking)
#  0 = pass
#  1 = fail
#  2 = missing data
#  9 = error

set restart = 0
set filearg = 0

if ( $#argv == 1 ) then
  # One Dir was passed, assume it's for restart test
  set restart = 1
  set base_dir = $argv[1]
  set test_dir = $argv[1]
else if ( $#argv == 2 ) then
  if ( -f $argv[1] ) then
    # Two Files were passed
    set filearg = 1
    set base_data = $argv[1]
    set test_data = $argv[2]
  else
    # Two Directories were passed, assume it's to compare restart files
    set base_dir = $argv[1]
    set test_dir = $argv[2]
  endif
else
  echo "Error in ${0}"
  echo "Usage: ${0} <base_file> <test_file>"
  echo "   does binary diff of two files"
  echo "Usage: ${0} <restart_dir>"
  echo "   looks for base_iced and iced binary files for comparison"
  echo "Usage: ${0} <base_dir> <test_dir>"
  echo "   looks for iced binary files in both directories for comparison"
  exit 9
endif

set failure = 0

if ($filearg == 1) then
  echo "base_data: $base_data"
  echo "test_data: $test_data"
  if ( -e ${base_data} && -e ${test_data}) then
    if ( { cmp -s ${base_data} ${test_data} } ) then
      echo "  compare OK"
    else
      set failure = 1
      echo "  compare FAIL"
    endif
  else
    set failure = 2
    echo "  missing data"
  endif

else
  set end_date = `ls -t1 $test_dir | head -1 | sed 's|^.*\([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]-[0-9][0-9][0-9][0-9][0-9]\).*|\1|'`
  # set nonomatch so that if foreach does not find anything it does not end the script
  set nonomatch
  foreach test_data (${test_dir}/iced*${end_date}*)
    set test_file = "${test_data:t}"
    if ($restart == 1) then
      set base_data = ${base_dir}/base_${test_file}
    else
      set base_data = ${base_dir}/${test_file}
    endif
    echo "base_data: ${base_data}"
    echo "test_data: ${test_data}"
    if ( -e ${base_data} && -e ${test_data}) then
      if ( { cmp -s ${base_data} ${test_data} } ) then
        echo "  compare OK"
      else
        set failure = 1
        echo "  compare FAIL"
      endif
    else
      if ($failure == 0) set failure = 2
      echo "  missing data"
    endif
  end
  unset nonomatch
endif

exit ${failure}

