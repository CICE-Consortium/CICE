#!/bin/csh -f

# Compare prognostic output in two log files
#-----------------------------------------------------------

# usage: comparelog.csh base_file test_file
#     does diff of two files
#
# Return Codes (depends on quality of error checking)
#  0 = pass
#  1 = fail
#  2 = missing data
#  9 = error

set filearg = 0
if ( $#argv == 2 ) then
  set filearg = 1
  set base_data = $argv[1]
  set test_data = $argv[2]
else
  echo "Error in ${0}"
  echo "Usage: ${0} <base_file> <test_file>"
  echo "   does diff of two files"
  exit 9
endif

set failure = 0
set base_out = "comparelog_base_out_file.log"
set test_out = "comparelog_test_out_file.log"

if ($filearg == 1) then
  echo "base_data: $base_data"
  echo "test_data: $test_data"
  if ( -f ${base_data} && -f ${test_data}) then
    if (${base_data} == ${test_data}) then
      set failure = 9
      echo "  input data are same"
    else

      touch ${base_out}
      cat ${base_data} | grep -A 99999999 istep1: | grep -e istep1: -e = | grep -iv "min, max, sum" >&! ${base_out}
      touch ${test_out}
      cat ${test_data} | grep -A 99999999 istep1: | grep -e istep1: -e = | grep -iv "min, max, sum" >&! ${test_out}

      set basenum = `cat ${base_out} | wc -l`
      set testnum = `cat ${base_out} | wc -l`
      set filediff = `diff -w ${base_out} ${test_out} | wc -l`

      if (${basenum} > 0 && ${testnum} > 0) then
        if ($filediff == 0) then
          echo "  compare OK"
        else
          set failure = 1
          echo "  compare FAIL"
        endif
      else
        set failure = 9
        echo "  compare on no output"
      endif
      rm ${base_out}
      rm ${test_out}
    endif
  else
    set failure = 2
    echo "  missing data"
  endif
endif

exit ${failure}

