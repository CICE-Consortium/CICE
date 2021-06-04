#!/bin/csh -f

# Compare prognostic output in two log files
#-----------------------------------------------------------

# usage: comparelog.csh base_file test_file [notcicefile]
#     does diff of two files
#     optional 3rd argument indicates the file is not a cice file so diff entire thing
#
# Return Codes (depends on quality of error checking)
#  0 = pass
#  1 = fail
#  2 = missing data
#  9 = error

set filearg = 0
set cicefile = 0
set notcicefile = "notcicefile"
if ( $#argv == 2 ) then
  set cicefile = 1
  set filearg = 1
  set base_data = $argv[1]
  set test_data = $argv[2]
  if ("$argv[1]" == "${notcicefile}") set filearg = 0
  if ("$argv[2]" == "${notcicefile}") set filearg = 0
else if ( $#argv == 3 ) then
  set cicefile = 0
  set filearg = 1
  set base_data = $argv[1]
  set test_data = $argv[2]
  if ("$argv[3]" != "${notcicefile}") set filearg = 0
endif

if (${filearg} == 0) then
  echo "Error in ${0}"
  echo "Usage: ${0} <base_file> <test_file> [notcicefile]"
  echo "   does diff of two files"
  exit 9
endif

set failure = 0
set base_out = "comparelog_base_out_file.log"
set test_out = "comparelog_test_out_file.log"

if (${filearg} == 1) then
  echo "base_data: $base_data"
  echo "test_data: $test_data"
  if ( -f ${base_data} && -f ${test_data}) then
    if (${base_data} == ${test_data}) then
      set failure = 9
      echo "  input data are same"
    else

      touch ${base_out}
      touch ${test_out}

      if (${cicefile} == 1) then
        cat ${base_data} | grep -A 99999999 "total ice area  (km^2)" | grep -e istep1: -e = | grep -iv "min, max, sum" >&! ${base_out}
        cat ${test_data} | grep -A 99999999 "total ice area  (km^2)" | grep -e istep1: -e = | grep -iv "min, max, sum" >&! ${test_out}
      else
        cp -f ${base_data} ${base_out}
        cp -f ${test_data} ${test_out}
      endif

      set basenum = `cat ${base_out} | wc -l`
      set testnum = `cat ${test_out} | wc -l`
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

