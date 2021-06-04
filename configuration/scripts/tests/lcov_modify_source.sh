#!/bin/bash

filelist=`find cicecore icepack -type f -name "*.F90"`
LCOV_EXCL="  ! LCOV_EXCL_LINE"

#echo $filelist

for file in $filelist; do

  echo $file
  ofile=${file}.orig
  nfile=${file}

  mv ${file} ${file}.orig

  # line by line making sure each line has a trailing newline (-n)
  # preserve whitespace (IFS)
  # and include backslashes (-r)
  IFS=''
  contblock=0
  cat $ofile | while read -r line || [[ -n $line ]]; do

    if [[ $contblock == 1 ]]; then
        # in a continuation block
        if [[ $line =~ ^.*"&".*$ ]]; then
            # found another continuation line, add exclude string and write out line
            echo ${line} ${LCOV_EXCL} >> ${nfile}
        else
            # continuation block ends, write out line
            contblock=0
            echo ${line} >> ${nfile}
        fi
    else
        # not in a continuation block, write out line
        echo ${line} >> ${nfile}
        if [[ $line =~ ^\s*.*"&".*$ && ! $line =~ ^\s*( if ).*$ ]]; then
            # new continuation block found
            contblock=1
        fi
    fi

  done

done
