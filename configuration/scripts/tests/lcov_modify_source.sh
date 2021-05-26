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
  cat $ofile | while read -r line || [[ -n $oline ]]; do

    if [[ $contblock == 1 ]]; then
        if [[ $line =~ ^.*"&".*$ ]]; then
#            echo ${line} ${LCOV_EXCL}
            echo ${line} ${LCOV_EXCL} >> ${nfile}
        else
            contblock=0
#            echo ${line}
            echo ${line} >> ${nfile}
        fi
    else
        echo ${line} >> ${nfile}
        if [[ $line =~ ^\s*.*"&".*$ && ! $line =~ ^\s*( if ).*$ ]]; then
            contblock=1
#            echo ${line}
        fi
    fi

  done

done
