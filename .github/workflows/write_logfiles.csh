#!/bin/csh 

#echo "hello"

foreach logfile (case*/logs/cice.runlog* testsuite.*/*/logs/cice.runlog*)
  echo "### ${logfile} ###"
  tail -200 $logfile
  echo " "
end

