#!/bin/csh 

#echo "hello"

foreach logfile (case*/logs/cice.runlog* testsuite.*/*/logs/cice.runlog*)
  echo "### ${logfile} ###"
  tail -20 $logfile
  echo " "
end

