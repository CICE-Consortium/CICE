#!/bin/csh

echo " "
set tmpfile = create_fails.tmp
set outfile = fails.ts
set runfile = rerun.csh

set delim = `pwd | rev | cut -d . -f 1 | rev`

./results.csh >& /dev/null

#fails, both "run" and "test" that failed
# treat decomp special
cat results.log | grep ' run\| test' | grep -v "#" | grep -v PASS | cut -f 2 -d " " | grep -v _decomp_ | sort -u >! $tmpfile
cat results.log | grep ' run\| test' | grep -v "#" | grep -v PASS | cut -f 2 -d " " | grep _decomp_  | rev | cut -d _ -f 2- | rev | sort -u >> $tmpfile

echo "# Test  Grid  PEs  Sets" >! $outfile
foreach line ( "`cat $tmpfile`" )
  #echo $line
  set test = `echo $line | cut -d "_" -f 3`
  set grid = `echo $line | cut -d "_" -f 4`
  set pes  = `echo $line | cut -d "_" -f 5`
  set opts = `echo $line | cut -d "_" -f 6- | sed 's/_/,/g'`
  echo "$test  $grid  $pes  $opts" >> $outfile
end
rm $tmpfile

#rerun, only "run" that failed
# treat decomp special
cat results.log | grep ' run' | grep -v "#" | grep -v PASS | cut -f 2 -d " " | grep -v _decomp_ | sort -u >! $tmpfile
cat results.log | grep ' run' | grep -v "#" | grep -v PASS | cut -f 2 -d " " | grep _decomp_  | rev | cut -d _ -f 2- | rev | sort -u >> $tmpfile

echo "#/bin/csh" >! $runfile
foreach line ( "`cat $tmpfile`" )
  #echo $line
  echo "cd ${line}.${delim}; ./*.submit; cd ../; sleep 5" >> $runfile
end
chmod +x $runfile
rm $tmpfile

echo "$0 done"
echo " "
echo "Failed runs can be resubmitted by running $runfile"
echo " "
echo "Failed tests can be rerun with the test suite file...... $outfile"
echo "To run a new test suite, copy $outfile to the top directory and do something like"
echo "  ./cice.setup --suite $outfile ..."
