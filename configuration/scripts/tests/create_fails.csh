#!/bin/csh

echo " "
set tmpfile = create_fails.tmp
set outfile = fails.ts

./results.csh >& /dev/null
cat results.log | grep ' run\| test' | grep -v "#" | grep -v PASS | cut -f 2 -d " " | sort -u >! $tmpfile

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
echo "$0 done"
echo "Failed tests can be rerun with the test suite file...... $outfile"
echo "To run a new test suite, copy $outfile to the top directory and do something like"
echo "  ./cice.setup --suite $outfile ..."
