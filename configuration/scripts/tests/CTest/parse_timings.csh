#!/bin/csh -f

# This script parses the timings from the tests in the test suite and writes them to the CTest
# Test.xml file prior to submitting to CDash.

# Loop through each line of the Test.xml file
set CTEST_TAG="`head -n 1 Testing/TAG`"
set testfile="`ls Testing/${CTEST_TAG}/Test.xml`"
set outfile="Testing/${CTEST_TAG}/Test.xml.generated"
set save_time=0
foreach line ("`cat $testfile`")
  if ( "$line" =~ *"Test Status"* ) then
    if ( "$line" =~ *"passed"* ) then
      set save_time=1
    else
      set save_time=0
    endif
  endif
  if ( "$line" =~ *"FullName"* ) then
    if ( "$line" =~ *"_run<"* && $save_time == 1) then
      set save_time=1
      # Grab the case name
      set casename=`echo $line | grep -oP '(?<=<FullName>\.\/).*?(?=</FullName>)' | rev | cut -c 5- | rev`
    else
      set save_time=0
    endif
  endif
  if ( "$line" =~ *"Execution Time"* && $save_time == 1 ) then
    # Find the case runlog
    set runlog=`ls ./${casename}.*/logs/*runlog*`
    foreach line1 ("`cat $runlog`")
      if ( "$line1" =~ *"Timer   2:"*) then
        set runtime=`echo $line1 | grep -oP "\d+\.(\d+)?" | sort -n | tail -1`
      endif
    end
    set local_runtime=`echo $line | grep -oP "\d+\.(\d+)?" | sort -n | tail -1`
    # Grab the leading whitespace
    # Replace the timing in Test.xml with the timing from the runlog file
    set line=`echo "$line" | sed "s/^\(.*<Value>\)${local_runtime}\(<\/Value>\)/\1${runtime}<\/Value>/"`
    set save_time=0
  endif
  echo "$line" >> $outfile
end

mv $testfile ${testfile}.bak
mv $outfile $testfile
