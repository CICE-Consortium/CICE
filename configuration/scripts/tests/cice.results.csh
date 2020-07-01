
#--- cice.results.csh --- 

cat ./results.log | grep -iv "#machinfo" | grep -iv "#envinfo"
set pends     = `cat ./results.log | grep PEND | wc -l`
set misses    = `cat ./results.log | grep MISS | wc -l`
set failures  = `cat ./results.log | grep FAIL | wc -l`
set failbuild = `cat ./results.log | grep FAIL | grep " build " | wc -l`
set failrun   = `cat ./results.log | grep FAIL | grep " run " | wc -l`
set failtest  = `cat ./results.log | grep FAIL | grep " test " | wc -l`
set failcomp  = `cat ./results.log | grep FAIL | grep " compare " | wc -l`
set failbfbc  = `cat ./results.log | grep FAIL | grep " bfbcomp " | wc -l`
set failgen   = `cat ./results.log | grep FAIL | grep " generate " | wc -l`
set success   = `cat ./results.log | grep 'PASS\|COPY' | wc -l`
set comments  = `cat ./results.log | grep "#" | wc -l`
set alltotal  = `cat ./results.log | wc -l`
@ total = $alltotal - $comments
@ chkcnt = $pends + $misses + $failures + $success

echo "#------- " >> results.log
echo " " >> results.log
echo "#totl = $total total" >> results.log
echo "#chkd = $chkcnt checked" >> results.log
echo "#pass = $success" >> results.log
echo "#pend = $pends" >> results.log
echo "#miss = $misses" >> results.log
echo "#fail = $failures" >> results.log
echo "    #failbuild = $failbuild" >> results.log
echo "    #failrun   = $failrun" >> results.log
echo "    #failtest  = $failtest" >> results.log
echo "    #failcomp  = $failcomp" >> results.log
echo "    #failbfbc  = $failbfbc" >> results.log
echo "    #failgen   = $failgen" >> results.log

echo ""
echo "Descriptors:"
echo " PASS - successful completion"
echo " COPY - previously compiled code was copied for new test"
echo " MISS - comparison data is missing"
echo " PEND - status is undetermined; test may still be queued, running, or timed out"
echo " FAIL - test failed"
echo ""
echo "$chkcnt measured results of $total total results"
echo "$success of $chkcnt tests PASSED"
echo "$pends of $chkcnt tests PENDING"
echo "$misses of $chkcnt tests MISSING data"
echo "$failures of $chkcnt tests FAILED"
#echo "  $failbuild of $failures FAILED build"
#echo "  $failrun of $failures FAILED run"
#echo "  $failtest of $failures FAILED test"
#echo "  $failcomp of $failures FAILED compare"
#echo "  $failbfbc of $failures FAILED bfbcomp"
#echo "  $failgen of $failures FAILED generate"
exit $failures
