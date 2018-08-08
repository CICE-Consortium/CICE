# Generate the base case
echo "Generating base case"
export base_dir=`./cice.setup -m conrad --test smoke -s qc,long --testid qc_base --acct ARLAP96070PET | grep 'Test case dir' | awk '{print$NF}'`
echo "base dir = $base_dir"

# Generate the BFB case
echo "Generating bfb case"
export bfb_dir=`./cice.setup -m conrad --test smoke -s qc,long --testid qc_bfb --acct ARLAP96070PET | grep 'Test case dir' | awk '{print$NF}'`

# Generate the non-BFB but non-climate-changing run
echo "Generating nonbfb case"
export nonbfb_dir=`./cice.setup -m conrad --test smoke -s qc_nonbfb,long --testid qc_test --acct ARLAP96070PET | grep 'Test case dir' | awk '{print$NF}'`

# Generate the non-BFB and climate changing test
echo "Generating fail case"
export fail_dir=`./cice.setup -m conrad --test smoke -s alt02,qc,long --testid qc_fail --acct ARLAP96070PET | grep 'Test case dir' | awk '{print$NF}'`

echo "$base_dir" > qc_dirs.txt
echo "$bfb_dir" >> qc_dirs.txt
echo "$nonbfb_dir" >> qc_dirs.txt
echo "$fail_dir" >> qc_dirs.txt

# cd to each directory, build, and submit
cd $base_dir && ./cice.build && ./cice.submit && cd ../
cd $bfb_dir && ./cice.build && ./cice.submit && cd ../
cd $nonbfb_dir && ./cice.build && ./cice.submit && cd ../
cd $fail_dir && ./cice.build && ./cice.submit && cd ../
