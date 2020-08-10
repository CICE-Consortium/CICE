#!/bin/csh 

# This was a script on gordon
# This script should only be run on hardware with the gnu compiler
# This should be run interactively because git push will require login information

#PBS -N cice_test
#PBS -q standard
#PBS -A NRLSS03755018
#PBS -l application=Regional-Arctic-System-Model
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M anthony.p.craig@gmail.com
#PBS -m be

#set verbose
module load costinit git

set scrdir = "~"
set testdir = "~/cice_testing"
set date0 = `date -u "+%y%m%d"`
set date = ${date0}cc

mkdir -p ${testdir}

cd ${testdir}

# Check out current cice master
echo " "
echo "*** checkout current cice master ***"
git clone --depth=1 https://github.com/cice-consortium/cice cice.master.${date} --recursive
cd cice.master.${date}
set hash = `git rev-parse --short HEAD `
cd ../

# Check out test_cice_icepack, remove all code and copy from cice master
# Need to be careful about dot files, particularly .git* files
# This copies in all files via standard file expansion (except dot files at root)
# This also copies in all dot file at the root that do not start with .g (ie. .git*)
echo " "
echo "*** checkout current test_cice_master ***"
git clone --depth=1 https://github.com/apcraig/test_cice_icepack test_cice_icepack.${date}
cd test_cice_icepack.${date}
echo " "
echo "*** remove current files and copy in cice master files ***"
set verbose
git rm -r *  >& /dev/null
cp -p -r ../cice.master.${date}/* .
cp -p    ../cice.master.${date}/.[a-f,h-z]* .

# Clean up icepack .git stuff and commit
rm -r -f icepack/.git*
git add .
unset verbose
echo " "
echo "*** git status of changes ***"
git status
echo " "
echo "*** commit changes ***"
git commit -m "update test_cice_icepack master to ${hash}"

# Push test_cice_icepack 
echo " "
echo "*** push changes to test_cice_icepack ***"
git push origin master

# Run test suite
echo " "
echo "*** run test suite ***"
./cice.setup --suite first_suite,base_suite,travis_suite,decomp_suite,reprosum_suite,io_suite,quick_suite -m gordon -e gnu --testid T${date} --coverage --queue standard

# The test suite will wait until all jobs are complete then run report_codecov.csh or report_lcov.csh
# If that fails, you can run report_codecov.csh or report_lcov.csh manually after all jobs are done
