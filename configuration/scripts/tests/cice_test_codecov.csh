#!/bin/csh 

# This was a script on gordon
# Should run this interactively with cut and paste due to manual work needed to generate 
#   test_cice_icepack repository

#PBS -N cice_test
#PBS -q standard
#PBS -A NRLSS03755018
#PBS -l application=Regional-Arctic-System-Model
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M anthony.p.craig@gmail.com
#PBS -m be

module load costinit git

set scrdir = "~"
set testdir = "~/cice_testing"
set date0 = `date -u "+%y%m%d"`
set date = ${date0}cc

mkdir -p ${testdir}

cd ${testdir}

# Check out current cice master
git clone https://github.com/cice-consortium/cice cice.master.${date} --recursive
cd cice.master.${date}
set hash = `git show --oneline -s | cut -f 1 -d " "`
cd ../

# Check out test_cice_icepack and update from cice master
git clone https://github.com/apcraig/test_cice_icepack test_cice_icepack.${date}
cd test_cice_icepack.${date}
git rm -r *
cp -p -r ../cice.master.${date}/* .

# Manually copy missed files if needed (should be just dot files, do not copy in .gitmodules)
diff -r ../cice.master.${date} . --exclude .git
cp as needed

# Clean up icepack .git stuff and Commit
rm -r -f icepack/.git*
git add .
git commit -m "update test_cice_icepack master to ${hash}"

# Push test_cice_icepack 
git push origin master

# Run test suite
./cice.setup --suite first_suite,base_suite,travis_suite,decomp_suite,reprosum_suite,quick_suite -m gordon -e gnu --testid T${date} --codecov --queue standard

