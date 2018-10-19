#!/bin/csh -f

set initargv = ( $argv[*] ) 

set helpheader = 0
set dash = "-"
set spval = "UnDeFiNeD"
set machcomp = ${spval}
set machine = ${spval}
set compiler = intel
set test = ${spval}
set grid = gx3
set pesx = 4x1
set sets = ""
set testid = ${spval}
set queue = ${spval}
set acct = ${spval}

if ($#argv < 1) then
  set helpheader = 1
endif

set argv = ( $initargv[*] )
# check for -h
while (1)
  if ($#argv < 1) break;
  if ("$argv[1]" =~ "-h" || "$argv[1]" =~ "--help") then
    set helpheader = 2
    if ("$argv[1]" == "-h") then
      set helpheader = 1
    endif
  endif
  shift argv
end

#------------------------------------------------------------
# Help output

if ( $helpheader > 0) then
cat << EOF1

NAME   
  run_qc.csh

SYNOPSIS
    -h || --help

    -m MACH 
        [-e ENV][-p MxN][-g GRID][-s SET1,SET2][--acct ACCT][--queue QUEUE]

DESCRIPTION
    --help, -h : help
    --mach, -m : machine, machine name (required)
    --env,  -e : compiler
    --pes,  -p : tasks x threads [x blocksize_x x blocksize_y [x maxblocks]]
    --acct     : account number for the batch submission
    --grid, -g : grid, grid
    --queue    : queue for the batch submission
    --testid   : test ID, user-defined id for testing
    --set,  -s : case option setting(s), comma separated (default = " ")

EXAMPLES
    run_qc.csh -m conrad -e intel -p 8x4
EOF1

if ($helpheader > 1) then
cat << EOF1

      Available --mach and --env combinations are in configuration/scripts/machines and include:
EOF1
      set soptions1 = `ls -1 configuration/scripts/machines | grep Macros | sed 's/Macros.//g' `
      set soptions = `echo $soptions1 | fmt -1 | sort `
      foreach sopt ($soptions)
        echo "             $sopt"
      end
cat << EOF1

      Available --set options are in configuration/scripts/options and include:
EOF1
      set soptions1 = `ls -1 configuration/scripts/options | grep set_ | sed 's/set_nml.//g' | sed 's/set_env.//g' | sed 's/set_files.//g' `
      set soptions = `echo $soptions1 | fmt -1 | sort -u `
      foreach sopt ($soptions)
        echo "             $sopt"
      end

endif
exit -1
endif

#------------------------------------------------------------
# Read in command line arguments

echo " "
echo "${0}:"

set argv = ( $initargv[*] )
# check for --version
while (1)
  if ($#argv < 1) break;
  if ("$argv[1]" =~ "--version" ) then
    echo "${0}: This is ${ICE_VERSION}"
    exit -1
  endif
  shift argv
end

set argv = ( $initargv[*] )
# read in all options
while (1)
  if ( $#argv < 1 ) break;
  set option = $argv[1];

  shift argv
  if ( $#argv < 1 ) then
    echo "${0}: ERROR1 in $option"
    exit -1
  endif
  if ($argv[1] =~ $dash* ) then
    echo "${0}: ERROR2 in $option"
    exit -1
  endif

  if ("$option" =~ --mach* || "$option" == "-m") then
    set machine = $argv[1]
  else if ("$option" =~ --env* || "$option" == "-e") then
    set compiler = $argv[1]
  else if ("$option" == "--grid" || "$option" == "-g") then
    set grid = $argv[1]
  else if ("$option" == "--queue") then
    set queue = $argv[1]
  else if ("$option" == "--pes" || "$option" == "-p") then
    set pesx = $argv[1]
  else if ("$option" == "--acct") then
    set acct = $argv[1]
  else if ("$option" =~ --set*  || "$option" == "-s") then
    set sets = $argv[1]
  else if ("$option" == "--testid") then
    set testid = $argv[1]
  else
    echo "${0}: ERROR unknown option $option, use -h for help"
    exit -1
  endif

  shift argv

end

if (${machine} == ${spval}) then
  echo "${0}: ERROR in arguments, --mach required"
  exit -1
endif

if ("$compiler" =~ "*,*") then
  echo "${0}: ERROR in arguments, cannot set multiple compilers"
  exit -1
endif

###############################################################################
########## Generate the base data by cloning the master branch ################
###############################################################################

if ( $testid != $spval ) then
  set outdir = "./CICE_master.${testid}"
else
  set outdir = "./CICE_master"
endif

if ( -d $outdir ) then
  echo "$outdir directory already exists!"
  exit -1
endif

# Clone master
git clone --recursive https://github.com/CICE-Consortium/CICE.git $outdir
cd $outdir

# Create a variable that houses all of the arguments
set options = "-m $machine --test smoke"
if ($acct != $spval) then
  set options = "$options --acct $acct"
endif
if ($grid != $spval) then
  set options = "$options --grid $grid"
endif
if ($compiler != $spval) then
  set options = "$options -e $compiler"
endif
if ($pesx != $spval) then
  set options = "$options -p $pesx"
endif
if ($queue != $spval) then
  set options = "$options --queue $queue"
endif

# Build the base test case
if ( $testid == $spval ) then
  ./cice.setup $options -s qc,long --testid qc_master | tee log.txt
  set rc = $?
else
  ./cice.setup $options -s qc,long --testid ${testid}_master | tee log.txt
  set rc = $?
endif
if ( $rc != 0 ) then
  exit $rc
endif

# Get the casename from the log file
set base_casename = `grep 'Test case dir' log.txt | awk '{print $NF}'`
rm log.txt

# Build and run the case
cd $base_casename
./cice.build
set rc = $?
if ( $rc != 0 ) then
  exit $rc
endif
./cice.submit | tee log.txt
set rc = $?
if ( $rc != 0 ) then
  exit $rc
endif
set base_jobid = `grep -oP "\d+" log.txt | sort -n | tail -1`
rm log.txt

# Move back to base directory
cd ../../

###############################################################################
############################ Generate the test data ###########################
###############################################################################

# See if any additional sets were given.  If so, add to options
if ( $sets != "" ) then
  set options = "$options -s ${sets},qc,long"
else
  set options = "$options -s qc,long"
endif

# Build the test case
if ( $testid == $spval ) then
  ./cice.setup $options --testid qc_test | tee log.txt
  set rc = $?
else
  ./cice.setup $options --testid ${testid}_test | tee log.txt
  set rc = $?
endif
if ( $rc != 0 ) then
  exit $rc
endif

# Get the casename from the log file
set test_casename = `grep 'Test case dir' log.txt | awk '{print $NF}'`
rm log.txt

cd $test_casename
./cice.build
set rc = $?
if ( $rc != 0 ) then
  exit $rc
endif
./cice.submit | tee log.txt
set rc = $?
if ( $rc != 0 ) then
  exit $rc
endif
set test_jobid = `grep -oP "\d+" log.txt | sort -n | tail -1`
rm log.txt

cd ..

# Get the ICE_MACHINE_QSTAT variable from the env.* file
set env_file = `ls ${base_casename}/env.*`
source $env_file

# Get the rundir for each test
set basedir = `grep ' ICE_RUNDIR ' $base_casename/cice.settings | awk '{print $NF}'`
set testdir = `grep ' ICE_RUNDIR ' $test_casename/cice.settings | awk '{print $NF}'`
echo ""
echo "---"
echo ""
echo "Waiting for jobs to complete.  If this script is interrupted, run the "
echo "following command (after the jobs finish) to perform the QC test:"
echo ""
echo "`pwd`/configuration/scripts/test/QC/cice.t-test.py ${basedir} ${testdir}"
echo ""

# Wait for both jobs to finish
foreach job ($base_jobid $test_jobid)
  while (1)
    ${ICE_MACHINE_QSTAT} $job >&/dev/null
    if ($? != 0) then
      echo "Job $job completed"
      break
    endif
    echo "Waiting for $job to complete"
    sleep 300  # Sleep for 5 minutes, so as not to overwhelm the queue manager
  end
end

echo ""
echo "Running QC test"
./configuration/scripts/tests/QC/cice.t-test.py $basedir $testdir
