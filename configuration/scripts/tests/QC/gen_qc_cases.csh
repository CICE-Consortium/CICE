#!/bin/csh -f

set initargv = ( $argv[*] )

# Initialize variables that can be passed as arguments
set helpheader = 0
set dash = "-"
set spval = "UnDeFiNeD"
set machine = ${spval}
set compiler = ${spval}
set grid = ${spval}
set pesx = ${spval}
set queue = ${spval}
set acct = ${spval}
set testid = ${spval}

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
#------------------------------------------------------------

if ( $helpheader > 0) then
cat << EOF1

NAME   
  gen_qc_cases.csh

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

EXAMPLES
    gen_qc_cases.csh -m conrad -e intel -p 8x4

EOF1
endif

if ($helpheader > 1) then
cat << EOF1

      Available --mach and --env combinations are in configuration/scripts/machines and include:
EOF1
      set soptions1 = `ls -1 configuration/scripts/machines | grep Macros | sed 's/Macros.//g' `
      set soptions = `echo $soptions1 | fmt -1 | sort `
      foreach sopt ($soptions)
        echo "             $sopt"
      end
endif

#------------------------------------------------------------
# Read in command line arguments
#------------------------------------------------------------

echo " "
echo "${0}:"

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

#------------------------------------------------------------
# Generate the case directories
#------------------------------------------------------------

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

# Generate the base case
echo "Generating base case"
if ($testid != $spval) then
  set result = `./cice.setup $options -s qc,long --testid qc_base_$testid | grep 'Test case dir\|already exists'`
else
  set result = `./cice.setup $options -s qc,long --testid qc_base | grep 'Test case dir\|already exists'`
endif
set base_dir = `echo "$result" | awk '{print$NF}'`
if ($base_dir == "exists") then
  # Case already exists.  Exit
  echo "$result"
  exit -1
endif

# Generate the BFB case
echo "Generating bfb case"
if ($testid != $spval) then
  set result = `./cice.setup $options -s qc,long --testid qc_bfb_$testid | grep 'Test case dir\|already exists'`
else
  set result = `./cice.setup $options -s qc,long --testid qc_bfb | grep 'Test case dir\|already exists'`
endif
set bfb_dir = `echo "$result" | awk '{print$NF}'`
if ($bfb_dir == "exists") then
  # Case already exists.  Exit
  echo "$result"
  exit -1
endif

# Generate the non-BFB but non-climate-changing case
echo "Generating nonbfb case"
if ($testid != $spval) then
  set result = `./cice.setup $options -s qc_nonbfb,long --testid qc_test_$testid | grep 'Test case dir\|already exists'`
else
  set result = `./cice.setup $options -s qc_nonbfb,long --testid qc_test | grep 'Test case dir\|already exists'`
endif
set nonbfb_dir = `echo "$result" | awk '{print$NF}'`
if ($nonbfb_dir == "exists") then
  # Case already exists.  Exit
  echo "$result"
  exit -1
endif

# Generate the non-BFB and climate changing case
echo "Generating fail case"
if ($testid != $spval) then
  set result = `./cice.setup $options -s alt02,qc,long --testid qc_fail_$testid | grep 'Test case dir\|already exists'`
else
  set result = `./cice.setup $options -s alt02,qc,long --testid qc_fail | grep 'Test case dir\|already exists'`
endif
set fail_dir = `echo "$result" | awk '{print$NF}'`
if ($fail_dir == "exists") then
  # Case already exists.  Exit
  echo "$result"
  exit -1
endif

#------------------------------------------------------------
# Print case directories to file
#------------------------------------------------------------

set QC_DIR = "./qc_logs"
mkdir -p $QC_DIR
echo "$base_dir" > $QC_DIR/qc_dirs.txt
echo "$bfb_dir" >> $QC_DIR/qc_dirs.txt
echo "$nonbfb_dir" >> $QC_DIR/qc_dirs.txt
echo "$fail_dir" >> $QC_DIR/qc_dirs.txt

#------------------------------------------------------------
# cd to each directory, build, and submit
#------------------------------------------------------------

echo ""
echo "Building $base_dir and storing output in $QC_DIR/base_build.txt"
cd $base_dir 
./cice.build >& ../$QC_DIR/base_build.txt
./cice.submit
cd ../ 

echo ""
echo "Building $bfb_dir and storing output in $QC_DIR/bfb_build.txt"
cd $bfb_dir
./cice.build >& ../$QC_DIR/bfb_build.txt
./cice.submit
cd ../

echo ""
echo "Building $nonbfb_dir and storing output in $QC_DIR/nonbfb_build.txt"
cd $nonbfb_dir
./cice.build >& ../$QC_DIR/nonbfb_build.txt
./cice.submit
cd ../

echo ""
echo "Building $fail_dir and storing output in $QC_DIR/fail_build.txt"
cd $fail_dir
./cice.build >& ../$QC_DIR/fail_build.txt
./cice.submit
cd ../

echo ""
echo "Done building cases, and submitting jobs to queue."
echo "When all jobs are completed, run the compare_qc_cases.csh script to"
echo "validate the QC script"
