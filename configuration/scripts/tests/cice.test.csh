#! /bin/csh -f

echo ${0}

# Define this as a test case in the cice.settings file
cat >> ./cice.settings << EOF

# Add variable defining this case_dir as a test
setenv CICE_TEST true
EOF

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.test
set subfile = cice.submit

set ntasks = ${CICE_NTASKS}
set nthrds = ${CICE_NTHRDS}
set maxtpn = ${CICE_MACHINE_TPNODE}
set acct   = ${CICE_MACHINE_ACCT}

@ taskpernode = ${maxtpn} / $nthrds
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}

set ptile = $taskpernode
if ($ptile > ${maxtpn} / 2) @ ptile = ${maxtpn} / 2

#==========================================

# Create test script that runs cice.run, and validates
#==========================================

cat >! ${jobfile} << EOF0
#!/bin/csh -f 
EOF0

#==========================================

if (${CICE_MACHINE} =~ yellowstone*) then
cat >> ${jobfile} << EOFB
#BSUB -n ${ntasks}
#BSUB -R "span[ptile=${ptile}]"
#BSUB -q caldera
#BSUB -N
###BSUB -x
#BSUB -a poe
#BSUB -o poe.stdout.%J
#BSUB -e poe.stderr.%J
#BSUB -J ${CICE_CASENAME}
#BSUB -W 0:10
#BSUB -P ${acct}
EOFB

else if (${CICE_MACHINE} =~ thunder* || ${CICE_MACHINE} =~ gordon* || ${CICE_MACHINE} =~ conrad*) then
cat >> ${jobfile} << EOFB
#PBS -N ${CICE_CASENAME}
#PBS -q debug
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=0:10:00
#PBS -j oe
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${CICE_MACHINE} =~ cori*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${CICE_CASENAME}
#SBATCH -p debug
###SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -t 0:10:00
#SBATCH -L SCRATCH
#SBATCH -C haswell
###SBATCH -e filename
###SBATCH -o filename
###SBATCH --mail-type FAIL
###SBATCH --mail-user username@domain.com
EOFB

else if (${CICE_MACHINE} =~ wolf*) then
cat >> ${jobfile} << EOFB
#MSUB -A climateacme
#MSUB -l walltime=01:00:00
#MSUB -l nodes=1:ppn=4
#MSUB -N cice
EOFB
#cat >> ${jobfile} << EOFB
##SBATCH -J ${CICE_CASENAME}
####SBATCH -p tossdev
##SBATCH -t 0:45:00
##SBATCH -A ${acct}
##SBATCH -N ${nnodes}
####SBATCH -L SCRATCH
####SBATCH -C haswell
##SBATCH -e slurm%j.err
##SBATCH -o slurm%j.out
####SBATCH --mail-type FAIL
####SBATCH --mail-user username@domain.com
#EOFB

endif

cat >> ${jobfile} << EOF2
cd ${CICE_CASEDIR}
source ./cice.settings || exit 2
source ./env.\${CICE_MACHINE} || exit 2

# Check to see if executable exists in CICE_RUNDIR
if ( ! -f \${CICE_RUNDIR}/cice ) then
  echo "cice executable does not exist in \${CICE_RUNDIR}.  "
  echo "Please run cice.build before this test."
  exit 99
endif

# Run the CICE model
./cice.run

EOF2

if ($1 != "") then
cat >> ${jobfile} << EOF3
  # Get the final output filename
  foreach file (\${CICE_RUNDIR}/restart/*)
    set test_data = \$file
  end
  
  set baseline_data = $1/\$test_data:t

  echo "Performing binary comparison between files:"
  echo "baseline: \$baseline_data"
  echo "test:     \$test_data"
  if ( { cmp -s \$test_data \$baseline_data } ) then
    echo "PASS \${CICE_CASENAME} test" >> ${CICE_CASEDIR}/test_output
  else
    echo "FAIL \${CICE_CASENAME} test" >> ${CICE_CASEDIR}/test_output
  endif
EOF3
endif

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
