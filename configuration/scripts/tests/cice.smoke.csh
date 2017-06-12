#! /bin/csh -f

echo ${0}

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.smoke
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

# Create test script that runs cice.build, cice.run, and validates
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
#SBATCH -J ${CICE_CASENAME}
#SBATCH -t 0:10:00
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
#SBATCH --mail-type FAIL
#SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=low
EOFB

endif

cat >> ${jobfile} << EOF2
cd ${CICE_CASEDIR}
source ./cice.settings || exit 2
source ./env.\${CICE_MACHINE} || exit 2

# Compile the CICE code
./cice.build
set rc_build = \$?

# Run the CICE model
./cice.run
set rc_run = \$?

EOF2

if ($1 != "") then
cat >> ${jobfile} << EOF3
  # Get the final output filename
  foreach file (\${CICE_RUNDIR}/restart/*)
    set test_data = \$file
  end
  
  set baseline_data = $1/\$test_data:t

  if ( { cmp -s \$test_data \$baseline_data } ) then
    set rc_valid = 0
  else
    set rc_valid = 1
  endif
EOF3
endif

cat >> ${jobfile} << EOF4

if (\$rc_build == 0) then
  echo "Build:      [0;92mPASS[0;0m"
else
  echo "Build:      [0;41mFAIL[0;0m"
endif

if (\$rc_run == 0) then
  echo "Run:        [0;92mPASS[0;0m"
else
  echo "Run:        [0;41mFAIL[0;0m"
endif
EOF4

if ($1 != "") then
cat >> ${jobfile} << EOF5

if (\$rc_valid == 0) then
  echo "Validation: [0;92mPASS[0;0m"
else
  echo "Validation: [0;41mFAIL[0;0m"
endif
EOF5
endif

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.test

EOFS

chmod +x ${subfile}
