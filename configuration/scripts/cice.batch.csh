#! /bin/csh -f

if ( $1 != "" ) then
  echo ${0:t} ${1}
else
  echo ${0:t}
endif

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = $1

set ntasks = ${CICE_NTASKS}
set nthrds = ${CICE_NTHRDS}
set maxtpn = ${CICE_MACHINE_TPNODE}
set acct   = ${CICE_MACHINE_ACCT}

@ taskpernode = ${maxtpn} / $nthrds
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}
@ corespernode = ${taskpernodelimit} * ${nthrds}

set ptile = $taskpernode
if ($ptile > ${maxtpn} / 2) @ ptile = ${maxtpn} / 2

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
#BSUB -W ${CICE_RUNLENGTH}
#BSUB -P ${acct}
EOFB

#==========================================

else if (${CICE_MACHINE} =~ cheyenne*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
#PBS -m ae 
#PBS -V
#PBS -q regular
#PBS -N ${CICE_CASENAME}
#PBS -A ${CICE_ACCT}
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${CICE_RUNLENGTH}
EOFB

else if (${CICE_MACHINE} =~ thunder* || ${CICE_MACHINE} =~ gordon* || ${CICE_MACHINE} =~ conrad*) then
cat >> ${jobfile} << EOFB
#PBS -N ${CICE_CASENAME}
#PBS -q debug
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=${CICE_RUNLENGTH}
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
#SBATCH -t ${CICE_RUNLENGTH}
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
#SBATCH -t ${CICE_RUNLENGTH}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=low
EOFB

else if (${CICE_MACHINE} =~ pinto*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${CICE_CASENAME}
#SBATCH -t ${CICE_RUNLENGTH}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=standby
EOFB

else
  echo "${0} ERROR: ${CICE_MACHINE} unknown"
  exit -1
endif

exit 0
