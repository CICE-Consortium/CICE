#! /bin/csh -f

if ( $1 != "" ) then
  echo ${0:t} ${1}
else
  echo ${0:t}
endif

#source ./cice.settings
#source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = $1

set ntasks = ${ICE_NTASKS}
set nthrds = ${ICE_NTHRDS}
set maxtpn = ${ICE_MACHINE_TPNODE}
set acct   = ${ICE_ACCOUNT}

@ ncores = ${ntasks} * ${nthrds}
@ taskpernode = ${maxtpn} / $nthrds
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}
@ corespernode = ${taskpernodelimit} * ${nthrds}

set ptile = $taskpernode
if ($ptile > ${maxtpn} / 2) @ ptile = ${maxtpn} / 2

set queue = "${ICE_QUEUE}"
set batchtime = "00:15:00"
if (${ICE_RUNLENGTH} >  0) set batchtime = "00:29:00"
if (${ICE_RUNLENGTH} == 1) set batchtime = "00:59:00"
if (${ICE_RUNLENGTH} == 2) set batchtime = "2:00:00"
if (${ICE_RUNLENGTH} == 3) set batchtime = "3:00:00"
if (${ICE_RUNLENGTH} == 4) set batchtime = "4:00:00"
if (${ICE_RUNLENGTH} == 5) set batchtime = "5:00:00"
if (${ICE_RUNLENGTH} == 6) set batchtime = "6:00:00"
if (${ICE_RUNLENGTH} == 7) set batchtime = "7:00:00"
if (${ICE_RUNLENGTH} >= 8) set batchtime = "8:00:00"

set shortcase = `echo ${ICE_CASENAME} | cut -c1-15`

#==========================================

cat >! ${jobfile} << EOF0
#!/bin/csh -f 
EOF0

#==========================================

if (${ICE_MACHINE} =~ cheyenne*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
###PBS -m ae 
#PBS -V
#PBS -q ${queue}
#PBS -N ${ICE_CASENAME}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${batchtime}
EOFB

else if (${ICE_MACHINE} =~ hobart*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
###PBS -m ae 
#PBS -V
#PBS -q short
#PBS -N ${ICE_CASENAME}
#PBS -l nodes=1:ppn=24
EOFB

else if (${ICE_MACHINE} =~ thunder* || ${ICE_MACHINE} =~ gordon* || ${ICE_MACHINE} =~ conrad*) then
cat >> ${jobfile} << EOFB
#PBS -N ${shortcase}
#PBS -q ${queue}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=${batchtime}
#PBS -j oe
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ onyx*) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -q ${queue}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=${batchtime}
#PBS -j oe
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ cori*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -p ${queue}
###SBATCH -A ${acct}
#SBATCH -n ${ncores}
#SBATCH -t ${batchtime}
#SBATCH -L SCRATCH
#SBATCH -C haswell
###SBATCH -e filename
###SBATCH -o filename
###SBATCH --mail-type FAIL
###SBATCH --mail-user username@domain.com
EOFB

else if (${ICE_MACHINE} =~ badger*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -t ${batchtime}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=standby
EOFB

else if (${ICE_MACHINE} =~ fram*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -t ${batchtime}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=armnjfl@ec.gc.ca
#SBATCH --qos=standby
EOFB

else if (${ICE_MACHINE} =~ cesium*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -t ${batchtime}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=philippe.blain@canada.ca
#SBATCH --qos=standby
EOFB

else if (${ICE_MACHINE} =~ testmachine*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ travisCI*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else
  echo "${0} ERROR: ${ICE_MACHINE} unknown"
  exit -1
endif

exit 0
