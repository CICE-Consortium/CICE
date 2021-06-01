#!/bin/csh -f

if ( $1 != "" ) then
  echo "running cice.batch.csh (creating ${1})"
else
  echo "running cice.batch.csh"
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
if (${taskpernode} == 0) set taskpernode = 1
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}
@ corespernode = ${taskpernodelimit} * ${nthrds}

set ptile = $taskpernode
if ($ptile > ${maxtpn} / 2) @ ptile = ${maxtpn} / 2

set runlength = ${ICE_RUNLENGTH}
if ($?ICE_MACHINE_MAXRUNLENGTH) then
  if (${runlength} > ${ICE_MACHINE_MAXRUNLENGTH}) then
    set runlength = ${ICE_MACHINE_MAXRUNLENGTH}
  endif
endif

set queue = "${ICE_QUEUE}"
set batchtime = "00:15:00"
if (${runlength} == 0) set batchtime = "00:29:00"
if (${runlength} == 1) set batchtime = "00:59:00"
if (${runlength} == 2) set batchtime = "2:00:00"
if (${runlength} == 3) set batchtime = "3:00:00"
if (${runlength} == 4) set batchtime = "4:00:00"
if (${runlength} == 5) set batchtime = "5:00:00"
if (${runlength} == 6) set batchtime = "6:00:00"
if (${runlength} == 7) set batchtime = "7:00:00"
if (${runlength} >= 8) set batchtime = "8:00:00"

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

else if (${ICE_MACHINE} =~ izumi*) then
if (${runlength} > 2) set queue = "medium"
cat >> ${jobfile} << EOFB
#PBS -j oe 
###PBS -m ae 
#PBS -V
#PBS -q ${queue}
#PBS -N ${ICE_CASENAME}
#PBS -l nodes=${nnodes}:ppn=${taskpernode}
#PBS -l walltime=${batchtime}
EOFB

else if (${ICE_MACHINE} =~ gordon* || ${ICE_MACHINE} =~ conrad*  || ${ICE_MACHINE} =~ gaffney* || ${ICE_MACHINE} =~ koehr* || ${ICE_MACHINE} =~ mustang) then
cat >> ${jobfile} << EOFB
#PBS -N ${shortcase}
#PBS -q ${queue}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=${batchtime}
#PBS -j oe
#PBS -W umask=022
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ onyx*) then
# special for onyx with 44 cores per node and constraint on mpiprocs
set tpn1 = ${taskpernode}
if (${taskpernode} < 44) set tpn1 = 22
if (${taskpernode} < 22) set tpn1 = 11
if (${taskpernode} < 11) set tpn1 =  4
if (${taskpernode} <  4) set tpn1 =  2
if (${taskpernode} <  2) set tpn1 =  1
@ nn1 = ${ntasks} / ${tpn1}
if (${nn1} * ${tpn1} < ${ntasks}) @ nn1 = $nn1 + 1
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -q ${queue}
#PBS -A ${acct}
#PBS -l select=${nn1}:ncpus=${maxtpn}:mpiprocs=${tpn1}
#PBS -l walltime=${batchtime}
#PBS -j oe
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ cori*) then
@ nthrds2 = ${nthrds} * 2
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -A ${acct}
#SBATCH --qos ${queue}
#SBATCH --time ${batchtime}
#SBATCH --nodes ${nnodes}
#SBATCH --ntasks ${ntasks}
#SBATCH --cpus-per-task ${nthrds2}
#SBATCH --constraint haswell
###SBATCH -e filename
###SBATCH -o filename
###SBATCH --mail-type FAIL
###SBATCH --mail-user username@domain.com
EOFB

else if (${ICE_MACHINE} =~ compy*) then
if (${runlength} <= 2) set queue = "short"
cat >> ${jobfile} <<EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -A ${acct}
#SBATCH --qos ${queue}
#SBATCH --time ${batchtime}
#SBATCH --nodes ${nnodes}
#SBATCH --ntasks ${ntasks}
#SBATCH --cpus-per-task ${nthrds}
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

else if (${ICE_MACHINE} =~ millikan*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -t ${batchtime}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=amelie.bouchat@canada.ca
#SBATCH --qos=standby
EOFB

else if (${ICE_MACHINE} =~ daley* || ${ICE_MACHINE} =~ banting* ) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -j oe
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${batchtime}
EOFB

else if (${ICE_MACHINE} =~ freya* ) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -j oe
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${batchtime}
EOFB

else if (${ICE_MACHINE} =~ hera*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH --partition=hera
#SBATCH --qos=${queue}
#SBATCH -A ${acct}
#SBATCH --time=${batchtime}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=${taskpernodelimit}
#SBATCH --cpus-per-task=${nthrds}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
##SBATCH --mail-type FAIL
##SBATCH --mail-user=xxx@noaa.gov
EOFB

else if (${ICE_MACHINE} =~ orion*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH --partition=orion
#SBATCH --qos=${queue}
#SBATCH -A ${acct}
#SBATCH --time=${batchtime}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=${taskpernodelimit}
#SBATCH --cpus-per-task=${nthrds}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
##SBATCH --mail-type FAIL
##SBATCH --mail-user=xxx@noaa.gov
EOFB

else if (${ICE_MACHINE} =~ phase3*) then
if ( ${nnodes} > 15) then
  setenv p3tile 16
  setenv mem `expr 100 \* 1024 / $nnodes`
else
  setenv p3tile ${nnodes}
  setenv mem 8192
endif
echo mem = ${mem} nnodes and p3tiles ${nnodes} ${p3tile} p3tile must be le nnodes
cat >> ${jobfile} << EOFB
#BSUB -J ${ICE_CASENAME}
#BSUB -q "dev_shared"
#BSUB -P RTO-T2O
#BSUB -W `echo ${batchtime} | cut -f1-2 -d:`
#BSUB -n ${nnodes}
#BSUB -R "affinity[core]"
#BSUB -R "span[ptile=${p3tile}]"
#BSUB -R "rusage[mem=${mem}]"
#BSUB -o /u/Robert.Grumbine/${ICE_CASENAME}.out.%J
#BSUB -e /u/Robert.Grumbine/${ICE_CASENAME}.err.%J
EOFB

else if (${ICE_MACHINE} =~ high_Sierra*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ testmachine*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ travisCI*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ conda*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB


else
  echo "${0} ERROR: ${ICE_MACHINE} unknown"
  exit -1
endif

exit 0
