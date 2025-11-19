#!/bin/csh -f

if ( $1 != "" ) then
  echo "running cice.batch.csh (creating ${1})"
else
  echo "running cice.batch.csh"
endif

#source ./cice.settings
#source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = $1

source ${ICE_SCRIPTS}/setup_machparams.csh

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

else if (${ICE_MACHINE} =~ derecho*) then
set memstr = ""
if (${ncores} <= 8 && ${runlength} <= 1 && ${batchmem} <= 20) then
  set queue = "develop"
  set memstr = ":mem=${batchmem}GB"
endif
cat >> ${jobfile} << EOFB
#PBS -q ${queue}
#PBS -l job_priority=regular
#PBS -N ${ICE_CASENAME}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}${memstr}
#PBS -l walltime=${batchtime}
#PBS -j oe
#PBS -W umask=022
#PBS -o ${ICE_CASEDIR}
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ gadi*) then
if (${queue} =~ *sr) then #sapphire rapids
  @ memuse = ( $ncores * 481 / 100 )
  set corespernode = 52
else if (${queue} =~ *bw) then #broadwell
  @ memuse = ( $ncores * 457 / 100 )
  set corespernode = 28
else if (${queue} =~ *sl) then 
  @ memuse = ( $ncores * 6 )
  set corespernode = 32
else #normal queues
  @ memuse = ( $ncores * 395 / 100 )
  set corespernode = 48
endif
if (${ncores} > ${corespernode}) then
  # need to use a whole number of nodes
  @ ncores = ( ( $ncores + $corespernode - 1 ) / $corespernode ) * $corespernode
endif
if (${runlength} <= 0) then
  set batchtime = "00:30:00"
endif
cat >> ${jobfile} << EOFB
#PBS -q ${queue}
#PBS -P ${ICE_MACHINE_PROJ}
#PBS -N ${ICE_CASENAME}
#PBS -l storage=gdata/${ICE_MACHINE_PROJ}+scratch/${ICE_MACHINE_PROJ}+gdata/ik11
#PBS -l ncpus=${ncores}
#PBS -l mem=${memuse}gb
#PBS -l walltime=${batchtime}
#PBS -j oe 
#PBS -W umask=003
#PBS -o ${ICE_CASEDIR}
source /etc/profile.d/modules.csh
module use `echo ${MODULEPATH} | sed 's/:/ /g'` #copy the users modules
EOFB

else if (${ICE_MACHINE} =~ gust*) then
cat >> ${jobfile} << EOFB
#PBS -q ${queue}
#PBS -l job_priority=regular
#PBS -N ${ICE_CASENAME}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${batchtime}
#PBS -j oe 
#PBS -W umask=022
#PBS -o ${ICE_CASEDIR}
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

else if (${ICE_MACHINE} =~ gaffney* || ${ICE_MACHINE} =~ koehr* || ${ICE_MACHINE} =~ mustang* || ${ICE_MACHINE} =~ carpenter*) then
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

else if (${ICE_MACHINE} =~ blueback*) then
if (${runlength} > 0) set queue = "standard"
cat >> ${jobfile} << EOFB
#SBATCH --job-name=${ICE_CASENAME}
#SBATCH --account=${acct}
#SBATCH --qos=${queue}
#SBATCH --time=${batchtime}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks=${ntasks}
#SBATCH --ntasks-per-node=${taskpernode}
#SBATCH --constraint standard
###SBATCH -e filename
###SBATCH -o filename
###SBATCH --mail-type FAIL
###SBATCH --mail-user username@domain.com
EOFB

else if (${ICE_MACHINE} =~ narwhal*) then
if (${runlength} <= 0) then
  set batchtime = "00:29:59"
else
  set queue = "standard"
endif
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

else if (${ICE_MACHINE} =~ nrlssc*) then
# nrlssc queue system has nodes with different task per node
if (${taskpernode} <= 12) set tpnstr = 'twelve'
if (${taskpernode} == 20) set tpnstr = 'twenty'
if (${taskpernode} >= 24) set tpnstr = 'twentyfour'
#if (${taskpernode} == 28) set tpnstr = 'twentyeight'

cat >> ${jobfile} <<EOFB
#PBS -N ${shortcase}
#PBS -q ${queue}
#PBS -l nodes=${nnodes}:ppn=${taskpernode}:${tpnstr}
#PBS -l walltime=${batchtime}
#PBS -j oe
#PBS -W umask=022
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

else if (${ICE_MACHINE} =~ perlmutter*) then
@ nthrds2 = ${nthrds} * 2
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -A ${acct}
#SBATCH --qos=${queue}
#SBATCH --time=${batchtime}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks=${ntasks}
#SBATCH --cpus-per-task=${nthrds2}
#SBATCH --constraint cpu
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

else if (${ICE_MACHINE} =~ chicoma*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -t ${batchtime}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=eclare@lanl.gov
##SBATCH --qos=debug
#SBATCH --qos=standard
##SBATCH --qos=standby
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
#PBS -W umask=022
EOFB

else if (${ICE_MACHINE} =~ robert* || ${ICE_MACHINE} =~ underhill* || ${ICE_MACHINE} =~ ppp6* || ${ICE_MACHINE} =~ ppp5*) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -j oe
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}:mem=20gb
#PBS -l walltime=${batchtime}
#PBS -W umask=022
EOFB

else if (${ICE_MACHINE} =~ ppp3*) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -j oe
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}:mem=20gb:res_tmpfs=1000:res_image=eccc/eccc_all_ppp_ubuntu-18.04-amd64_latest
#PBS -l walltime=${batchtime}
#PBS -W umask=022
#PBS -q development
#PBS -o ${ICE_CASEDIR}
#PBS -S /bin/csh
EOFB

else if (${ICE_MACHINE} =~ gpsc3*) then
cat >> ${jobfile} << EOFB
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH -J ${ICE_CASENAME}
#SBATCH -A ${acct}
#SBATCH --partition ${queue}
#SBATCH --time ${batchtime}
#SBATCH --nodes ${nnodes}
#SBATCH --ntasks ${ntasks}
#SBATCH --cpus-per-task ${nthrds}
#SBATCH --mem-per-cpu=${batchmem}G
#SBATCH --comment="image=eccc/eccc_all_default_ubuntu-18.04-amd64_latest"
EOFB


else if (${ICE_MACHINE} =~ freya* ) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -j oe
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${batchtime}
EOFB

else if (${ICE_MACHINE} =~ boreas* ) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -j oe
#PBS -q ${queue}
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${batchtime}
EOFB

else if (${ICE_MACHINE} =~ gaeac5*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH --partition=batch
#SBATCH --qos=${queue}
#SBATCH --account=${acct}
#SBATCH --clusters=c5
#SBATCH --time=${batchtime}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=${taskpernodelimit}
#SBATCH --cpus-per-task=${nthrds}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
##SBATCH --mail-type FAIL
##SBATCH --mail-user=xxx@noaa.gov
EOFB

else if (${ICE_MACHINE} =~ gaeac6*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH --partition=batch
#SBATCH --qos=${queue}
#SBATCH --account=${acct}
#SBATCH --clusters=c6
#SBATCH --time=${batchtime}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=${taskpernodelimit}
#SBATCH --cpus-per-task=${nthrds}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
##SBATCH --mail-type FAIL
##SBATCH --mail-user=xxx@noaa.gov
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

else if (${ICE_MACHINE} =~ ursa*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH --partition=u1-compute
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

else if (${ICE_MACHINE} =~ wcoss2*) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -o ${ICE_CASENAME}
#PBS -j oe 
#PBS -A ICE-DEV
#PBS -l walltime=${batchtime}
##PBS -l select=${nnodes}:ncpus=${taskpernodelimit}
##PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l place=vscatter,select=${nnodes}:ncpus=${corespernode}:mpiprocs=${corespernode}:mem=256M
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
