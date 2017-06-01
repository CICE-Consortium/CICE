#! /bin/csh -f

echo ${0}

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.run
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

#==========================================

elseif (${CICE_MACHINE} =~ cheyenne*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
#PBS -m ae 
#PBS -V
#PBS -q regular
#PBS -N ${CICE_CASENAME}
#PBS -A ${CICE_ACCT}
#PBS -l select=${nnodes}:ncpus=${ntasks}:mpiprocs=${ntasks}
#PBS -l walltime=02:00:00
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
#SBATCH -t 0:45:00
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=low
EOFB

else if (${CICE_MACHINE} =~ pinto*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${CICE_CASENAME}
#SBATCH -t 0:45:00
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=standby
EOFB

endif

#==========================================

cat >> ${jobfile} << EOF1

#--------------------------------------------

cd ${CICE_CASEDIR}
source ./cice.settings || exit 2
source ./env.\${CICE_MACHINE} || exit 2

echo " "
echo "\${0}:"

set  stamp   = \`date '+%y%m%d-%H%M%S'\`
set CICE_RUNLOG_FILE = "cice.runlog.\${stamp}"

#--------------------------------------------

if !(-d \${CICE_RUNDIR}) mkdir -p \${CICE_RUNDIR}
if !(-d \${CICE_HSTDIR}) mkdir -p \${CICE_HSTDIR}
if !(-d \${CICE_RSTDIR}) mkdir -p \${CICE_RSTDIR}

if !(-e ice.restart_file) cp \${CICE_RSTPFILE} \${CICE_RUNDIR}

#--------------------------------------------
cd \${CICE_RUNDIR}

setenv OMP_NUM_THREADS ${nthrds}

cp -f \${CICE_CASEDIR}/ice_in \${CICE_RUNDIR}
echo " "
echo "CICE rundir is \${CICE_RUNDIR}"
echo "CICE log file is \${CICE_RUNLOG_FILE}"
echo "CICE run started : \`date\`"

EOF1

#==========================================

if (${CICE_MACHINE} =~ yellowstone*) then
cat >> ${jobfile} << EOFR
setenv MP_TASK_AFFINITY core:\${OMP_NUM_THREADS}
mpirun.lsf ./cice >&! \$CICE_RUNLOG_FILE
EOFR

elseif (${CICE_MACHINE} =~ cheyenne*) then
cat >> ${jobfile} << EOFR
mpiexec_mpt -n ${ntasks} ./cice >&! \$CICE_RUNLOG_FILE
EOFR

else if (${CICE_MACHINE} =~ thunder*) then
cat >> ${jobfile} << EOFR
mpiexec_mpt -np ${ntasks} omplace ./cice >&! \$CICE_RUNLOG_FILE
EOFR

else if (${CICE_MACHINE} =~ gordon* || ${CICE_MACHINE} =~ conrad*) then
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${ntasks} -d ${nthrds} ./cice >&! \$CICE_RUNLOG_FILE
EOFR

else if (${CICE_MACHINE} =~ cori*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$CICE_RUNLOG_FILE
EOFR

else if (${CICE_MACHINE} =~ wolf*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$CICE_RUNLOG_FILE
EOFR

else if (${CICE_MACHINE} =~ pinto*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$CICE_RUNLOG_FILE
EOFR
#cat >> ${jobfile} << EOFR
#srun -n ${ntasks} -c ${nthrds} ./cice >&! \$CICE_RUNLOG_FILE
#EOFR

endif

#==========================================

cat >> ${jobfile} << EOFE
echo "CICE run finished: \`date\`"
echo " "

#--------------------------------------------

if !(-d \${CICE_LOGDIR}) mkdir -p \${CICE_LOGDIR}
cp -p \${CICE_RUNLOG_FILE} \${CICE_LOGDIR}

grep ' CICE COMPLETED SUCCESSFULLY' \${CICE_RUNLOG_FILE}  || echo "CICE run did not complete - see \${CICE_LOGDIR}/\${CICE_RUNLOG_FILE}" && exit -1

echo "\`date\` \${0}: \${CICE_CASENAME} run completed \${CICE_RUNLOG_FILE}"  >> \${CICE_CASEDIR}/README.case
echo "done \${0}"

EOFE

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${CICE_MACHINE_SUBMIT} ${jobfile}
echo "\`date\` \${0}: ${CICE_CASENAME} job submitted"  >> ${CICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
