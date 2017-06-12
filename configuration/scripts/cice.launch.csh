#! /bin/csh -f

echo ${0}

source ./cice.settings

set jobfile = $1

set ntasks = ${CICE_NTASKS}
set nthrds = ${CICE_NTHRDS}

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
