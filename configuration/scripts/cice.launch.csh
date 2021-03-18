#!/bin/csh -f

#echo ${0}
echo "running cice.launch.csh"

#source ./cice.settings
#source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = $1

set ntasks = ${ICE_NTASKS}
set nthrds = ${ICE_NTHRDS}
set maxtpn = ${ICE_MACHINE_TPNODE}

@ ncores = ${ntasks} * ${nthrds}
@ taskpernode = ${maxtpn} / $nthrds
if (${taskpernode} == 0) set taskpernode = 1
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}
@ corespernode = ${taskpernodelimit} * ${nthrds}

#==========================================
if (${ICE_MACHINE} =~ cheyenne*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpiexec_mpt -np ${ntasks} omplace ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ hobart* || ${ICE_MACHINE} =~ izumi*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpiexec -n ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ gaffney* || ${ICE_MACHINE} =~ koehr* || ${ICE_MACHINE} =~ mustang*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpiexec_mpt -np ${ntasks} omplace ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ onyx*) then
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ gordon* || ${ICE_MACHINE} =~ conrad*) then
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ cori*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
srun --cpu-bind=cores ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ compy*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
srun --mpi=pmi2 --kill-on-bad-exit --cpu-bind=cores ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ badger*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ fram*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ cesium*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ millikan*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ daley* || ${ICE_MACHINE} =~ banting*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ freya*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
aprun -n 1 -N 1 -d 1 ./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHINE} =~ hera*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ orion*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ high_Sierra*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
#./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ phase3*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
#./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ testmachine*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHINE} =~ travisCI*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif
#cat >> ${jobfile} << EOFR
#srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
#EOFR

#=======
else if (${ICE_MACHINE} =~ conda*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif


#=======
else
  echo "${0} ERROR ${ICE_MACHINE} unknown"
  exit -1
endif
#=======

exit 0
