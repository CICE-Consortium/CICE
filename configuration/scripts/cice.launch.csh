#! /bin/csh -f

#echo ${0}
echo "running cice.launch.csh"

source ./cice.settings

set jobfile = $1

set ntasks = ${ICE_NTASKS}
set nthrds = ${ICE_NTHRDS}

#==========================================

if (${ICE_MACHINE} =~ cheyenne*) then
cat >> ${jobfile} << EOFR
mpiexec_mpt -n ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ thunder*) then
cat >> ${jobfile} << EOFR
mpiexec_mpt -np ${ntasks} omplace ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ onyx*) then
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${ntasks} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ gordon* || ${ICE_MACHINE} =~ conrad*) then
cat >> ${jobfile} << EOFR
aprun -n ${ntasks} -N ${ntasks} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ cori*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ wolf*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ pinto*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ fram*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ testmachine*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR

else if (${ICE_MACHINE} =~ travisCI*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#cat >> ${jobfile} << EOFR
#srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
#EOFR

else
  echo "${0} ERROR ${ICE_MACHINE} unknown"
  exit -1
endif

exit 0
