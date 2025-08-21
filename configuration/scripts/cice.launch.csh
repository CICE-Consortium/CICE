#!/bin/csh -f

#echo ${0}
echo "running cice.launch.csh"

#source ./cice.settings
#source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = $1

source ${ICE_SCRIPTS}/setup_machparams.csh

#==========================================
if (${ICE_MACHCOMP} =~ cheyenne*) then
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
else if (${ICE_MACHCOMP} =~ derecho*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpiexec --cpu-bind depth -n ${ntasks} -ppn ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHCOMP} =~ gust*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpiexec --cpu-bind depth -n ${ntasks} -ppn ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHCOMP} =~ gadi*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpirun -n ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHCOMP} =~ hobart* || ${ICE_MACHCOMP} =~ izumi*) then
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
else if (${ICE_MACHCOMP} =~ gaffney* || ${ICE_MACHCOMP} =~ koehr* || ${ICE_MACHCOMP} =~ mustang*) then
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
else if (${ICE_MACHCOMP} =~ nrlssc*) then
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
else if (${ICE_MACHCOMP} =~ narwhal_*hpcx*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} -hostfile \$PBS_NODEFILE \${EXTRA_OMPI_SETTINGS} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ onyx* || ${ICE_MACHCOMP} =~ narwhal*) then
cat >> ${jobfile} << EOFR
aprun -q -n ${ntasks} -N ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ carpenter*) then 
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else

if (${ICE_ENVNAME} =~ intelimpi* || ${ICE_ENVNAME} =~ gnuimpi*) then
cat >> ${jobfile} << EOFR
mpiexec -n ${ntasks} -ppn ${taskpernodelimit} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
mpiexec --cpu-bind depth -n ${ntasks} -ppn ${taskpernodelimit} -d ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

endif

#=======
else if (${ICE_MACHCOMP} =~ cori* || ${ICE_MACHCOMP} =~ perlmutter* || ${ICE_MACHCOMP} =~ blueback*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
#./cice >&! \$ICE_RUNLOG_FILE
srun --cpu-bind=cores --kill-on-bad-exit ./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
srun --cpu-bind=cores -n ${ntasks} -c ${nthrds} --kill-on-bad-exit ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHCOMP} =~ compy*) then
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
else if (${ICE_MACHCOMP} =~ chicoma*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHCOMP} =~ discover*) then
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
else if (${ICE_MACHCOMP} =~ fram*) then
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
else if (${ICE_MACHCOMP} =~ cesium*) then
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
else if (${ICE_MACHCOMP} =~ millikan*) then
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
else if (${ICE_MACHCOMP} =~ daley* || ${ICE_MACHCOMP} =~ banting*) then
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
else if (${ICE_MACHCOMP} =~ ppp5* || ${ICE_MACHCOMP} =~ ppp6* || ${ICE_MACHCOMP} =~ robert* || ${ICE_MACHCOMP} =~ underhill*) then
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
else if (${ICE_MACHCOMP} =~ ppp3*) then
if (${ICE_COMMDIR} =~ serial*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR
else
cat >> ${jobfile} << EOFR
rumpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
EOFR
endif

#=======
else if (${ICE_MACHCOMP} =~ gpsc3*) then
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
else if (${ICE_MACHCOMP} =~ freya*) then
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
else if (${ICE_MACHCOMP} =~ boreas*) then
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
else if (${ICE_MACHCOMP} =~ gaea*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ hera* || ${ICE_MACHCOMP} =~ ursa*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ orion*) then
cat >> ${jobfile} << EOFR
srun -n ${ntasks} -c ${nthrds} ./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ high_Sierra*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
#./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ phase3*) then
cat >> ${jobfile} << EOFR
mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
#./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ testmachine*) then
cat >> ${jobfile} << EOFR
./cice >&! \$ICE_RUNLOG_FILE
EOFR

#=======
else if (${ICE_MACHCOMP} =~ travisCI*) then
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
else if (${ICE_MACHCOMP} =~ conda*) then
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
  echo "${0} ERROR ${ICE_MACHCOMP} unknown"
  exit -1
endif
#=======

exit 0
