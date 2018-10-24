#! /bin/csh -f

#echo ${0}
echo "running cice.run.setup.csh"

#source ./cice.settings
#source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = cice.run
set subfile = cice.submit

set nthrds = ${ICE_NTHRDS}

#==========================================

# Write the batch code into the job file
${ICE_SCRIPTS}/cice.batch.csh ${jobfile}
if ($status != 0) then
  echo "${0}: ERROR cice.batch.csh aborted"
  exit -1
endif

#==========================================

cat >> ${jobfile} << EOF1

#--------------------------------------------

cd ${ICE_CASEDIR}
source ./cice.settings || exit 2
source ./env.\${ICE_MACHCOMP} || exit 2

echo " "
echo "\${0}:"

set  stamp   = \`date '+%y%m%d-%H%M%S'\`
set ICE_RUNLOG_FILE = "cice.runlog.\${stamp}"

#--------------------------------------------

./setup_run_dirs.csh

#--------------------------------------------
cd \${ICE_RUNDIR}

setenv OMP_NUM_THREADS ${nthrds}

cp -f \${ICE_CASEDIR}/ice_in \${ICE_RUNDIR}
echo " "
echo "CICE rundir is \${ICE_RUNDIR}"
echo "CICE log file is \${ICE_RUNLOG_FILE}"
echo "CICE run started : \`date\`"

EOF1

#==========================================

# Write the job launching logic into the job file
${ICE_SCRIPTS}/cice.launch.csh ${jobfile}
if ($status != 0) then
  echo "${0}: ERROR cice.launch.csh aborted"
  exit -1
endif

#==========================================

cat >> ${jobfile} << EOFE
echo "CICE run finished: \`date\`"
echo " "

#--------------------------------------------

if !(-d \${ICE_LOGDIR}) mkdir -p \${ICE_LOGDIR}
cp -p \${ICE_RUNLOG_FILE} \${ICE_LOGDIR}

grep ' CICE COMPLETED SUCCESSFULLY' \${ICE_RUNLOG_FILE}
if ( \$status != 0 ) then
  echo "CICE run did not complete - see \${ICE_LOGDIR}/\${ICE_RUNLOG_FILE}"
  echo "\`date\` \${0}: \${ICE_CASENAME} run did NOT complete \${ICE_RUNLOG_FILE}"  >> \${ICE_CASEDIR}/README.case
  exit -1
endif

echo "\`date\` \${0}: \${ICE_CASENAME} run completed \${ICE_RUNLOG_FILE}"  >> \${ICE_CASEDIR}/README.case
echo "done \${0}"

EOFE

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${ICE_MACHINE_SUBMIT} ./${jobfile}
echo "\`date\` \${0}: ${ICE_CASENAME} job submitted"  >> ${ICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
exit 0
