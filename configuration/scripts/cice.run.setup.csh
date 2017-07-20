#! /bin/csh -f

echo ${0}

source ./cice.settings
source ${CICE_CASEDIR}/env.${CICE_MACHINE} || exit 2

set jobfile = cice.run
set subfile = cice.submit

set nthrds = ${CICE_NTHRDS}

#==========================================

# Write the batch code into the job file
${CICE_SCRDIR}/cice.batch.csh ${jobfile}
if ($? == -1) then
  exit -1
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

if !(-e \${CICE_RUNDIR}/ice.restart_file) cp \${CICE_RSTPFILE} \${CICE_RUNDIR}

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

# Write the job launching logic into the job file
${CICE_SCRDIR}/cice.launch.csh ${jobfile}
if ($? == -1) then
  exit -1
endif

#==========================================

cat >> ${jobfile} << EOFE
echo "CICE run finished: \`date\`"
echo " "

#--------------------------------------------

if !(-d \${CICE_LOGDIR}) mkdir -p \${CICE_LOGDIR}
cp -p \${CICE_RUNLOG_FILE} \${CICE_LOGDIR}

grep ' CICE COMPLETED SUCCESSFULLY' \${CICE_RUNLOG_FILE}
if ( \$? != 0 ) then
  echo "CICE run did not complete - see \${CICE_LOGDIR}/\${CICE_RUNLOG_FILE}"
  echo "\`date\` \${0}: \${CICE_CASENAME} run did NOT complete \${CICE_RUNLOG_FILE}"  >> \${CICE_CASEDIR}/README.case
  exit -1
endif

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
