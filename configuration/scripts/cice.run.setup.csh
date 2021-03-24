#!/bin/csh -f

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
cp -f \${ICE_CASEDIR}/env.${ICE_MACHCOMP} \${ICE_RUNDIR}
cp -f \${ICE_CASEDIR}/cice.settings \${ICE_RUNDIR}
set diagtype = \`grep -i diag_type \${ICE_CASEDIR}/ice_in | grep -i stdout | wc -l\`
set diagfile = \`grep -i diag_file \${ICE_CASEDIR}/ice_in | sed -e "s/.* = '\(.*\)'/\1/"\`

echo " "
echo "CICE case directory is \${ICE_CASEDIR}"
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

set checkfile = \${ICE_RUNLOG_FILE}
cp -p \${ICE_RUNLOG_FILE} \${ICE_LOGDIR}
echo "CICE output file is \${ICE_LOGDIR}/\${ICE_RUNLOG_FILE}"
echo "\`date\` \${0}: CICE output file is \${ICE_LOGDIR}/\${ICE_RUNLOG_FILE}" >> \${ICE_CASEDIR}/README.case

if ( \${diagtype} == 0) then
  set checkfile = \${diagfile}
  cp -p \${diagfile} \${ICE_LOGDIR}
  echo "CICE output file is \${ICE_LOGDIR}/\${diagfile}"
  echo "\`date\` \${0}: CICE output file is \${ICE_LOGDIR}/\${diagfile}" >> \${ICE_CASEDIR}/README.case
endif

grep ' CICE COMPLETED SUCCESSFULLY' \${checkfile}
if ( \$status == 0 ) then
  echo "CICE run completed successfully"
  echo "\`date\` \${0}: CICE run completed successfully"  >> \${ICE_CASEDIR}/README.case
else
  grep 'COMPLETED SUCCESSFULLY' \${checkfile}
  if ( \$status == 0 ) then
    echo "Run completed successfully"
    echo "\`date\` \${0}: Run completed successfully"  >> \${ICE_CASEDIR}/README.case
  else
    echo "CICE run did NOT complete"
    echo "\`date\` \${0}: CICE run did NOT complete"  >> \${ICE_CASEDIR}/README.case
    exit -1
  endif
endif

if ( \${diagtype} == 0) then
  exit -1
endif

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
