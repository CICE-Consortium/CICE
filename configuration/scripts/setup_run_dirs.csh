#! /bin/csh -f

source ./cice.settings
source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} -nomodules || exit 2

if !(-d ${ICE_RUNDIR}) then
  echo "mkdir ${ICE_RUNDIR}"
  mkdir -p ${ICE_RUNDIR}
endif
if !(-d ${ICE_HSTDIR}) mkdir -p ${ICE_HSTDIR}
if !(-d ${ICE_RSTDIR}) mkdir -p ${ICE_RSTDIR}

exit 0
