#! /bin/csh -f

source ./cice.settings
source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

if (${ICE_MACHINE} != ${ICE_MACHINE_ENVNAME}) then
  echo "WARNING, is ICE_MACHINE setting OK, ${ICE_MACHINE}, ${ICE_MACHINE_ENVNAME}"
endif

if (${ICE_COMPILER} != ${ICE_MACHINE_COMPILER}) then
  echo "WARNING, is ICE_COMPILER setting OK, ${ICE_COMPILER}, ${ICE_MACHINE_COMPILER}"
endif

echo " "

if !(-d ${ICE_OBJDIR}) mkdir -p ${ICE_OBJDIR}
cd ${ICE_OBJDIR}

setenv ICE_CPPDEFS " "

if (${ICE_IOTYPE} == 'netcdf') then
  set IODIR = io_netcdf
  setenv ICE_CPPDEFS "${ICE_CPPDEFS} -Dncdf"
else if (${ICE_IOTYPE} == 'pio') then
  set IODIR = io_pio
  setenv ICE_CPPDEFS "${ICE_CPPDEFS} -Dncdf"
else
  set IODIR = io_binary
endif

### List of source code directories (in order of importance).
cat >! Filepath << EOF
${ICE_SANDBOX}/cicecore/drivers/${ICE_DRVOPT}
${ICE_SANDBOX}/cicecore/cicedynB/dynamics
${ICE_SANDBOX}/cicecore/cicedynB/general
${ICE_SANDBOX}/cicecore/cicedynB/analysis
${ICE_SANDBOX}/cicecore/cicedynB/infrastructure
${ICE_SANDBOX}/cicecore/cicedynB/infrastructure/io/$IODIR
${ICE_SANDBOX}/cicecore/cicedynB/infrastructure/comm/${ICE_COMMDIR}
${ICE_SANDBOX}/cicecore/shared
${ICE_SANDBOX}/icepack/columnphysics
EOF

if !(-d ${ICE_RUNDIR}) mkdir -p ${ICE_RUNDIR}

echo "make $*"
${ICE_MACHINE_MAKE} -j ${ICE_MACHINE_BLDTHRDS} VPFILE=Filepath EXEC=${ICE_RUNDIR}/cice \
    -f  ${ICE_CASEDIR}/Makefile MACFILE=${ICE_CASEDIR}/Macros.${ICE_MACHCOMP} $*
