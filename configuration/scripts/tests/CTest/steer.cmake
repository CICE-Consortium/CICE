# -----------------------------------------------------------  
# -- Get environment
# -----------------------------------------------------------  

## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)

set(CTEST_SITE                          "${HOSTNAME}")

## -- Set site / build name
## --------------------------

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

find_program(GIT_CMD NAMES git)
exec_program(${GIT_CMD} ARGS rev-parse --short HEAD OUTPUT_VARIABLE GIT_COMMIT_HASH)

find_program(IFORT_CMD NAMES ifort)
exec_program(${IFORT_CMD} ARGS --version | head -n 1 | awk '{print $3}' OUTPUT_VARIABLE COMPILER_VERSION)

set(CTEST_BUILD_NAME        "${osname}-${cpu}-intel${COMPILER_VERSION}-${GIT_COMMIT_HASH}")

message("build name = ${CTEST_BUILD_NAME}")

set(CTEST_DASHBOARD_ROOT   "$ENV{PWD}")
set(CTEST_SOURCE_DIRECTORY "$ENV{PWD}")
set(CTEST_BINARY_DIRECTORY "$ENV{PWD}")
message("source directory = ${CTEST_SOURCE_DIRECTORY}")

ctest_start(${MODEL} TRACK ${MODEL})
ctest_test( BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

ctest_submit(           RETURN_VALUE res)
