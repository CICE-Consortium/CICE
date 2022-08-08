#!/bin/csh -f

# inputs
# mpi tasks
set ntasks = ${ICE_NTASKS}
# threads
set nthrds = ${ICE_NTHRDS}
# max tasks per node
set maxtpn = ${ICE_MACHINE_TPNODE}
# batch charge account
set acct   = ${ICE_ACCOUNT}

# compute total cores needed and distribution of cores on nodes
# ncores = total cores needed (tasks * threads)
# taskpernode = number of MPI tasks per node based on size of node and threads
# nodes = number of total nodes needed based on tasks/threads
# taskpernodelimit = max(taskpernode, ntasks), when using less than 1 node
# corespernode = number of cores per node used
@ ncores = ${ntasks} * ${nthrds}
@ taskpernode = ${maxtpn} / $nthrds
if (${taskpernode} == 0) set taskpernode = 1
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}
@ corespernode = ${taskpernodelimit} * ${nthrds}

set runlength = ${ICE_RUNLENGTH}
if ($?ICE_MACHINE_MAXRUNLENGTH) then
  if (${runlength} > ${ICE_MACHINE_MAXRUNLENGTH}) then
    set runlength = ${ICE_MACHINE_MAXRUNLENGTH}
  endif
endif

set memuse = ${ICE_MEMUSE}
if ($?ICE_MACHINE_MAXMEMUSE) then
  if (${memuse} > ${ICE_MACHINE_MAXMEMUSE}) then
    set memuse = ${ICE_MACHINE_MAXMEMUSE}
  endif
endif

set queue = "${ICE_QUEUE}"
set batchtime = "00:15:00"
if (${runlength} == 0) set batchtime = "00:29:00"
if (${runlength} == 1) set batchtime = "00:59:00"
if (${runlength} == 2) set batchtime = "2:00:00"
if (${runlength} == 3) set batchtime = "3:00:00"
if (${runlength} == 4) set batchtime = "4:00:00"
if (${runlength} == 5) set batchtime = "5:00:00"
if (${runlength} == 6) set batchtime = "6:00:00"
if (${runlength} == 7) set batchtime = "7:00:00"
if (${runlength} >= 8) set batchtime = "8:00:00"
set batchmem = "5"
if (${memuse} == 1) set batchmem = "5"
if (${memuse} == 2) set batchmem = "10"
if (${memuse} == 3) set batchmem = "15"
if (${memuse} == 4) set batchmem = "20"
if (${memuse} == 5) set batchmem = "50"
if (${memuse} == 6) set batchmem = "100"
if (${memuse} == 7) set batchmem = "150"
if (${memuse} >= 8) set batchmem = "200"

set shortcase = `echo ${ICE_CASENAME} | cut -c1-15`

