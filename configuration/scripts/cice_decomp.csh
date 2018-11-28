#!/bin/csh -f

#--- inputs ---

#echo "${0:t}  input ICE_DECOMP_GRID  = $ICE_DECOMP_GRID"
#echo "${0:t}  input ICE_DECOMP_NTASK = $ICE_DECOMP_NTASK"
#echo "${0:t}  input ICE_DECOMP_NTHRD = $ICE_DECOMP_NTHRD"
#echo "${0:t}  input ICE_DECOMP_BLCKX = $ICE_DECOMP_BLCKX"
#echo "${0:t}  input ICE_DECOMP_BLCKY = $ICE_DECOMP_BLCKY"
#echo "${0:t}  input ICE_DECOMP_MXBLCKS = $ICE_DECOMP_MXBLCKS"

set grid = $ICE_DECOMP_GRID
set task = $ICE_DECOMP_NTASK
set thrd = $ICE_DECOMP_NTHRD

if (${task} <= 0 || ${thrd} <= 0) then
  echo "${0:t}: ERROR task and thread must be gt 0"
  exit -9
endif

#--- computation ---

@ cicepes = ${task} * ${thrd}

if (${grid} == 'col') then
  set nxglob = 5
  set nyglob = 5
  if (${cicepes} <= 1) then
    set blckx = 5; set blcky = 5
  else 
    set blckx = 1; set blcky = 1
  endif

else if (${grid} == 'gbox128') then
  set nxglob = 128
  set nyglob = 128
  if (${cicepes} <= 1) then
    set blckx = 128; set blcky = 128
  else if (${cicepes} <= 8) then
    set blckx = 32; set blcky = 32
  else if (${cicepes} <= 32) then
    set blckx = 16; set blcky = 16
  else
    set blckx = 8; set blcky = 8
  endif

else if (${grid} == 'gbox80') then
  set nxglob = 80
  set nyglob = 80
  if (${cicepes} <= 1) then
    set blckx = 80; set blcky = 80
  else if (${cicepes} <= 8) then
    set blckx = 20; set blcky = 20
  else
    set blckx = 8; set blcky = 8
  endif

else if (${grid} == 'gx3') then
  set nxglob = 100
  set nyglob = 116
  if (${cicepes} <= 1) then
    set blckx = 100; set blcky = 116
  else if (${cicepes} <= 8) then
    set blckx = 25; set blcky = 29
  else if (${cicepes} <= 32) then
    set blckx = 5; set blcky = 29
  else
    set blckx = 5; set blcky = 4
  endif

else if (${grid} == 'gx1') then
  set nxglob = 320
  set nyglob = 384
  if (${cicepes} <= 16) then
    set blckx = 40; set blcky = 48
  else if (${cicepes} <= 20) then
    set blckx = 32; set blcky = 48
  else if (${cicepes} <= 32) then
    set blckx = 20; set blcky = 24
  else if (${cicepes} <= 40) then
    set blckx = 16; set blcky = 24
  else if (${cicepes} < 80) then
    set blckx = 10; set blcky = 16
  else if (${cicepes} == 80) then
    set blckx = 8;  set blcky = 16
  else
    set blckx = 10; set blcky = 12
  endif

else if (${grid} == 'tx1') then
  set nxglob = 360
  set nyglob = 240
  if (${cicepes} <= 16) then
    set blckx = 90; set blcky = 60
  else if (${cicepes} <= 64) then
    set blckx = 20; set blcky = 20
  else
    set blckx = 10; set blcky = 10
  endif

else
  echo "${0:t}: ERROR unknown grid ${grid}"
  exit -9
endif

# check and override
if (${ICE_DECOMP_BLCKX} > 0 && ${ICE_DECOMP_BLCKY} > 0) then
  set blckx = ${ICE_DECOMP_BLCKX}
  set blcky = ${ICE_DECOMP_BLCKY}
else if (${ICE_DECOMP_BLCKX} < 1 && ${ICE_DECOMP_BLCKY} < 1) then
  # continue, use values computed above
else
  echo "${0:t}: ERROR user defined blocksize illegal"
  exit -9
endif

@ bx = $nxglob / ${blckx}
if ($bx * ${blckx} != $nxglob) @ bx = $bx + 1
@ by = $nyglob / ${blcky}
if ($by * ${blcky} != $nyglob) @ by = $by + 1

@ m = ($bx * $by) / ${task}
if ($m * ${task} != $bx * $by) @ m = $m + 1
set mxblcks = $m

# override
if (${ICE_DECOMP_MXBLCKS} > 0) set mxblcks = ${ICE_DECOMP_MXBLCKS}

set decomp = 'cartesian'
set dshape = 'slenderX2'
if (${nxglob} % ${cicepes} != 0) set decomp = 'roundrobin'
if (${mxblcks} * ${blcky} * 2 < ${nyglob}) set decomp = 'roundrobin'

#--- outputs ---

setenv ICE_DECOMP_NXGLOB $nxglob
setenv ICE_DECOMP_NYGLOB $nyglob
setenv ICE_DECOMP_BLCKX  $blckx
setenv ICE_DECOMP_BLCKY  $blcky
setenv ICE_DECOMP_MXBLCKS $mxblcks
setenv ICE_DECOMP_DECOMP $decomp
setenv ICE_DECOMP_DSHAPE $dshape

#echo "${0:t} output ICE_DECOMP_NXGLOB   = $ICE_DECOMP_NXGLOB"
#echo "${0:t} output ICE_DECOMP_NYGLOB   = $ICE_DECOMP_NYGLOB"
#echo "${0:t} output ICE_DECOMP_BLCKX    = $ICE_DECOMP_BLCKX"
#echo "${0:t} output ICE_DECOMP_BLCKY    = $ICE_DECOMP_BLCKY"
#echo "${0:t} output ICE_DECOMP_MXBLCKS  = $ICE_DECOMP_MXBLCKS"
#echo "${0:t} output ICE_DECOMP_DECOMP   = $ICE_DECOMP_DECOMP"
#echo "${0:t} output ICE_DECOMP_DSHAPE   = $ICE_DECOMP_DSHAPE"

exit 0
