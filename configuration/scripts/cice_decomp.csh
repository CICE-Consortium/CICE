#!/bin/csh -f

#--- inputs ---

echo "${0:t}  input CICE_DECOMP_GRID  = $CICE_DECOMP_GRID"
echo "${0:t}  input CICE_DECOMP_NTASK = $CICE_DECOMP_NTASK"
echo "${0:t}  input CICE_DECOMP_NTHRD = $CICE_DECOMP_NTHRD"

set grid = $CICE_DECOMP_GRID
set task = $CICE_DECOMP_NTASK
set thrd = $CICE_DECOMP_NTHRD

#--- computation ---

set gridfile = unknown
set kmtfile  = unknown
set initfile = unknown
set rstpfile = unknown

@ cicepes = ${task} * ${thrd}

if (${grid} == 'col') then
  set nxglob = 5
  set nyglob = 5
  if (${cicepes} <= 1) then
    set blckx = 5; set blcky = 5
  else 
    set blckx = 1; set blcky = 1
  endif

else if (${grid} == 'gx3') then
  if ($?CICE_SANDBOX) then
    set initfile = ${CICE_SANDBOX}/configuration/data/gx3/iced_gx3_v5.nc
    set gridfile = ${CICE_SANDBOX}/configuration/data/gx3/global_gx3.grid
    set kmtfile  = ${CICE_SANDBOX}/configuration/data/gx3/global_gx3.kmt
    set rstpfile = ${CICE_SANDBOX}/configuration/data/gx3/ice.restart_file
  endif
  set nxglob = 100
  set nyglob = 116
  if (${cicepes} <= 8) then
    set blckx = 25; set blcky = 29
  else if (${cicepes} <= 32) then
    set blckx = 5; set blcky = 29
  else
    set blckx = 5; set blcky = 5
  endif

else if (${grid} == 'gx1') then
  if ($?CICE_SANDBOX) then
    set initfile = ${CICE_SANDBOX}/configuration/data/gx1/iced_gx1_v5.nc
    set gridfile = ${CICE_SANDBOX}/configuration/data/gx1/global_gx1.grid
    set kmtfile  = ${CICE_SANDBOX}/configuration/data/gx1/global_gx1.kmt
    set rstpfile = ${CICE_SANDBOX}/configuration/data/gx1/ice.restart_file
  endif
  set nxglob = 320
  set nyglob = 384
  if (${cicepes} <= 16) then
    set blckx = 40; set blcky = 48
  else if (${cicepes} <= 64) then
    set blckx = 20; set blcky = 24
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

@ bx = $nxglob / ${blckx}
if ($bx * ${blckx} != $nxglob) @ bx = $bx + 1
@ by = $nyglob / ${blcky}
if ($by * ${blcky} != $nyglob) @ by = $by + 1

@ m = ($bx * $by) / ${task}
if ($m * ${task} != $bx * $by) @ m = $m + 1
set mxblcks = $m

set decomp = 'cartesian'
set dshape = 'slenderX2'
if (${nxglob} % ${cicepes} != 0) set decomp = 'roundrobin'

#--- outputs ---

setenv CICE_DECOMP_NXGLOB $nxglob
setenv CICE_DECOMP_NYGLOB $nyglob
setenv CICE_DECOMP_BLCKX  $blckx
setenv CICE_DECOMP_BLCKY  $blcky
setenv CICE_DECOMP_MXBLCKS $mxblcks
setenv CICE_DECOMP_DECOMP $decomp
setenv CICE_DECOMP_DSHAPE $dshape
setenv CICE_DECOMP_GRIDFILE $gridfile
setenv CICE_DECOMP_KMTFILE  $kmtfile
setenv CICE_DECOMP_INITFILE $initfile
setenv CICE_DECOMP_RSTPFILE $rstpfile

echo "${0:t} output CICE_DECOMP_NXGLOB   = $CICE_DECOMP_NXGLOB"
echo "${0:t} output CICE_DECOMP_NYGLOB   = $CICE_DECOMP_NYGLOB"
echo "${0:t} output CICE_DECOMP_BLCKX    = $CICE_DECOMP_BLCKX"
echo "${0:t} output CICE_DECOMP_BLCKY    = $CICE_DECOMP_BLCKY"
echo "${0:t} output CICE_DECOMP_MXBLCKS  = $CICE_DECOMP_MXBLCKS"
echo "${0:t} output CICE_DECOMP_DECOMP   = $CICE_DECOMP_DECOMP"
echo "${0:t} output CICE_DECOMP_DSHAPE   = $CICE_DECOMP_DSHAPE"
echo "${0:t} output CICE_DECOMP_GRIDFILE = $CICE_DECOMP_GRIDFILE"
echo "${0:t} output CICE_DECOMP_KMTFILE  = $CICE_DECOMP_KMTFILE"
echo "${0:t} output CICE_DECOMP_INITFILE = $CICE_DECOMP_INITFILE"
echo "${0:t} output CICE_DECOMP_RSTPFILE = $CICE_DECOMP_RSTPFILE"

exit 0
