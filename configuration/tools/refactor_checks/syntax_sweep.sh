#!/bin/bash
# gfortran -fsyntax-only sweep over the CICE source.
#
# This is not a build: there is no compiler on this machine that matches the
# target, and no way to run the model here.  What it does catch is everything
# the Fortran front end can see -- interface mismatches, rank and type errors,
# undeclared names, bad USE lists -- which is most of what a refactor breaks.
#
# ifort still catches four classes this misses, every one of them found the
# hard way on boreas: a constant out-of-bounds subscript (gfortran warns and
# exits 0), a variable subscript that evaluates to 0, a duplicated PUBLIC, and
# a reference to an unallocated allocatable.  A clean sweep here is necessary,
# never sufficient.
#
# usage: ./syntax_sweep.sh [comm] [io]      e.g. ./syntax_sweep.sh x io_binary
set -u
COMM=${1:-serial}
IO=${2:-io_binary}

# This script lives in <cice>/configuration/tools/refactor_checks.
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
MODS=$WORK/mods; mkdir -p "$MODS"

DIRS="
$HERE
$ROOT/icepack/columnphysics
$ROOT/cicecore/shared
$ROOT/cicecore/cicedyn/infrastructure/comm
$ROOT/cicecore/cicedyn/infrastructure/io/$IO
$ROOT/cicecore/cicedyn/general
$ROOT/cicecore/cicedyn/analysis
$ROOT/cicecore/cicedyn/dynamics
$ROOT/cicecore/cicedyn/infrastructure
$ROOT/cicecore/drivers/standalone/cice
"

SRC=""
for d in $DIRS; do
   [ -d "$d" ] || continue
   SRC="$SRC $(find "$d" -maxdepth 1 -name '*.F90' | sort)"
done

TOTAL=$(echo $SRC | wc -w)
REMAIN="$SRC"
OK=0

# Repeat until a pass compiles nothing new: this resolves the dependency order
# without having to know it.
while : ; do
   PROGRESS=0
   NEXT=""
   for f in $REMAIN; do
      if gfortran -fsyntax-only -ffree-line-length-none -fallow-argument-mismatch \
            -J"$MODS" -I"$MODS" "$f" > "$WORK/$(basename "$f").err" 2>&1; then
         OK=$((OK+1)); PROGRESS=1
      else
         NEXT="$NEXT $f"
      fi
   done
   REMAIN="$NEXT"
   [ "$PROGRESS" -eq 0 ] && break
   [ -z "$REMAIN" ] && break
done

echo "=== sweep comm=$COMM io=$IO : $OK / $TOTAL compiled ==="
if [ -n "${REMAIN// /}" ]; then
   for f in $REMAIN; do
      echo "--- FAILED: $(basename "$f")"
      # gfortran puts the line number, the source line and the caret ABOVE the
      # Error: line, so grepping for "Error" alone throws away the only part
      # that says where. Show the whole first few diagnostics instead.
      head -40 "$WORK/$(basename "$f").err"
   done
   exit 1
fi
