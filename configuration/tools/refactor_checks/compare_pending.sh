#!/bin/bash
#==============================================================================
# Recover the result lines for test cases whose run finished but was never
# recorded.
#
# The suite writes "PEND <test> run" at submission and replaces it with
# PASS/FAIL when the test script finishes.  If the batch job is killed after
# the model completes but before that script appends its results -- a
# wall-clock SIGTERM, typically -- the run is perfectly good, yet test_output
# still reads PEND and no "compare" line is ever written.  Those tests then
# look identical to ones that never ran.
#
# The restart data is already on disk, so the comparisons can simply be redone.
# This replicates what configuration/scripts/tests/baseline.script does, for
# exactly those cases, without re-running the model.
#
# Two things are recovered per case:
#   test    - restart vs continuous self-check   (comparebfb.csh <dir>)
#   compare - vs the baseline                    (comparebfb.csh <base> <test>)
#
# Usage, from the testsuite.<name> directory:
#   ./compare_pending.sh                      report only (writes nothing)
#   ./compare_pending.sh --apply              also append the result lines to
#                                             each case's test_output
#   ./compare_pending.sh --run-root DIR       look for restart data under
#                                             DIR/<case>/restart instead of
#                                             the ICE_RUNDIR in cice.settings
#
# Default is report-only on purpose: --apply edits your result files.
#
# --run-root exists because the run directories are not always where
# cice.settings says: the data can be relocated after the job finishes, so the
# recorded ICE_RUNDIR no longer resolves.  Point it at the directory holding
# the per-case run dirs, e.g.  --run-root /data/dktar/cice_original/run
#==============================================================================

set -u

apply=0
run_root=""
while [ $# -gt 0 ]; do
    case "$1" in
        --apply)    apply=1; shift ;;
        --run-root) run_root="${2:-}"; shift 2 ;;
        *) echo "usage: $0 [--apply] [--run-root DIR]" >&2; exit 2 ;;
    esac
done
if [ -n "$run_root" ] && [ ! -d "$run_root" ]; then
    echo "error: --run-root '$run_root' is not a directory" >&2; exit 2
fi

# comparebfb.csh exit codes
status_name () {
    case "$1" in
        0) echo PASS ;;
        1) echo FAIL ;;
        2) echo MISS ;;
        *) echo FAIL ;;   # 9 = usage error, reported as FAIL like baseline.script
    esac
}

n_seen=0; n_done=0; n_skip=0
declare -A tally=( [PASS]=0 [FAIL]=0 [MISS]=0 )

for case_dir in */; do
    case_dir=${case_dir%/}
    out="$case_dir/test_output"
    [ -f "$out" ] || continue

    # a case needing recovery has PEND ... run and no PASS/FAIL ... run
    pend=$(awk '$1=="PEND" && $3=="run" {print $2; exit}' "$out")
    [ -n "$pend" ] || continue
    if grep -qE "^(PASS|FAIL) ${pend} run" "$out"; then continue; fi

    n_seen=$((n_seen+1))

    settings="$case_dir/cice.settings"
    if [ ! -f "$settings" ]; then
        echo "SKIP $pend  (no cice.settings)"; n_skip=$((n_skip+1)); continue
    fi

    # cice.settings is csh.  Source it there and hand the values back, rather
    # than trying to parse setenv lines that reference each other.
    vals=$(csh -f -c "source '$settings' >& /dev/null; printf '%s\n' \
           \"\$ICE_RUNDIR\" \"\$ICE_BASELINE\" \"\$ICE_BASECOM\" \
           \"\$ICE_TESTNAME\" \"\$ICE_CASEDIR\" \"\$ICE_SPVAL\" \"\$ICE_TEST\"" 2>/dev/null)
    { read -r RUNDIR; read -r BASELINE; read -r BASECOM
      read -r TESTNAME; read -r CASEDIR; read -r SPVAL; read -r TEST
    } <<< "$vals"

    if [ -z "${RUNDIR:-}" ]; then
        echo "SKIP $pend  (could not read cice.settings)"; n_skip=$((n_skip+1)); continue
    fi

    cmp_csh="$CASEDIR/casescripts/comparebfb.csh"
    [ -x "$cmp_csh" ] || cmp_csh="$case_dir/casescripts/comparebfb.csh"
    if [ ! -f "$cmp_csh" ]; then
        echo "SKIP $pend  (no comparebfb.csh)"; n_skip=$((n_skip+1)); continue
    fi

    # What actually matters is whether the restart data still exists -- that
    # is what comparebfb compares.  The runlog is only informational, and it
    # lives in two places: $ICE_RUNDIR while the job runs, and <case>/logs
    # afterwards, since CICE copies it out at the end.  baseline.script only
    # ever looks in the rundir because it runs inside the job.
    # Prefer the recorded ICE_RUNDIR; fall back to --run-root/<case>/restart
    # when the data has been relocated since the job ran.
    rst="$RUNDIR/restart"
    nrst=$(ls -1 "$rst" 2>/dev/null | wc -l)
    if [ "$nrst" -eq 0 ] && [ -n "$run_root" ]; then
        alt="$run_root/$case_dir/restart"
        if [ -d "$alt" ]; then
            rst="$alt"
            nrst=$(ls -1 "$rst" 2>/dev/null | wc -l)
        fi
    fi
    if [ "$nrst" -eq 0 ]; then
        nlog=$(ls -1 "$RUNDIR"/cice.runlog* "$case_dir"/logs/cice.runlog* 2>/dev/null | wc -l)
        if [ "$nlog" -eq 0 ]; then
            echo "SKIP $pend"
            echo "     no restart data and no runlog -- this case never ran"
        else
            echo "SKIP $pend"
            echo "     runlog present but no restart data in $rst"
            echo "     rundir has been cleaned or purged: must re-run, cannot recover"
        fi
        n_skip=$((n_skip+1)); continue
    fi

    lines=()

    # --- restart self-check: the "test" line -------------------------------
    if [ "$TEST" = "restart" ]; then
        csh -f "$cmp_csh" "$rst" >/dev/null 2>&1
        s=$(status_name $?)
        lines+=("$s $TESTNAME test")
        tally[$s]=$(( ${tally[$s]} + 1 ))
    fi

    # --- baseline comparison: the "compare" line ---------------------------
    if [ -n "${BASECOM:-}" ] && [ "$BASECOM" != "$SPVAL" ]; then
        base_dir="$BASELINE/$BASECOM/$TESTNAME/restart"
        csh -f "$cmp_csh" "$base_dir" "$rst" >/dev/null 2>&1
        rc=$?
        s=$(status_name $rc)
        if [ "$s" != "PASS" ]; then
            echo "     base: $base_dir  ($(ls -1 "$base_dir" 2>/dev/null | wc -l) files)"
            echo "     test: $rst  ($(ls -1 "$rst" 2>/dev/null | wc -l) files)"
        fi
        # "recovered" marks the line's PROVENANCE -- reconstructed after the
        # fact rather than written by the job -- so it goes on every such line.
        # It used to be appended only when the comparison did not pass, which
        # read as "successfully recovered" and meant the exact opposite.  The
        # verdict is the leading status word; this word is not a verdict.
        lines+=("$s $TESTNAME compare $BASECOM -1 -1 -1 recovered")
        tally[$s]=$(( ${tally[$s]} + 1 ))
    fi

    if [ ${#lines[@]} -eq 0 ]; then
        echo "SKIP $pend  (nothing to compare -- no baseline requested)"
        n_skip=$((n_skip+1)); continue
    fi

    for l in "${lines[@]}"; do
        echo "$l"
        [ $apply -eq 1 ] && echo "$l" >> "$out"
    done
    n_done=$((n_done+1))
done

echo
echo "cases needing recovery : $n_seen"
echo "  compared             : $n_done"
echo "  skipped              : $n_skip"
echo "results: PASS=${tally[PASS]}  FAIL=${tally[FAIL]}  MISS=${tally[MISS]}"
if [ $apply -eq 0 ]; then
    echo
    echo "(report only -- nothing written.  Re-run with --apply to append"
    echo " these lines to each case's test_output.)"
fi
