# refactor_checks

Three checks that caught real faults during the compressed-layout refactor of
the dynamics, kept because the faults they catch are ones the compiler and the
test suite both miss, and because the remaining C/CD and VP conversions are
exposed to exactly the same classes.

None of these replace running the test suite. They are cheap things to do
*before* spending a cycle on a machine.

## syntax_sweep.sh

A `gfortran -fsyntax-only` sweep over the whole source, in dependency order.

    ./syntax_sweep.sh [comm] [io]        # e.g. ./syntax_sweep.sh x io_binary

It compiles nothing and runs nothing; it catches what the Fortran front end can
see — interface mismatches, rank and type errors, undeclared names, bad USE
lists — which is most of what a refactor breaks. `io_binary` covers every
module under `cicedyn/`; `io_netcdf` and `io_pio2` additionally need `netcdf`
and `pio` module stubs, which are not provided here.

`mpi_stub.F90` is a stub `mpi` module, needed because `comm/` uses `use mpi`
unconditionally. Without it the sweep stops at `ice_communicate` and silently
leaves every module downstream of it unchecked — which is most of the code a
dynamics refactor touches. It declares only the names the source references
and is never linked or run.

**A clean sweep is necessary, not sufficient.** ifort catches four classes this
misses, each found the hard way on a machine: a constant out-of-bounds
subscript (gfortran warns and exits 0), a variable subscript that evaluates to
zero, a duplicated PUBLIC, and a reference to an unallocated allocatable.

## arg_align.py

Checks call sites against subroutine dummy lists, positionally.

    ./arg_align.py [--all] FILE.F90 ...

It exists because of this, in `ice_dyn_evp1d`:

    dummy list:  ... G_forcexU , G_forceyU , G_Tbu   , G_umass , G_uvel , G_vvel
    call site:   ... G_forcexU , G_forceyU , G_umass , G_Tbu   , G_uvel , G_vvel

Both are `real(dbl_kind), dimension(:,:), intent(in)`, so the swap is
type-correct and **no compiler can see it** — not gfortran, not ifort with
every check enabled. `umass` then received `Tbu`, which is zero in deep water,
and the 1-D solver ran with zero ice mass at every active cell. It took three
machine cycles and a purpose-built runtime audit to find.

By default it reports only *transpositions* — the actual's own name appearing
at a different position in the dummy list — and argument-count mismatches.
Deliberate renames (`nActive` → `na0`) are suppressed, because there are
hundreds of them and they bury the signal; `--all` shows them. Count
mismatches are usually omitted OPTIONAL dummies, which is legal and common
here, so check before believing one.

It is a naming heuristic, not a type checker: it only works where both sides
use related names, which in this code base they nearly always do.

## compare_pending.sh

Recovers result lines for test cases whose run finished but was never recorded
— a batch job killed after the model completed but before the test script
appended its results leaves `PEND ... run` in `test_output`, indistinguishable
from a test that never ran. The restart data is on disk, so the comparisons
can simply be redone.

    ./compare_pending.sh                 # from testsuite.<name>, report only
    ./compare_pending.sh --apply         # also append the recovered lines

**Read a recovered FAIL with care.** The script cannot tell a job killed
*after* the model finished — where the restart is good and the result is
genuinely recoverable — from one killed *before* it finished, which leaves a
restart from the wrong model time with the right filename and a plausible
size. Check the runlog ended normally before trusting a recovered result, and
prefer re-running the case when it did not. Lines it writes are marked
`recovered` to record that provenance; the verdict is the leading status word.
