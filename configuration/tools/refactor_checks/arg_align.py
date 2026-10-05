#!/usr/bin/env python3
"""Check Fortran call sites against subroutine dummy lists, positionally.

Written after a transposed argument pair cost three boreas cycles.  In
ice_dyn_evp1d the dummy list read

    ... G_forcexU , G_forceyU , G_Tbu   , G_umass , G_uvel , G_vvel

while the call passed

    ... G_forcexU , G_forceyU , G_umass , G_Tbu   , G_uvel , G_vvel

Both are real(dbl_kind), dimension(:,:), intent(in), so the swap is
type-correct and NO compiler can see it -- not gfortran, not ifort, not with
every check enabled.  umass then received Tbu, which is zero in deep water,
and the 1d solver ran with zero ice mass at every active cell.

This is a naming heuristic, not a type checker: it pairs the Nth actual with
the Nth dummy and complains when the names disagree once a leading L_/G_ and
case are ignored.  So it only works where the two sides use related names --
which, in this code base, they almost always do.  Expect benign hits where
they deliberately differ (nActive -> na0); the ones to read are two dummies of
the SAME kind appearing in swapped order.

By default only transpositions and argument-count mismatches are reported.
Pass --all to see plain renames too.  Count mismatches are usually just
OPTIONAL dummies being omitted, which is legal and common here -- the
declarations are not parsed, so check before believing one.

usage: ./arg_align.py [--all] FILE.F90 [FILE.F90 ...]
"""
import re
import sys


def statements(lines):
    """Yield (index, logical line), joining Fortran free-form continuations.

    Comments are stripped first, and blank lines inside a continuation are
    tolerated -- gfortran accepts those and this code base contains some.
    """
    i = 0
    while i < len(lines):
        start = i
        buf = []
        while i < len(lines):
            code = lines[i].split('!')[0]
            buf.append(code)
            if not code.rstrip().endswith('&'):
                break
            i += 1
            while i < len(lines) and not lines[i].split('!')[0].strip():
                i += 1   # blank line inside a continuation
        yield start, ' '.join(buf)
        i += 1


def arglist(stmt, name):
    """The argument names of a call or definition of `name`, in order."""
    m = re.search(re.escape(name) + r'\s*\(', stmt, re.I)
    if not m:
        return None
    depth, out, cur = 0, [], ''
    for ch in stmt[m.end():]:
        if ch == '(':
            depth += 1
        elif ch == ')':
            if depth == 0:
                break
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    out.append(cur)
    return [a.replace('&', '').strip() for a in out if a.replace('&', '').strip()]


def norm(a):
    return re.sub(r'^[LG]_', '', a).lower()


def main(paths, show_all=False):
    bad = 0
    for path in paths:
        lines = open(path, errors='replace').read().split('\n')
        stmts = list(statements(lines))
        defs = {}
        for _, s in stmts:
            m = re.match(r'\s*(?:recursive\s+)?subroutine\s+(\w+)\s*\(', s, re.I)
            if m:
                defs[m.group(1).lower()] = arglist(s, m.group(1))
        for ln, s in stmts:
            m = re.search(r'\bcall\s+(\w+)\s*\(', s, re.I)
            if not m:
                continue
            name = m.group(1).lower()
            if name not in defs:
                continue           # defined elsewhere; nothing to compare
            actual, dummy = arglist(s, m.group(1)), defs[name]
            if actual is None or dummy is None:
                continue
            if any('=' in a for a in actual):
                continue           # keyword arguments are not positional
            if len(actual) != len(dummy):
                print(f'{path}:{ln+1}: {name}: {len(actual)} actual vs '
                      f'{len(dummy)} dummy arguments')
                bad += 1
                continue
            ndummy = [norm(d) for d in dummy]
            for k, (a, d) in enumerate(zip(actual, dummy), 1):
                if norm(a) == norm(d):
                    continue
                # A plain rename (nActive -> na0) is deliberate and everywhere.
                # A TRANSPOSITION is not: the actual's own name sits at some
                # other position in the dummy list, which means two arguments
                # of compatible type were written in the wrong order.  That is
                # the case no compiler can catch, so report only that unless
                # asked for everything.
                where = ndummy.index(norm(a)) + 1 if norm(a) in ndummy else None
                if where is None and not show_all:
                    continue
                note = f' -- but dummy {where} is {a!r}: TRANSPOSED?' if where else ''
                print(f'{path}:{ln+1}: {name}: position {k}: '
                      f'actual {a!r} -> dummy {d!r}{note}')
                bad += 1
    return 1 if bad else 0


if __name__ == '__main__':
    argv = sys.argv[1:]
    show_all = '--all' in argv
    argv = [a for a in argv if a != '--all']
    if not argv:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(argv, show_all))
