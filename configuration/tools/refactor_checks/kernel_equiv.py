#!/usr/bin/env python3
"""Check that a compressed kernel is arithmetically identical to its block original.

The compressed kernels in this refactor are mechanical rewrites: the same
expressions with block subscripts replaced by compressed slots. That is exactly
the kind of edit where a single wrong neighbour -- isW where isS belongs --
compiles cleanly, runs cleanly, and is wrong only in the last digits of a field
nobody prints. The suite catches it, but a suite cycle is hours and this is
seconds.

So: extract both loop bodies, substitute the slot forms back to their block
equivalents, strip whitespace, and compare. Anything other than IDENTICAL is
either a real divergence or a substitution this script does not know about --
both worth looking at.

    ./kernel_equiv.py FILE.F90 block_name compressed_name [-m SLOT=BLOCK ...]

The default map covers the directions used so far. Add more with -m, e.g.

    ./kernel_equiv.py ice_dyn_shared.F90 strain_rates_Tdt \\
        strain_rates_Tdt_compressed -m '_s(isSW)=(ijb-nx_block-1)'

This proves the ARITHMETIC matches. It cannot tell you the neighbour table is
right -- that nbr(is,nbr_W) really is the cell at ijb-1 -- which is a property
of init_ragged_layout, not of the kernel.
"""
import argparse
import re
import sys

# slot form -> block form, for the stencils converted so far
DEFAULT_MAP = [
    ('_s(is)',    '(ijb)'),
    ('_s(isW)',   '(ijb-1)'),
    ('_s(isE)',   '(ijb+1)'),
    ('_s(isS)',   '(ijb-nx_block)'),
    ('_s(isN)',   '(ijb+nx_block)'),
    ('_s(isSW)',  '(ijb-nx_block-1)'),
    ('_s(isSE)',  '(ijb-nx_block+1)'),
    ('_s(isNW)',  '(ijb+nx_block-1)'),
    ('_s(isNE)',  '(ijb+nx_block+1)'),
]

# lines that exist only in the compressed form and have no block counterpart
DERIVATIONS = re.compile(
    r'is\w*\s*=\s*(Xmap\(ijb\)|nbr\(is,\s*nbr_\w+\s*\))')


def body(src, name):
    """The executable body of `name`, comments stripped, as one string."""
    m = re.search(r'subroutine\s+' + re.escape(name) + r'\s*\(', src)
    if not m:
        sys.exit(f'error: no subroutine {name}')
    # anchor the end on a line of its own: 'end subroutine foo' is a prefix of
    # 'end subroutine foobar', and in this code base that really happens
    e = re.search(r'end subroutine\s+' + re.escape(name) + r'\s*$',
                  src[m.start():], re.M)
    if not e:
        sys.exit(f'error: no end subroutine {name}')
    text = src[m.start():m.start() + e.start()]
    loop = re.search(r'^\s*do ij\b.*$', text, re.M)
    if not loop:
        sys.exit(f'error: no "do ij" loop in {name}')
    text = text[loop.start():]
    return ' '.join(l.split('!')[0].strip() for l in text.split('\n'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('block')
    ap.add_argument('compressed')
    ap.add_argument('-m', '--map', action='append', default=[],
                    metavar='SLOT=BLOCK', help='extra substitution')
    a = ap.parse_args()

    src = open(a.file, errors='replace').read()
    blk = body(src, a.block)
    cmp_ = body(src, a.compressed)

    subs = list(DEFAULT_MAP)
    for extra in a.map:
        if '=' not in extra:
            sys.exit(f'error: -m wants SLOT=BLOCK, got {extra!r}')
        k, v = extra.split('=', 1)
        subs.insert(0, (k, v))

    norm = cmp_
    for slot, block in subs:
        norm = norm.replace(slot, block)
    norm = DERIVATIONS.sub('', norm)

    squash = lambda t: re.sub(r'\s+', '', t)
    A, B = squash(blk), squash(norm)

    if A == B:
        print(f'IDENTICAL  {a.block} == {a.compressed}')
        return 0

    print(f'DIFFERS    {a.block} != {a.compressed}')
    for i, (x, y) in enumerate(zip(A, B)):
        if x != y:
            print(f'  first difference at character {i}')
            print(f'    block      ...{A[max(0, i-60):i+60]}')
            print(f'    compressed ...{B[max(0, i-60):i+60]}')
            break
    else:
        print(f'  common prefix, lengths differ: {len(A)} vs {len(B)}')
        j = min(len(A), len(B))
        print(f'    block      tail: {A[j:][:120]}')
        print(f'    compressed tail: {B[j:][:120]}')
    return 1


if __name__ == '__main__':
    sys.exit(main())
