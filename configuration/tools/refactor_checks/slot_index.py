#!/usr/bin/env python3
"""Find dummies that receive a compressed (_slot) actual but are read block-wise.

The third fault class of this refactor, after transposed arguments and
dead-slot arithmetic, and the one no existing check could see.

In strain_rates_U the dummies uvelU and vvelU were declared

    real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: uvelU, vvelU

and read as uvelU(ijb).  The driver was later changed to pass uvel_slot(...)
-- a compressed array indexed by SLOT, not by block position.  An
explicit-shape dummy is not bounds-checked against its actual, so the read ran
off the end of the slice in silence, DeltaU came out NaN, and the abort landed
three routines away inside visc_replpress.

arg_align.py cannot catch it: uvel_slot -> uvelU is a rename, not a
transposition.  The mismatch is one of LAYOUT, not of order.

The rule this encodes: if a call passes `something_slot` at position N, then
the dummy at position N must be indexed by a slot variable (is, isW, isNE, ...)
and never by ijb or (i,j).

    ./slot_index.py FILE.F90 [FILE.F90 ...]        layout check (above)
    ./slot_index.py --pack FILE.F90 [...]          pack/unpack symmetry

The --pack mode catches a sibling: a slot array UNPACKED after a subcycle but
never PACKED before it. 59334cd switched strain_rates_U to read uvel_slot and
vvel_slot without adding them to the C/CD pack lists, so the first subcycle
pass read whatever the previous timestep had left there -- the right answer in
a continuous run, zeros after a restart. The restart test stopped reproducing
and the cause was four commits upstream of the symptom.

Pure outputs legitimately appear unpacked-but-not-packed (strint*, taub*): they
are written every pass before anything reads them. The check lists what it
finds; the question to ask of each is "is this read before it is written?".

Both the call and the routine must live in the files given; cross-file calls
are skipped silently, so pass the whole dynamics directory to cover them.
"""
import re
import sys

SLOT_IX = re.compile(r'\b(is|is[A-Z][A-Za-z]*)\b')
BLOCK_IX = re.compile(r'\(\s*ijb\b|\(\s*i\s*,\s*j\s*\)')


def statements(lines):
    i = 0
    while i < len(lines):
        start, buf = i, []
        while i < len(lines):
            code = lines[i].split('!')[0]
            buf.append(code)
            if not code.rstrip().endswith('&'):
                break
            i += 1
            while i < len(lines) and not lines[i].split('!')[0].strip():
                i += 1
        yield start, ' '.join(buf)
        i += 1


def arglist(stmt, name):
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


def main(paths):
    src = {}
    for p in paths:
        src[p] = open(p, errors='replace').read().split('\n')

    # index every routine: name -> (path, dummy list, body text)
    routines = {}
    for p, lines in src.items():
        stmts = list(statements(lines))
        for n, (ln, st) in enumerate(stmts):
            m = re.match(r'\s*(?:recursive\s+)?subroutine\s+(\w+)\s*\(', st, re.I)
            if not m:
                continue
            name = m.group(1)
            end = re.compile(r'^\s*end subroutine\s+' + re.escape(name) + r'\s*$',
                             re.I)
            body = []
            for l in lines[ln:]:
                body.append(l)
                if end.match(l):
                    break
            routines[name.lower()] = (p, arglist(st, name), '\n'.join(body))

    bad = 0
    for p, lines in src.items():
        for ln, st in statements(lines):
            m = re.search(r'\bcall\s+(\w+)\s*\(', st, re.I)
            if not m:
                continue
            key = m.group(1).lower()
            if key not in routines:
                continue
            rp, dummies, body = routines[key]
            actuals = arglist(st, m.group(1))
            if not actuals or not dummies or len(actuals) != len(dummies):
                continue
            if any('=' in a for a in actuals):
                continue
            for a, d in zip(actuals, dummies):
                if '_slot' not in a.lower():
                    continue
                # how is this dummy read inside the routine?
                uses = re.findall(r'\b' + re.escape(d) + r'\s*\(([^()]*)\)', body)
                blockwise = [u for u in uses if BLOCK_IX.search('(' + u + ')')]
                slotwise = [u for u in uses if SLOT_IX.search(u)]
                if blockwise and not slotwise:
                    print(f'{p}:{ln+1}: {m.group(1)}: {a.split("(")[0]} -> dummy '
                          f'{d!r} is read block-wise: {d}({blockwise[0].strip()})')
                    bad += 1
    return 1 if bad else 0


def pack_check(paths):
    """Slot arrays unpacked after a subcycle but never packed before it."""
    bad = 0
    for p in paths:
        lines = open(p, errors='replace').read().split('\n')
        starts = [i for i, l in enumerate(lines) if 'do ksub' in l]
        ends = [i for i, l in enumerate(lines)
                if re.match(r'\s+enddo\s+! sub ?cycling', l)]
        for a, b in zip(starts, ends):
            packed = {re.search(r'pack_to_slots\(\s*(\w+)\s*\(', lines[i]).group(1).lower()
                      for i in range(max(0, a - 80), a)
                      if re.search(r'pack_to_slots\(\s*(\w+)\s*\(', lines[i])}
            unpacked = {re.search(r'unpack_from_slots\((\w+)_slot', lines[i]).group(1).lower()
                        for i in range(b, min(b + 45, len(lines)))
                        if re.search(r'unpack_from_slots\((\w+)_slot', lines[i])}
            missing = sorted(unpacked - packed)
            if missing:
                print(f'{p}:{a+1}: subcycle unpacks but never packs: '
                      + ', '.join(missing))
                bad += 1
    return 1 if bad else 0


if __name__ == '__main__':
    argv = sys.argv[1:]
    if argv and argv[0] == '--pack':
        if len(argv) < 2:
            print(__doc__)
            sys.exit(2)
        sys.exit(pack_check(argv[1:]))
    if not argv:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(argv))
