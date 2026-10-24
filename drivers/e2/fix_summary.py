"""Validate the grid's summary rows, and repair the ones a core mask split.

The driver writes a cell's core mask into a CSV field unquoted. At one thread the
mask is a single number and the row is well formed; at four it is `9,11,13,15`,
which occupies four fields instead of one and shifts every later column left by
three. The header has 15 fields and such a row has 18.

What makes this worth a tool rather than an eyeball is that the mis-read values
are *plausible*. On the first threaded row a column-indexed read returns
`exit=11`, `wall=13`, `rss=15`, `points=15` — a failed cell that ran for thirteen
seconds, which is a readable and entirely wrong story. Nothing about it looks like
corruption, and `awk -F,` reports no error.

The repair is only applied where it is unambiguous: the row must have exactly
`15 + k` fields, and fields 7 through `7 + k` must all be integers that are valid
core numbers on this host. Then they are one mask and the row is rewritten with it
quoted. Anything else is reported and left alone, because a tool that guesses at a
malformed results row is worse than one that refuses.

The driver cannot be fixed while it runs — bash is positioned inside its read
loop and editing the file would shift the offsets underneath it — so this runs
after the fact, and is idempotent so it can run as often as wanted.
"""

import argparse
import csv
import io
import os
import shutil
import sys

CORES = set(range(28))          # ada: 28 cores, 0-27


def header_width(path):
    with open(path) as fh:
        return len(next(csv.reader(fh)))


def csv_width(line):
    """Field count as a CSV PARSER sees it, so a quoted mask counts as one.

    Detecting with `line.split(",")` is not idempotent: a row this tool has
    already repaired to `"9,11,13,15"` still splits into four, so the tool
    refuses its own output. That went unnoticed because the idempotence check
    read only the last two lines of the report, below where the refusal printed.
    """
    try:
        return len(next(csv.reader(io.StringIO(line))))
    except Exception:
        return -1


def repair_row(fields, width, core_at=6):
    """Return (repaired_fields, mask) or (None, reason) if not unambiguous."""
    surplus = len(fields) - width
    if surplus <= 0:
        return None, "not wide"
    span = fields[core_at:core_at + surplus + 1]
    if not all(f.isdigit() and int(f) in CORES for f in span):
        return None, "surplus fields are not core numbers: %s" % ",".join(span)
    mask = ",".join(span)
    out = fields[:core_at] + ['"%s"' % mask] + fields[core_at + surplus + 1:]
    if len(out) != width:
        return None, "repair would give %d fields, not %d" % (len(out), width)
    return out, mask


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path", nargs="?",
                    default="/scratch/xinyugu4/paper-e3/grid/summary.csv")
    ap.add_argument("--apply", action="store_true",
                    help="rewrite the file; without it, only report")
    args = ap.parse_args()

    width = header_width(args.path)
    with open(args.path) as fh:
        lines = fh.read().splitlines()

    out, fixed, refused = [], [], []
    for i, line in enumerate(lines):
        if i == 0 or not line.strip():
            out.append(line)
            continue
        if csv_width(line) == width:
            out.append(line)          # already well formed, quoted mask included
            continue
        fields = line.split(",")
        repaired, info = repair_row(fields, width)
        if repaired is None:
            refused.append((i + 1, len(fields), info, line))
            out.append(line)
        else:
            fixed.append((i + 1, len(fields), info))
            out.append(",".join(repaired))

    print("%s: header %d fields, %d data row(s)" % (args.path, width, len(lines) - 1))
    print("  well formed : %d" % sum(1 for l in lines[1:]
                                     if l.strip() and csv_width(l) == width))
    print("  repairable  : %d" % len(fixed))
    for n, got, mask in fixed:
        print("     line %d had %d fields; core mask %s" % (n, got, mask))
    print("  refused     : %d" % len(refused))
    for n, got, why, _ in refused:
        print("     line %d had %d fields: %s" % (n, got, why))

    if not fixed:
        print("\nnothing to do")
        return 1 if refused else 0

    if not args.apply:
        print("\nreport only; pass --apply to rewrite")
        return 0

    backup = args.path + ".unrepaired"
    if not os.path.exists(backup):
        shutil.copy2(args.path, backup)
        print("\nkept the original at %s" % backup)
    tmp = args.path + ".partial"
    with open(tmp, "w") as fh:
        fh.write("\n".join(out) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.rename(tmp, args.path)
    print("rewrote %s with %d row(s) repaired" % (args.path, len(fixed)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
