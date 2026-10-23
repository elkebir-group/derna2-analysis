#!/usr/bin/env python3
"""The binding screen's table: what binds, what is merely spellable, what neither.

Joins the two halves. `spellable` (from `screen_lattice.py`) is what the codon
lattice admits at all; `occurs` is what the designs on the baseline front F0
actually spell, and only those occurrences can move the front. The distinction
is the point of the screen:

  occurs > 0                 the motif is on F0, so forbidding it must move the
                             front — a price measured here is a real price
  occurs = 0, spellable > 0  the lattice admits the motif but the optimum avoids
                             it; forbidding it may still bind off the optimum
  spellable = 0              the motif forbids nothing on this protein, and a
                             zero price against it is a price against no
                             constraint at all

A set is credited only for entries that appear, so a set never inherits binding
from a sibling motif. Reports censored proteins as censored rather than as zero.

The screen runs over more than one corpus and its front records need not sit in
the worktree, so each half is a path. The defaults are the ≤ 500 aa screen's.

Usage: `python3 screen_occurs.py [--fronts DIR] [--lattice PATH] [--out PATH]`,
from this directory. Stdlib only. Unlike `screen_lattice.py` this reads nothing
from the cloned corpus — only the front records derna wrote and the lattice
table — so it runs without `-I` and imports its sibling normally.
"""

import argparse
import csv
import glob
import itertools
import json
import os

from screen_lattice import IUPAC, motif_sets

HERE = os.path.dirname(os.path.abspath(__file__))
SCREEN = os.path.join(HERE, "screen")


def literals(pattern):
    """An IUPAC pattern as the set of literal RNA strings it matches."""
    return {"".join(p) for p in itertools.product(*(IUPAC[c] for c in pattern))}


def front_sequences(path):
    """Every design on the front, and how many of them sit on the hull."""
    seqs, hull = [], 0
    for line in open(path):
        rec = json.loads(line)
        if "sequence" in rec:
            seqs.append(rec["sequence"].replace("T", "U"))
            hull += bool(rec.get("on_hull"))
    return seqs, hull


def count(seqs, lits, m):
    """Occurrences across the front, and how many designs carry at least one."""
    total = designs = 0
    for s in seqs:
        hits = sum(1 for i in range(len(s) - m + 1) if s[i:i + m] in lits)
        total += hits
        designs += hits > 0
    return total, designs


def options(argv=None):
    """Where the join reads its two halves and writes its table."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fronts", default=os.path.join(SCREEN, "f0"),
                    help="per-protein front records, with summary.csv beside them")
    ap.add_argument("--lattice", default=os.path.join(SCREEN, "lattice.jsonl"),
                    help="the lattice half, as screen_lattice.py writes it")
    ap.add_argument("--out", default=os.path.join(SCREEN, "binding.jsonl"),
                    help="the joined table")
    return ap.parse_args(argv)


def censored(fronts):
    """Proteins whose F0 did not complete, by accession.

    A protein can hold both a censored row and a later successful one: the
    screen is re-runnable, and a cell lost to something other than the budget
    — an external signal, a host restart — is recomputed on a later pass. The
    last row for a protein is the one that stands, so a success supersedes an
    earlier censoring rather than sitting beside it. Reading this as a set of
    names would report a completed protein as censored forever.
    """
    path = os.path.join(fronts, "summary.csv")
    if not os.path.exists(path):
        return {}
    last = {}
    for row in csv.DictReader(open(path)):
        last[row["protein"]] = row["exit"]
    return {p: code for p, code in last.items() if code != "0"}


def lattice(path):
    """(protein, motif) -> spellable count, and the per-protein summaries."""
    spell, summary = {}, {}
    for line in open(path):
        rec = json.loads(line)
        if "motif" in rec:
            spell[(rec["protein"], rec["motif"])] = rec["spellable"]
        elif "pp_adjacent" in rec:
            summary[rec["protein"]] = rec
    return spell, summary


def main():
    args = options()
    sets = motif_sets()
    spell, summary = lattice(args.lattice)
    skipped = censored(args.fronts)

    rows = []
    for path in sorted(glob.glob(os.path.join(args.fronts, "*.jsonl"))):
        acc = os.path.basename(path).rsplit(".", 1)[0]
        seqs, hull = front_sequences(path)
        for label, pats in sets.items():
            for pat in pats:
                dna = pat.replace("U", "T")
                occ, car = count(seqs, literals(pat), len(pat))
                rows.append({
                    "protein": acc,
                    "aa": summary.get(acc, {}).get("aa"),
                    "set": label,
                    "motif": dna,
                    "occurs": occ,
                    "designs_carrying": car,
                    "front_designs": len(seqs),
                    "hull_designs": hull,
                    "spellable": spell.get((acc, dna)),
                })

    out = args.out
    with open(out, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
        fh.flush()
        os.fsync(fh.fileno())

    done = sorted({r["protein"] for r in rows})
    print(f"{len(done)} protein(s) with a complete F0; "
          f"{len(skipped)} censored; {len(summary)} in the corpus")
    if skipped:
        print("  censored: " + ", ".join(f"{k} (exit {v})" for k, v in skipped.items()))

    print("\nper set, over the proteins whose F0 completed:")
    for label in sets:
        sub = [r for r in rows if r["set"] == label]
        on_front = {r["protein"] for r in sub if r["occurs"]}
        latent = {r["protein"] for r in sub if not r["occurs"] and r["spellable"]}
        print(f"  {label:14s} occurs on {len(on_front):2d}/{len(done)}  "
              f"spellable-but-absent on {len(latent):2d}")

    print("\nproteins where every set occurs on F0 (the selection pool):")
    pool = [p for p in done
            if all(any(r["occurs"] for r in rows
                       if r["protein"] == p and r["set"] == s) for s in sets)]
    for p in pool:
        s = summary.get(p, {})
        print(f"  {p}  {s.get('aa')} aa  GC {s.get('gc_propensity')}  "
              f"PP {s.get('pp_adjacent')}")
    if not pool:
        print("  (none yet — the corpus is still running)")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
