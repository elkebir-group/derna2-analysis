#!/usr/bin/env python3
"""Which motifs each protein's codon lattice can spell, and where PP sits.

The binding screen has two halves. This is the half that needs no design run:
whether a forbidden motif has any placement the lattice admits at all. It
separates *this protein's optimum happens to avoid the motif* from *no design
could contain it* — a motif with no spellable placement forbids nothing, so a
zero price measured against it is a price against no constraant at all.

A placement is decided residue by reside. A motif occupying positions
[s, s + m) of the coding sequence overlaps some run of codons, and imposes on
each one a slice of itself; the codons are chosen independently, so the motif
is spellable at s exactly when every touched residue has at least one
synonymous codon matching its slice. `spellable` counts the start positions
that pass, and `forced` (from `forced_motifs.py`) is the opposite extreme.

Reads the codon table and the motif files; writes one JSON record per
(protein, motif). Stdlib only, no runs, exact.

Usage: `python3 -I screen_lattice.py <corpus-dir> [more dirs...]`
"""

import csv
import glob
import json
import os
import sys

IUPAC = {
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "U",
    "R": "AG", "Y": "CU", "S": "CG", "W": "AU", "K": "GU", "M": "AC",
    "B": "CGU", "D": "AGU", "H": "ACU", "V": "ACG", "N": "ACGU",
}

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CODONS = os.path.join(REPO, "packages/derna-codon-dp/data/codon_usage_human.csv")


def codon_menus(path):
    """amino acid -> its synonymous codons, as RNA."""
    menus = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            aa = row["amino_acid"].strip()
            if aa in ("*", "Stop"):
                continue
            menus.setdefault(aa, []).append(row["codon"].strip().replace("T", "U"))
    return menus


def motif_sets():
    """set label -> its IUPAC patterns, as RNA, comments stripped."""
    sets = {}
    for path in sorted(glob.glob(os.path.join(HERE, "motifs", "set?-*.txt"))):
        label = os.path.basename(path).split("-", 1)[1].rsplit(".", 1)[0]
        pats = []
        for line in open(path):
            pat = line.split("#")[0].strip().upper()
            if pat:
                pats.append(pat.replace("T", "U"))
        sets[label] = pats
    return sets


def read_fasta(path):
    """The one record's sequence, uppercased, whitespace removed."""
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip().upper())
    return "".join(seq)


def matches(codon, pattern_slice, offset):
    """Does `codon` match `pattern_slice` placed at `offset` within it?"""
    return all(
        codon[offset + i] in IUPAC[ch]
        for i, ch in enumerate(pattern_slice)
    )


def spellable_starts(peptide, pattern, menus):
    """The coding-sequence offsets at which the lattice admits `pattern`.

    Each touched residue is checked against the slice of the pattern that lands
    inside its codon. The choices are independent, so one matching codon per
    touched residue is enough for the whole placement.
    """
    m = len(pattern)
    starts = []
    for s in range(0, 3 * len(peptide) - m + 1):
        ok = True
        pos = s
        while pos < s + m and ok:
            residue = pos // 3
            offset = pos % 3
            take = min(3 - offset, s + m - pos)
            sl = pattern[pos - s:pos - s + take]
            ok = any(matches(c, sl, offset) for c in menus[peptide[residue]])
            pos += take
        if ok:
            starts.append(s)
    return starts


def main():
    menus = codon_menus(CODONS)
    sets = motif_sets()
    known = set(menus)

    for root in sys.argv[1:]:
        for fa in sorted(glob.glob(os.path.join(root, "*.fasta"))):
            acc = os.path.basename(fa).rsplit(".", 1)[0]
            pep = read_fasta(fa)
            if not pep or set(pep) - known:
                print(json.dumps({
                    "protein": acc, "skipped": "non-standard residues",
                    "residues": sorted(set(pep) - known),
                }), flush=True)
                continue
            for label, pats in sets.items():
                for pat in pats:
                    starts = spellable_starts(pep, pat, menus)
                    print(json.dumps({
                        "protein": acc,
                        "aa": len(pep),
                        "set": label,
                        "motif": pat.replace("U", "T"),
                        "spellable": len(starts),
                        "first_start": starts[0] if starts else None,
                    }), flush=True)
            print(json.dumps({
                "protein": acc,
                "aa": len(pep),
                "pp_adjacent": "PP" in pep,
                "pp_count": sum(
                    1 for i in range(len(pep) - 1) if pep[i:i + 2] == "PP"
                ),
                "gc_propensity": round(sum(
                    max(c.count("G") + c.count("C") for c in menus[a]) / 3
                    for a in pep
                ) / len(pep), 4),
            }), flush=True)


if __name__ == "__main__":
    main()
