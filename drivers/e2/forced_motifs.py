#!/usr/bin/env python3
"""Which motifs the genetic code forces, and on which peptides.

A forbidden motif is *forced* by a peptide when every synonymous encoding of
that peptide contains the motif: no design avoids it, so any protein carrying
the peptide is infeasible under a set containing that entry. The search is over
peptides of length 1..L, which is sufficient because a motif of length m can
span at most ceil((m + 4) / 3) residues.

Reads the codon table and the motif files; writes one JSONL record per entry.
Stdlib only.
"""

import csv
import glob
import itertools
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


def expand(pattern):
    """An IUPAC pattern as the set of literal RNA strings it matches."""
    return {"".join(p) for p in itertools.product(*(IUPAC[c] for c in pattern.upper()))}


def contains(text, literals, m):
    return any(text[i:i + m] in literals for i in range(len(text) - m + 1))


def forced_peptides(pattern, menus, max_len):
    """The shortest peptides every encoding of which contains `pattern`.

    Returns a list of peptides, and stops at the first length that yields any:
    a longer peptide containing a forced one is forced for the same reason, so
    only the minimal ones carry information.
    """
    literals = expand(pattern)
    m = len(pattern)
    aas = sorted(menus)
    for length in range(1, max_len + 1):
        hits = []
        for pep in itertools.product(aas, repeat=length):
            if all(
                contains("".join(cs), literals, m)
                for cs in itertools.product(*(menus[a] for a in pep))
            ):
                hits.append("".join(pep))
        if hits:
            return length, hits
    return None, []


def main():
    menus = codon_menus(os.path.join(REPO, "packages/derna-codon-dp/data/codon_usage_human.csv"))
    entries = []
    for path in sorted(glob.glob(os.path.join(HERE, "motifs", "set?-*.txt"))):
        label = os.path.basename(path).split("-", 1)[1].rsplit(".", 1)[0]
        for line in open(path):
            pat = line.split("#")[0].strip()
            if pat:
                entries.append((label, pat))
    # Controls: the two thresholds the plan reports as infeasible.
    entries += [("control", "GGGG"), ("control", "CCNCCN")]

    for label, pat in entries:
        # A motif of length m spans at most ceil((m + 4) / 3) residues.
        max_len = -(-(len(pat) + 4) // 3)
        length, peps = forced_peptides(pat, menus, max_len)
        rec = {
            "set": label,
            "motif": pat,
            "forced": bool(peps),
            "peptide_len": length,
            "peptides": peps[:20],
            "peptide_count": len(peps),
            "searched_to": max_len,
        }
        print(json.dumps(rec))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
