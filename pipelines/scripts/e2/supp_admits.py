"""Which motifs each protein can spell at all, over the codon lattice.

    E2_DEPOSIT=... E2_MOTIF_DIR=... E2_SUPP=... python supp_admits.py

A protein ADMITS a motif when at least one synonymous coding sequence of it contains the
motif. The check slides a window of as many residues as the motif can span, enumerates that
window's synonymous codon strings and looks for the pattern, so it is exact rather than a
sample.

A motif no window admits cannot appear in any design of that protein. Forbidding it there
removes nothing: the feasible set is unchanged and the front is identical, so a zero cost
measured against it is the absence of a constraint and not a price paid.
"""

import itertools
import re
import sys

import supp_common as common

MOTIFS = common.motifs()
RX = {m: re.compile("".join(common.IUPAC.get(c, c) for c in m)) for _, m in MOTIFS}


def admits(protein, motif):
    """The first residue window whose synonymous codons spell the motif, or None.

    A motif of length L starting at any phase spans at most ceil((L + 2) / 3) residues, so
    a window of that many residues is enough for the scan to be exact.
    """
    w = -(-(len(motif) + 2) // 3)
    rx = RX[motif]
    for i in range(len(protein) - w + 1):
        win = protein[i:i + w]
        if any(a not in common.CODONS for a in win):
            continue
        for combo in itertools.product(*(common.CODONS[a] for a in win)):
            if rx.search("".join(combo)):
                return {"pos": i, "peptide": win, "codons": "".join(combo)}
    return None


def main():
    seqs, _ = common.sequences()
    # Panel (d) reports over the corpus the rest of the paper does, not over whatever
    # unconstrained fronts happen to be on disk: that scan found 103, twenty-one of them
    # outside the corpus, and a count over it would put a denominator in front of a reader
    # that appears nowhere else in the paper.
    pop = common.corpus()
    if pop is not None:
        missing = sorted(pop - set(seqs))
        assert not missing, f"no sequence for {missing}"
        seqs = {c: v for c, v in seqs.items() if c in pop}
    rows = []
    for cell, seq in sorted(seqs.items()):
        hits = {m: admits(seq, m) for _, m in MOTIFS}
        rows.append({"cell": cell, "aa": len(seq),
                     "admits": {m: h is not None for m, h in hits.items()},
                     "where": {m: h for m, h in hits.items() if h}})
        print(cell, sum(rows[-1]["admits"].values()), "of", len(MOTIFS), flush=True)
    # No digest of the finished set here. This file's population is the proteins with an
    # unconstrained front, which is not the set of closed cells: digesting the latter moved
    # the pin on every grid arrival while every number in the file stayed as it was. The
    # members are named, because a count of 103 cannot report a one-in-one-out swap.
    common.write("admits.json", {
        "motifs": MOTIFS,
        "rows": rows,
        "provenance": {
            "n_proteins": len(rows),
            "members": sorted(seqs),
            "inputs": common.digest_of(*common.motif_files()),
            "scripts": common.digest_of(__file__, common.__file__),
        },
    })


if __name__ == "__main__":
    sys.exit(main())
