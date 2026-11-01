"""The band owner's reported set, as recorded in CELL-SET.md — one source of truth.

Exclusion is **per rung**, not per protein. CELL-SET.md excludes only Q9HC38's F4 and retains
its F1, F2 and F3, so a protein-level filter would drop three fronts the owner kept. A rung's
population is therefore its own: every count printed against a rung carries the denominator of
the rows that rung is reported over, and the excluded members are named rather than summarized
as a composition.

Two kinds of exclusion live here and they expire differently. A *cost* judgement is taken
against a particular build and expires when the instrument changes rather than when the data
does. A *corpus* judgement — a sequence that is not a member of the population the figures
report over — holds for as long as that population is what it is, and is listed over every
rung of the protein so that `rungs` returns nothing and the protein leaves the finished set
whole, rather than appearing as one reported over fewer depths.

Re-read CELL-SET.md whenever a count disagrees with the band owner's: a rung whose exclusion
has been reversed stays excluded here until this list is edited, and the symptom is a finished
count one *higher* than the owner's, with the rung simply absent from that protein's `rungs`
— which reads as a protein reported over fewer depths rather than as an error, and so
survives being looked at.

Excluded cells are absent from FINISHED.tsv but their files remain on disk, so a directory
scan picks them up and nothing downstream would notice.
"""

_CORPUS = ("not a human protein; the figures report over the human corpus")

EXCLUDE = {
    ("Q9HC38", "F4"): "excluded per CELL-SET.md (two rounds, share 0.9775)",
    # eGFP is a reporter and RBD a viral domain: the only two entries of the finished set
    # that are not UniProt accessions. Every rung, so each leaves the set whole.
    ("eGFP", "F1"): _CORPUS, ("eGFP", "F2"): _CORPUS, ("eGFP", "F3"): _CORPUS,
    ("eGFP", "F4"): _CORPUS,
    ("RBD", "F1"): _CORPUS, ("RBD", "F2"): _CORPUS, ("RBD", "F3"): _CORPUS,
    ("RBD", "F4"): _CORPUS,
}


def excluded(cell, nest):
    return (cell, nest) in EXCLUDE


def rungs(cell, nests):
    """The rungs of `nests` this cell is reported over."""
    return [n for n in nests if not excluded(cell, n)]
