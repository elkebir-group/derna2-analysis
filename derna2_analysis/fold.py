"""The free fold of one given sequence.

A deposited front records each design's energy and not always its structure, and a
design another tool emitted carries whatever that tool reported. Both are re-folded
here, by the same Turner 2004 implementation the designs were made under: `mfe_fold`
reads one sequence per line and writes `<structure>\t<energy>`, the energy as an
integer hundredth of a kcal/mol.
"""

import subprocess

#: `mfe_fold` writes an integer and the rows everywhere else are kcal/mol. The factor is
#: 100 whatever the field is called at the Rust end: a 234-nt design the deposit records
#: at -148.7 kcal/mol comes back as -14870.
UNITS_PER_KCAL = 100.0


def free_fold(mfe_fold, sequence):
    """One sequence's own minimum-free-energy structure and energy in kcal/mol.

    Free of any constraint: what the molecule does on its own, which is what a
    design's claim about itself is checked against.
    """
    return free_folds(mfe_fold, [sequence])[0]


def free_folds(mfe_fold, sequences):
    """`free_fold` over several sequences in one call, in the order given."""
    sequences = list(sequences)
    if not sequences:
        return []
    for s in sequences:
        if not s or not s.isalpha():
            raise ValueError(f"not a sequence: {s!r}")
    got = subprocess.run(
        [str(mfe_fold)],
        input="\n".join(sequences) + "\n",
        capture_output=True,
        text=True,
        check=True,
    )
    lines = [ln for ln in got.stdout.splitlines() if ln.strip()]
    if len(lines) != len(sequences):
        raise ValueError(f"{mfe_fold}: {len(lines)} fold(s) for {len(sequences)} sequence(s)")
    out = []
    for sequence, line in zip(sequences, lines, strict=True):
        structure, energy = line.split("\t")
        if len(structure) != len(sequence):
            raise ValueError(
                f"{mfe_fold}: a structure of {len(structure)} for a sequence of {len(sequence)}"
            )
        out.append({"structure": structure, "mfe_kcal_mol": int(energy) / UNITS_PER_KCAL})
    return out
