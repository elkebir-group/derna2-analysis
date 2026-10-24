#!/usr/bin/env python3
"""Convert the canonical human codon table into DERNA v1's `-c` format.

v1 reads the table with `fill_codon`: lines come in pairs, a header line whose
first field is empty and whose remaining 3-character fields name the codons,
then a value line whose first field is the amino acid letter and whose
remaining fields are that codon's usage, positionally matched to the header.

The values need no normalization. v1 computes
`codon_cai[aa][codon] = log(usage / max_usage_for_that_aa)`, so any positive
scaling of one amino acid's row gives identical weights — frequencies and raw
counts both work.

Stop codons are dropped: v1 indexes 20 amino acids and appends the stop itself.

Usage: derna1_table.py CANONICAL_CSV OUT_CSV
"""
import sys
from collections import OrderedDict

STOP = "*"


def read_canonical(path):
    """`CODON,AA,freq` triples -> {aa: OrderedDict(codon -> freq)}."""
    table = OrderedDict()
    with open(path, encoding="utf-8-sig") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3:
                continue
            codon, aa, freq = parts[0], parts[1], parts[2]
            if len(codon) != 3 or aa == STOP:
                continue
            try:
                value = float(freq)
            except ValueError:
                continue
            table.setdefault(aa, OrderedDict())[codon.replace("T", "U")] = value
    return table


def write_derna1(table, path):
    lines = []
    for aa, codons in table.items():
        lines.append("," + ",".join(codons))
        lines.append(aa + "," + ",".join(f"{v:g}" for v in codons.values()))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(src, dst):
    table = read_canonical(src)
    if len(table) != 20:
        raise SystemExit(f"expected 20 amino acids, parsed {len(table)}: "
                         f"{sorted(table)}")
    total = sum(len(c) for c in table.values())
    if total != 61:
        raise SystemExit(f"expected 61 sense codons, parsed {total}")
    for aa, codons in table.items():
        if max(codons.values()) <= 0:
            raise SystemExit(f"{aa} has no positive usage: {dict(codons)}")
    write_derna1(table, dst)
    print(f"wrote {dst}: 20 amino acids, {total} sense codons")
    for aa in list(table)[:3]:
        shown = " ".join(f"{c}={v:g}" for c, v in table[aa].items())
        print(f"  {aa}  {shown}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
