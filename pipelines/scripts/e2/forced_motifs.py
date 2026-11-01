"""Which motif entries the genetic code forces, recomputed and checked against the record.

    python forced_motifs.py DRIVER MOTIF_DIR CODON_TABLE RECORDED OUT

Runs the vendored drivers/e2/forced_motifs.py unchanged, pointed at this repository's motif
sets and codon table rather than at the derna2 checkout it was written in, and fails unless
every per-entry record equals the one the plan reports (`e3-probe/forced.jsonl`). It needs no
binary: the answer follows from the codon table alone.

The record also holds four rows for the union of every set (`joint_len`), which the vendored
script has no mode to produce; they are counted and left unchecked.
"""

import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path


def load(driver):
    spec = importlib.util.spec_from_file_location("forced_motifs", driver)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    driver, motif_dir, table, recorded, out = sys.argv[1:]
    mod = load(driver)
    # The driver finds the sets at HERE/motifs and the table under its repository root.
    mod.HERE = str(Path(motif_dir).resolve().parent)
    read_menus = mod.codon_menus
    mod.codon_menus = lambda _path: read_menus(table)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mod.main()
    got = [json.loads(line) for line in buf.getvalue().splitlines() if line.strip()]
    rows = [json.loads(line) for line in Path(recorded).read_text().splitlines() if line.strip()]
    want = [r for r in rows if "set" in r]
    unchecked = len(rows) - len(want)
    Path(out).write_text(buf.getvalue())
    if got != want:
        diff = [(g, w) for g, w in zip(got, want, strict=False) if g != w]
        sys.exit(f"{len(got)} records against {len(want)} recorded; first difference {diff[:1]}")
    forced = [r["motif"] for r in got if r["forced"]]
    print(
        f"{len(got)} entries match the record; forced: {forced or 'none'}; "
        f"{unchecked} union rows not reproducible from the driver, unchecked"
    )


if __name__ == "__main__":
    main()
