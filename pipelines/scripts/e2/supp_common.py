"""Shared inputs for the motif supplement: the finished set, sequences and motif sets.

Paths come from the environment, as the other E2 steps take them, because two of these
steps are run as a prefix of another and argv cannot reach them.

The finished table is digested PARSED rather than as a file. `finished_tsv.py` writes an
ISO timestamp into the header's first line and the table is regenerated hourly, so a digest
over the whole file changes while the finished set does not, and a pin over it could only
ever fail.
"""

import hashlib
import json
import os

from cellset import excluded

E = os.environ["E2_DEPOSIT"]
MOTIF_DIR = os.environ["E2_MOTIF_DIR"]
OUT = os.environ["E2_SUPP"]
# names.json keys the patterns as the set files spell them, with T for U.
NAMES_JSON = os.environ.get("E2_MOTIF_NAMES", f"{MOTIF_DIR}/names.json")
# E1's auxiliary tree of screened instances. Three members of the corpus have no
# unconstrained front in E2's own screen -- P22301 and P41159 were added to E1 after the
# snapshot E2 screened, and PRT-03 is only there -- and their fronts, which is to say their
# sequences, are here. Read after E2's own directories, so a protein E2 screened itself is
# always taken from E2's copy.
E1_EXTRA = os.environ.get("E2_E1_EXTRA", "")
CORPUS_JSON = os.environ.get("E2_CORPUS", "")

# The set files, in the order config.yaml lists them: the levels add them in that order and
# a set's Roman numeral is its position. Taken from the environment rather than written here,
# because the same list names the rules' inputs and two copies differ the first time one of
# them is edited.
SET_FILES = os.environ["E2_MOTIF_SETS"].split(",")
SETS = dict(zip(["I", "II", "III"], SET_FILES, strict=True))
LEVELS = ["F0", "F1", "F2", "F3"]
NESTS = ["F1", "F2", "F3"]
# IUPAC codes a motif set may use. T is read as U throughout, as the set files state.
IUPAC = {"R": "[AG]", "Y": "[CU]", "N": "[ACGU]", "W": "[AU]", "S": "[CG]"}
CODONS = {
    "A": ["GCU", "GCC", "GCA", "GCG"], "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAU", "AAC"], "D": ["GAU", "GAC"], "C": ["UGU", "UGC"],
    "Q": ["CAA", "CAG"], "E": ["GAA", "GAG"], "G": ["GGU", "GGC", "GGA", "GGG"],
    "H": ["CAU", "CAC"], "I": ["AUU", "AUC", "AUA"], "M": ["AUG"],
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"], "K": ["AAA", "AAG"],
    "F": ["UUU", "UUC"], "P": ["CCU", "CCC", "CCA", "CCG"],
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "T": ["ACU", "ACC", "ACA", "ACG"], "W": ["UGG"], "Y": ["UAU", "UAC"],
    "V": ["GUU", "GUC", "GUA", "GUG"],
}


def motifs():
    """[(set, pattern)] in the cumulative order the levels add them."""
    out = []
    for s, name in SETS.items():
        for line in open(f"{MOTIF_DIR}/{name}.txt"):
            m = line.split("#")[0].strip().upper().replace("T", "U")
            if m:
                out.append((s, m))
    return out


def names():
    """{pattern (U-spelled): display name}, in the main plate's abbreviation."""
    raw = json.load(open(NAMES_JSON))
    return {k.replace("T", "U"): v.replace(", reverse complement", " (rc)")
            for k, v in raw.items()}


def finished():
    """{cell: {nest: round_file}} for every cell the table reports as done, and its digest.

    A cell counts when its status is `complete` (exact front) or `retired` (bounded); an
    `open` cell has no final front to read. The digest is over the parsed rows only.

    Rows that `cellset` excludes are dropped here rather than filtered by each step, so the
    supplement and the composite plate report over one set. A protein excluded at every
    rung leaves with no rung of it remaining, and the steps below, which require every
    required level, then never see it.
    """
    rows = []
    for line in open(f"{E}/FINISHED.tsv"):
        if line.startswith("#") or line.startswith("cell\t") or not line.strip():
            continue
        rows.append(line.rstrip("\n").split("\t"))
    rows = [f for f in rows if not excluded(*f[0].rsplit("_", 1))]
    tsv = {}
    for f in rows:
        cell, nest = f[0].rsplit("_", 1)
        tsv.setdefault(cell, {})[nest] = {"path": f[3], "status": f[2], "band": f[1]}
    # The round file is keyed by its BASENAME. The table spells it as an absolute path, so
    # a digest over the column moves when the deposit is mounted anywhere else, and three
    # pinned caches would change with no number in them changing. The basename is what
    # identifies the round -- cell, level and round index -- which is what this digest is
    # for.
    keyed = sorted((f[0], f[2], os.path.basename(f[3])) for f in rows
                   if f[2] in ("complete", "retired"))
    digest = hashlib.sha256(json.dumps(keyed).encode()).hexdigest()
    return tsv, digest


def sequences():
    """{cell: residue sequence} from each protein's unconstrained front header.

    CR is stripped explicitly and a trailing stop is dropped: the corpus fastas are mixed
    CRLF, so a reader that does not strip carries one `\\r` per sequence line into the
    length, which scales as n / 60.
    """
    seqs, paths = {}, {}
    dirs = [f"{E}/screen/{d}" for d in ("f0", "f0-l", "f0-five", "f0-xl")]
    if E1_EXTRA:
        dirs.append(f"{E1_EXTRA}/f0")
    for p in dirs:
        if not os.path.isdir(p):
            continue
        for fn in sorted(os.listdir(p)):
            if not fn.endswith(".jsonl") or fn[:-6] in seqs:
                continue
            with open(f"{p}/{fn}") as fh:
                h = json.loads(fh.readline())
            seqs[fn[:-6]] = h["protein"]["sequence"].replace("\r", "").rstrip("*")
            paths[fn[:-6]] = f"{p}/{fn}"
    return seqs, paths


def corpus():
    """The proteins the figures report over, or None when no list is given.

    One list, committed, rather than two derivations that agree until they do not: the
    composite plate derives the same set from E1's sweep and asserts it equals this file.
    """
    if not CORPUS_JSON:
        return None
    return set(json.load(open(CORPUS_JSON))["members"])


def points(path):
    """[(cai, mfe, sequence)] for the Pareto optimal designs a front file records."""
    out = []
    for line in open(path):
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("record") == "point":
            out.append((r["cai"], r["mfe"], r.get("sequence")))
    return out


def header(path):
    with open(path) as fh:
        return json.loads(fh.readline())


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def motif_files():
    """The set files the levels forbid, in the order they are added."""
    return [f"{MOTIF_DIR}/{name}.txt" for name in SET_FILES]


def digest_of(*paths):
    """{basename: sha256} -- what a cache records of the code and the inputs it read."""
    return {os.path.basename(p): sha(p) for p in paths}


def write(name, payload):
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/{name}", "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)
    print(f"wrote {OUT}/{name}")
