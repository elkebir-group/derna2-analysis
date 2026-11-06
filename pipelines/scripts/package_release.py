"""Pack every deposit into the one archive a data repository takes.

    python package_release.py OUT_TAR [--prefix NAME]

`package_deposit.py` writes one archive per experiment, which is what a rule asks for. A
published record wants the opposite: one file, one identifier, one download, as the sibling
paper's deposit does (`h2opt.tar`, unpacked and pointed at by a config key). This writes
that, laid out so a reader sets one path:

    <prefix>/README.txt            what the archive is and how to use it
    <prefix>/pins/results-*.sha256 every manifest, so the archive checks itself
    <prefix>/results/<exp>/...     the deposits

and then `results_root: /path/to/<prefix>/results` in `pipelines/config/host.yaml`.

`e1_fig4` is deliberately not a directory here. All 171 files it names are byte-identical to
`e1`'s at the same relative paths, so its manifest verifies against `results/e1` and shipping
the tree twice would duplicate 1.20 GiB for nothing. Its manifest travels, because it is what
says which 171 of E1's 368 files Figure 4 was drawn from.

Every member is read back out of the finished archive and hashed against the manifest it came
from, for the reason `package_deposit.py` does it: an archive that was never read back is an
artifact nobody has checked, and it is the only thing a reader gets.
"""

import hashlib
import io
import sys
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PINS = ROOT / "pins"
RESULTS = ROOT / "data" / "results"

# The deposits that ship as trees. e1_fig4 is a subset of e1; see the module docstring.
TREES = ["e1", "e1_extra", "e2", "e3", "fig1", "fig1_fronts"]
SUBSETS = {"e1_fig4": "e1"}


def manifest(exp):
    out = []
    for line in (PINS / f"results-{exp}.sha256").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            digest, rel = line.split(maxsplit=1)
            out.append((digest, rel.strip()))
    return out


def readme(prefix):
    lines = [
        "DERNA2 analysis — deposited results",
        "",
        "The recorded results of the paper's experiments. The code that reads them is at",
        "https://github.com/elkebir-group/derna2-analysis, and the tool that produced them at",
        "https://github.com/elkebir-group/derna2.",
        "",
        "To use: unpack this archive and name its `results` directory in the analysis",
        "repository's host configuration (`pipelines/config/host.yaml`):",
        "",
        f"    results_root: /path/to/{prefix}/results",
        "",
        "Then `pipelines/run.sh` rebuilds every figure and table from these files and checks",
        "each one against the digests in pins/.",
        "",
        "Contents",
        "--------",
    ]
    total = 0
    for exp in TREES:
        entries = manifest(exp)
        size = sum((RESULTS / exp / rel).stat().st_size for _, rel in entries)
        total += size
        lines.append(f"  results/{exp:<12} {len(entries):>5} files  {size / 2**30:8.3f} GiB")
    lines += [
        f"  {'':<21}{sum(len(manifest(e)) for e in TREES):>5} files  {total / 2**30:8.3f} GiB",
        "",
        "Checksums",
        "---------",
        "pins/results-<exp>.sha256 names every file of that deposit, in sha256sum format with",
        "paths relative to results/<exp>. From the unpacked archive:",
        "",
        "    cd results/e2 && sha256sum -c ../../pins/results-e2.sha256",
        "",
        "Every file in this archive was verified against its manifest when the archive was",
        "written, and again by reading it back out of the finished archive.",
        "",
        "pins/results-e1_fig4.sha256 has no directory of its own on purpose. All 171 files it",
        "names are byte-identical to results/e1's at the same relative paths — Figure 4 was",
        "drawn from that subset of E1's 368 — so it verifies against results/e1 and the tree",
        "is not shipped twice:",
        "",
        "    cd results/e1 && sha256sum -c ../../pins/results-e1_fig4.sha256",
        "",
        "Provenance",
        "----------",
        "Each deposit's pins/<exp>-PROVENANCE.md in the analysis repository records what the",
        "deposit is, which run produced it, what was scrubbed from it and what the figures",
        "built from it were checked against. None of these files names a host or a person.",
        "",
        f"Written {time.strftime('%Y-%m-%d')}.",
        "",
    ]
    return "\n".join(lines).encode()


def add_bytes(tar, arcname, data):
    info = tarfile.TarInfo(arcname)
    info.size = len(data)
    info.mtime = int(time.time())
    info.mode = 0o444
    tar.addfile(info, io.BytesIO(data))


def main():
    args = sys.argv[1:]
    prefix = "derna2-analysis-results"
    if "--prefix" in args:
        i = args.index("--prefix")
        prefix = args[i + 1]
        args = args[:i] + args[i + 2:]
    (out,) = args

    want = {}
    with tarfile.open(out, "w", dereference=True) as tar:
        add_bytes(tar, f"{prefix}/README.txt", readme(prefix))
        for exp in TREES + list(SUBSETS):
            src = PINS / f"results-{exp}.sha256"
            add_bytes(tar, f"{prefix}/pins/{src.name}", src.read_bytes())
        for exp in TREES:
            for digest, rel in manifest(exp):
                arc = f"{prefix}/results/{exp}/{rel}"
                tar.add(RESULTS / exp / rel, arcname=arc, recursive=False)
                want[arc] = digest

    bad, seen = [], set()
    with tarfile.open(out, "r") as tar:
        for member in tar:
            if member.isdir() or member.name not in want:
                continue
            if not member.isfile():
                bad.append(f"{member.name}: stored as {'a symlink' if member.issym() else 'not a file'}")
                continue
            handle = tar.extractfile(member)
            digest = hashlib.sha256()
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
            if digest.hexdigest() != want[member.name]:
                bad.append(f"{member.name}: sha256 {digest.hexdigest()[:12]}, "
                           f"manifest expects {want[member.name][:12]}")
            seen.add(member.name)
    missing = sorted(set(want) - seen)
    if missing:
        bad.append(f"{len(missing)} manifest file(s) absent: {missing[:3]}")
    if bad:
        sys.exit(f"{out}: does not match the manifests:\n  " + "\n  ".join(bad[:10]))
    print(f"{out}: {len(want)} files from {len(TREES)} deposits, all verified inside the archive")
    print(f"  plus README.txt and {len(TREES) + len(SUBSETS)} manifests")


if __name__ == "__main__":
    main()
