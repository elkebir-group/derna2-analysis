"""Pack the deposit of one experiment: exactly the files its manifest names, nothing else.

    python package_deposit.py MANIFEST ROOT OUT_TAR PREFIX

The result directories also hold run logs and working scripts, which carry the host's paths
and are not part of the record; packing from the manifest leaves them out. Every file is
checked against its digest first.

A deposit is usually reached through symlinks -- `data/results/<exp>` points at wherever the
tree is mounted, and some deposits are a directory of per-file links. The archive must hold
the *bytes*, so it is opened with `dereference=True`: a symlink stored as a symlink arrives
at the reader as a dangling link to a path on a machine they do not have. That is not
hypothetical. Built without it, `e1.tar` stored `sweep.jsonl` as a link and carried 0 bytes
of the 557 MB it is meant to publish, and `fig1.tar` lost 4 of its 27 files, while both
reported success -- because `verify` above checks the source tree, which is intact, and
nothing checked the archive.

So this also checks what it wrote: every member is read back out of the finished archive and
hashed against the same manifest. It roughly doubles the work and is worth it, because the
failure it catches is silent, and the artifact it catches it in is the one thing a reader
gets.
"""

import sys
import tarfile
from pathlib import Path

from derna2_analysis.provenance import read_manifest, sha256, verify


def check_archive(out: str, manifest: str, prefix: str) -> int:
    """Hash every member of the finished archive against the manifest it was built from."""
    want = {rel: digest for digest, rel in read_manifest(manifest)}
    seen, bad = set(), []
    with tarfile.open(out, "r") as tar:
        for member in tar:
            if member.isdir():
                continue
            rel = member.name[len(prefix) + 1:] if member.name.startswith(prefix + "/") else member.name
            if not member.isfile():
                bad.append(f"{rel}: stored as {'a symlink' if member.issym() else 'not a file'}")
                continue
            handle = tar.extractfile(member)
            if handle is None:
                bad.append(f"{rel}: no contents")
                continue
            import hashlib

            digest = hashlib.sha256()
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
            if digest.hexdigest() != want.get(rel):
                bad.append(f"{rel}: sha256 {digest.hexdigest()[:12]}, manifest expects "
                           f"{str(want.get(rel))[:12]}")
            seen.add(rel)
    missing = sorted(set(want) - seen)
    if missing:
        bad.append(f"{len(missing)} file(s) of the manifest are not in the archive: {missing[:3]}")
    if bad:
        sys.exit(f"{out}: the archive does not match {manifest}:\n  " + "\n  ".join(bad[:10]))
    return len(seen)


def main():
    manifest, root, out, prefix = sys.argv[1:]
    files = verify(manifest, root)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    # dereference: store the bytes a link points at, never the link. See the module docstring.
    with tarfile.open(out, "w", dereference=True) as tar:
        for rel in files:
            tar.add(Path(root) / rel, arcname=f"{prefix}/{rel}", recursive=False)
    checked = check_archive(out, manifest, prefix)
    print(f"{out}: {len(files)} files from {manifest}, {checked} verified inside the archive")


if __name__ == "__main__":
    main()
