"""Checksums of binaries and deposited results against the manifests under pins/."""

import hashlib
from pathlib import Path


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def read_manifest(path):
    """`sha256sum` format: '<hex>  <relative path>' per line."""
    entries = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            digest, rel = line.split(maxsplit=1)
            entries.append((digest, rel.strip()))
    return entries


def verify(manifest, root):
    """Every file the manifest names, under root, with the digest it names.

    Returns the checked relative paths; raises on the first missing or changed file, since
    a figure drawn from a different build or a different deposit is a different result.
    """
    root = Path(root)
    checked = []
    for digest, rel in read_manifest(manifest):
        path = root / rel
        if not path.is_file():
            raise FileNotFoundError(f"{path}: named in {manifest}, missing")
        got = sha256(path)
        if got != digest:
            raise ValueError(f"{path}: sha256 {got[:12]}, {manifest} expects {digest[:12]}")
        checked.append(rel)
    return checked
