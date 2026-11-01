"""Checksums of binaries and deposited results against the manifests under pins/."""

import hashlib
from pathlib import Path


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def fnv1a64(path, chunk=1 << 20):
    """FNV-1a over the whole file, which is the `binary_fnv64` a run stamps on its rows.

    The stamp names the binary that wrote a row, so rows can be tied to a pinned file even
    where the source commit was never recorded.
    """
    h = 0xCBF29CE484222325
    with open(path, "rb") as f:
        while block := f.read(chunk):
            for byte in block:
                h = ((h ^ byte) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return f"{h:016x}"


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
    entries = read_manifest(manifest)
    if not entries:
        # A manifest with no entries is a freeze that has not happened yet; passing it would
        # stamp every file under root as pinned.
        raise ValueError(f"{manifest}: names no files")
    root = Path(root)
    checked = []
    for digest, rel in entries:
        path = root / rel
        if not path.is_file():
            raise FileNotFoundError(f"{path}: named in {manifest}, missing")
        got = sha256(path)
        if got != digest:
            raise ValueError(f"{path}: sha256 {got[:12]}, {manifest} expects {digest[:12]}")
        checked.append(rel)
    return checked
