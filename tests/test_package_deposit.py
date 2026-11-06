"""Packing a deposit: the archive holds the bytes, only the manifest's files, and is checked.

The first test is a regression. A deposit is reached through symlinks -- `data/results/<exp>`
points at wherever the tree is mounted, and `e1` and `fig1` are directories of per-file
links -- and `tar.add` stores a symlink as a symlink unless told otherwise. Built that way,
`e1.tar` carried 0 bytes of its 557 MB `sweep.jsonl` and `fig1.tar` lost 4 of its 27 files,
and both reported success: the script verified the source tree, which was intact, and
nothing read the archive back. These tests read the archive back.
"""

import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "pipelines" / "scripts" / "package_deposit.py"


def build(tmp_path, files, link=None):
    """A deposit of `files` {relpath: bytes}, its manifest, and optionally one entry that is
    a symlink to a file living outside the tree -- which is how every real deposit is read."""
    from derna2_analysis.provenance import sha256

    store, root = tmp_path / "store", tmp_path / "deposit"
    store.mkdir(), root.mkdir()
    for rel, data in files.items():
        (store / rel).parent.mkdir(parents=True, exist_ok=True)
        (store / rel).write_bytes(data)
        target = root / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if link == rel:
            target.symlink_to(store / rel)
        else:
            target.write_bytes(data)
    manifest = tmp_path / "m.sha256"
    manifest.write_text("".join(f"{sha256(root / rel)}  {rel}\n" for rel in sorted(files)))
    return root, manifest


def pack(tmp_path, root, manifest, prefix="exp"):
    out = tmp_path / "out.tar"
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(manifest), str(root), str(out), prefix],
        capture_output=True, text=True, cwd=ROOT,
    ), out


def test_a_symlinked_file_is_packed_as_its_bytes(tmp_path):
    files = {"sweep.jsonl": b'{"a":1}\n' * 100, "hull/x.json": b"{}"}
    root, manifest = build(tmp_path, files, link="sweep.jsonl")
    assert (root / "sweep.jsonl").is_symlink(), "the fixture must exercise a symlink"
    got, out = pack(tmp_path, root, manifest)
    assert got.returncode == 0, got.stderr
    with tarfile.open(out) as tar:
        members = {m.name: m for m in tar if not m.isdir()}
        assert not any(m.issym() for m in members.values()), "a symlink reached the archive"
        assert members["exp/sweep.jsonl"].size == len(files["sweep.jsonl"])
        assert tar.extractfile("exp/sweep.jsonl").read() == files["sweep.jsonl"]


def test_only_the_manifests_files_are_packed(tmp_path):
    root, manifest = build(tmp_path, {"kept.jsonl": b"1\n"})
    (root / "run.log").write_text("/home/someone/ran this\n")   # beside it, not in the manifest
    got, out = pack(tmp_path, root, manifest)
    assert got.returncode == 0, got.stderr
    with tarfile.open(out) as tar:
        assert sorted(m.name for m in tar if not m.isdir()) == ["exp/kept.jsonl"]


def test_the_archive_check_rejects_a_member_that_does_not_match(tmp_path):
    # The mirror: an archive check that cannot fail is not a check. Pack, then corrupt one
    # member and run the checker over it directly.
    root, manifest = build(tmp_path, {"a.jsonl": b"one\n", "b.jsonl": b"two\n"})
    got, out = pack(tmp_path, root, manifest)
    assert got.returncode == 0, got.stderr
    sys.path.insert(0, str(SCRIPT.parent))
    import package_deposit

    bad = tmp_path / "bad.tar"
    with tarfile.open(out) as src, tarfile.open(bad, "w") as dst:
        for m in src:
            if m.isfile() and m.name.endswith("a.jsonl"):
                import io

                data = src.extractfile(m).read() + b"x"
                m.size = len(data)
                dst.addfile(m, io.BytesIO(data))
            elif m.isfile():
                dst.addfile(m, src.extractfile(m))
    try:
        package_deposit.check_archive(str(bad), str(manifest), "exp")
    except SystemExit as exc:
        assert "a.jsonl" in str(exc), str(exc)
    else:
        raise AssertionError("the checker passed a corrupted archive")
