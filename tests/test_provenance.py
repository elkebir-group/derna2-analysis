import hashlib

import pytest

from derna2_analysis.provenance import verify


def write_manifest(tmp_path, files):
    lines = [f"{hashlib.sha256(data).hexdigest()}  {name}" for name, data in files.items()]
    (tmp_path / "pins.sha256").write_text("\n".join(lines) + "\n")
    return tmp_path / "pins.sha256"


def test_verify_passes_on_the_pinned_files(tmp_path):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "tool").write_bytes(b"build 1")
    manifest = write_manifest(tmp_path, {"tool": b"build 1"})
    assert verify(manifest, tmp_path / "bin") == ["tool"]


def test_verify_refuses_a_rebuilt_binary(tmp_path):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "tool").write_bytes(b"build 2")
    manifest = write_manifest(tmp_path, {"tool": b"build 1"})
    with pytest.raises(ValueError, match="tool"):
        verify(manifest, tmp_path / "bin")


def test_verify_refuses_a_missing_file(tmp_path):
    (tmp_path / "bin").mkdir()
    manifest = write_manifest(tmp_path, {"tool": b"build 1"})
    with pytest.raises(FileNotFoundError):
        verify(manifest, tmp_path / "bin")


def test_verify_refuses_a_manifest_that_names_nothing(tmp_path):
    # A freeze that has not happened yet must not stamp a directory as pinned.
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "tool").write_bytes(b"any build")
    manifest = tmp_path / "pins.sha256"
    manifest.write_text("# not frozen yet\n")
    with pytest.raises(ValueError, match="names no files"):
        verify(manifest, tmp_path / "bin")
