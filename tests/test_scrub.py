"""The host-path scrub: it reports before it writes, writes only what it was given, and
fails on a path no substitution covers rather than publishing it."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRUB = ROOT / "pipelines" / "scripts" / "scrub_host_paths.py"
WORKTREE = "/home/someone/checkout/.hidden/worktrees/a-name"


def deposit(tmp_path, text):
    (tmp_path / "raw").mkdir()
    out = tmp_path / "raw" / "run.log"
    out.write_text(text)
    out.chmod(0o444)
    return out


def scrub(tmp_path, *args):
    return subprocess.run(
        [sys.executable, str(SCRUB), str(tmp_path), *args], capture_output=True, text=True
    )


def test_it_reports_without_writing(tmp_path):
    log = deposit(tmp_path, f"read {WORKTREE}/fixtures/x.fasta\n")
    got = scrub(tmp_path, f"{WORKTREE}/=derna2/")
    assert got.returncode == 0, got.stderr
    assert "raw/run.log: 1x" in got.stdout
    assert "would rewrite" in got.stdout
    assert WORKTREE in log.read_text()


def test_apply_rewrites_the_file_and_leaves_it_read_only(tmp_path):
    log = deposit(tmp_path, f"read {WORKTREE}/fixtures/x.fasta\n")
    got = scrub(tmp_path, f"{WORKTREE}/=derna2/", "--apply")
    assert got.returncode == 0, got.stderr
    assert log.read_text() == "read derna2/fixtures/x.fasta\n"
    assert log.stat().st_mode & 0o222 == 0


def test_it_fails_on_a_host_path_no_substitution_covers(tmp_path):
    deposit(tmp_path, f"read {WORKTREE}/x.fasta\nwrote /scratch/elsewhere/out.jsonl\n")
    got = scrub(tmp_path, f"{WORKTREE}/=derna2/")
    assert got.returncode == 1
    assert "/scratch/elsewhere/out.jsonl" in got.stderr


def test_it_refuses_an_argument_that_is_not_a_substitution(tmp_path):
    deposit(tmp_path, "nothing to do\n")
    assert scrub(tmp_path, "derna2/").returncode != 0
