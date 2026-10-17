"""Tracked files carry no paths into a temporary or per-run working directory: such a path
names one machine's one run, means nothing to a reader, and breaks on every other host."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = ["/tmp/", "scratchpad"]


def test_no_tracked_file_points_into_a_temporary_directory():
    tracked = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    hits = []
    for rel in tracked:
        path = ROOT / rel
        if rel == "tests/test_tracked_files.py" or not path.is_file():
            continue
        text = path.read_text(errors="ignore")
        hits += [f"{rel}: {word}" for word in FORBIDDEN if word in text]
    assert not hits, hits
