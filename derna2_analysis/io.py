"""Result rows: JSON lines, one row per run, as the benchmark drivers write them."""

import json
from pathlib import Path


def read_jsonl(path):
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def write_jsonl(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")


def rows_of_pin(rows, pin):
    """The rows one build produced. A results file accumulates across re-pins, and rows of
    two builds are not one experiment, so a reader picks one pin. Rows written before the
    driver stamped them carry no pin and belong to none."""
    return [r for r in rows if r.get("pin") == pin]


def expect_count(rows, expected, what):
    """Refuse to read a results file whose row count is not the planned one: a short file
    reads exactly like a small grid."""
    if len(rows) != expected:
        raise ValueError(f"{what}: {len(rows)} rows, expected {expected}")
    return rows
