"""Report, and on request rewrite, the host paths a deposit's raw outputs carry.

    python scrub_host_paths.py ROOT OLD=NEW [OLD=NEW ...] [--apply]

A run's logs and per-round reports name the checkout they ran in, which is a property of
the host and not of the result. Without `--apply` nothing is written: the report names
every file and how many occurrences each substitution has, so the substitutions can be
read before they are made. With `--apply` the files are rewritten in place and the deposit
is re-checked, so a path that no substitution covers is an error rather than a surprise in
the published tarball. A rewritten deposit's manifest no longer matches it; regenerate the
manifest afterwards.
"""

import re
import sys
from pathlib import Path

# Any absolute path under a user's home or the scratch volume, which a published file
# should not carry whatever produced it.
HOST_PATH = re.compile(r"/(?:home|scratch)/[^\s\"',)]+")


def read(path: Path) -> str | None:
    """The file's text, or None if it is not text."""
    try:
        return path.read_text()
    except (UnicodeDecodeError, OSError):
        return None


def counts(text: str, subs: dict[str, str]) -> dict[str, int]:
    return {old: text.count(old) for old in subs if old in text}


def substitute(text: str, subs: dict[str, str]) -> str:
    for old, new in subs.items():
        text = text.replace(old, new)
    return text


def report(root: Path, subs: dict[str, str], apply: bool) -> list[str]:
    """Print one line per file a substitution covers; rewrite them if asked. Returns the
    host paths still present afterwards."""
    left: set[str] = set()
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        text = read(path)
        if text is None:
            continue
        found = counts(text, subs)
        if found:
            totals = " ".join(f"{n}x" for n in found.values())
            print(f"{path.relative_to(root)}: {totals}")
            text = substitute(text, subs)
            if apply:
                path.chmod(0o644)
                path.write_text(text)
                path.chmod(0o444)
        left.update(HOST_PATH.findall(text))
    return sorted(left)


def parse(args: list[str]) -> tuple[Path, dict[str, str], bool]:
    apply = "--apply" in args
    rest = [a for a in args if a != "--apply"]
    if len(rest) < 2 or not all("=" in a for a in rest[1:]):
        sys.exit(__doc__)
    subs = dict(a.split("=", 1) for a in rest[1:])
    return Path(rest[0]), subs, apply


def main():
    root, subs, apply = parse(sys.argv[1:])
    left = report(root, subs, apply)
    print(f"{'rewrote' if apply else 'would rewrite'} under {root}")
    if left:
        print(f"host paths no substitution covers: {len(left)}", file=sys.stderr)
        for path in left[:20]:
            print(f"  {path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
