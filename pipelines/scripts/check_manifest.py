"""Check files against a sha256 manifest and write a stamp listing them.

    python check_manifest.py MANIFEST ROOT STAMP

A rule that reads binaries or deposited results depends on the stamp, so it never runs on
a build or a deposit other than the one pinned.
"""

import sys
from pathlib import Path

from derna2_analysis.provenance import verify


def main():
    manifest, root, stamp = sys.argv[1:]
    checked = verify(manifest, root)
    Path(stamp).parent.mkdir(parents=True, exist_ok=True)
    Path(stamp).write_text("".join(f"{rel}\n" for rel in checked))
    print(f"{manifest}: {len(checked)} files match under {root}")


if __name__ == "__main__":
    main()
