"""Pack the deposit of one experiment: exactly the files its manifest names, nothing else.

    python package_deposit.py MANIFEST ROOT OUT_TAR PREFIX

The result directories also hold run logs and working scripts, which carry the host's paths
and are not part of the record; packing from the manifest leaves them out. Every file is
checked against its digest first.
"""

import sys
import tarfile
from pathlib import Path

from derna2_analysis.provenance import verify


def main():
    manifest, root, out, prefix = sys.argv[1:]
    files = verify(manifest, root)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(out, "w") as tar:
        for rel in files:
            tar.add(Path(root) / rel, arcname=f"{prefix}/{rel}", recursive=False)
    print(f"{out}: {len(files)} files from {manifest}")


if __name__ == "__main__":
    main()
