"""What reaching each constrained front costs against the unconstrained run.

    E2_DEPOSIT=... E2_MOTIF_DIR=... E2_SUPP=... python supp_cost.py

Per cell: the wall summed over every round, the peak RSS maxed over them, and the number of
Pareto optimal designs, each as a ratio to the protein's unconstrained run.

Two things decide what those ratios mean.

A level started from the front of the level below did not pay for that front, so its own
wall is only the increment. The cost of REACHING a front therefore carries the chain: the
wall adds the level below's and the peak RSS is the larger of the two. Read without that,
a level that narrowed from the one below shows a wall of four ten-thousandths of the
unconstrained run, which is the increment and not a cost.

A cell's rounds are named `<cell>_<nest>.rN.jsonl` beside an unsuffixed first file, and
some cells carry only the unsuffixed one; reading just the `.rN` files drops those
silently, which is 63 of 163 cells on the table this was written against.
"""

import glob
import json
import os
import sys

import supp_common as common

# What the two hosts of the grid are called in the rows and the table.
FIRST_HOST = os.environ.get("E2_FIRST_HOST", "host-a")
SECOND_HOST = os.environ.get("E2_SECOND_HOST", "host-b")


def builds(cell, nest, d):
    """{round file: (host, narrowed-from)} from the cell's .builds log.

    The grid ran on two hosts and the early rounds did not record which. A line with no
    `host=` is attributed by its environment: one of the two wrote a home directory the
    other does not have, which `E2_SECOND_HOST_MARKER` names, and a line without that
    marker is the first host's. Both labels and the marker are the run's, so they come
    from the environment rather than being spelled here.
    """
    out = {}
    p = f"{d}/{cell}_{nest}.builds"
    if os.path.exists(p):
        for line in open(p):
            f = line.split()
            if len(f) < 2:
                continue
            kv = dict(x.split("=", 1) for x in f[2:] if "=" in x)
            host = kv.get("host")
            if host is None:
                marker = os.environ.get("E2_SECOND_HOST_MARKER", "")
                second = marker and marker in kv.get("env", "")
                host = SECOND_HOST if second else FIRST_HOST
            out[f[1]] = (host, kv.get("narrowed"))
    return out


def chain(cell, nest, rec):
    """Summed wall, max peak RSS, rounds, hosts and whether the cell narrowed."""
    d = os.path.dirname(rec["path"])
    files = set(glob.glob(f"{d}/{cell}_{nest}.r*.jsonl")) | {rec["path"]}
    plain = f"{d}/{cell}_{nest}.jsonl"
    if os.path.exists(plain):
        files.add(plain)
    b = builds(cell, nest, d)
    wall, rss, rounds, hosts, narrowed = 0.0, 0, 0, set(), False
    for p in sorted(files):
        try:
            h = common.header(p)
        except (json.JSONDecodeError, OSError, ValueError):
            continue
        if h.get("wall_s") is None:
            continue
        wall += h["wall_s"]
        rss = max(rss, h.get("peak_rss_kb") or 0)
        rounds += 1
        host, nar = b.get(os.path.basename(p), (FIRST_HOST, None))
        hosts.add(host)
        narrowed |= bool(nar)
    return wall, rss, rounds, hosts, narrowed


IMPORTED_JSON = os.environ.get("E2_F0_TIMINGS",
                               f"{os.path.dirname(common.E)}/f0_timings_imported.json")


def imported_timings():
    """{cell: (wall_s, peak_rss_kb)} for fronts E2 imported rather than ran itself.

    An imported front carries its points and no timing, so a ratio against it has no
    denominator. The file supplies the one the front was measured with where the run it
    came from recorded one. It is a different build on a different host, which the file
    states and which every row using it carries, because a ratio formed against it prices
    the constraint together with whatever differs between the two.
    """
    raw = json.load(open(IMPORTED_JSON))
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def main():
    tsv, digest = common.finished()
    imported = imported_timings()
    _, f0path = common.sequences()
    done = sorted(c for c, v in tsv.items()
                  if all(v.get(n, {}).get("status") in ("complete", "retired")
                         for n in common.NESTS) and c in f0path)
    rows, no_f0, borrowed = [], [], []
    for c in done:
        f0 = common.header(f0path[c])
        n_f0 = len(common.points(f0path[c]))
        wall0, rss0 = f0.get("wall_s"), f0.get("peak_rss_kb")
        source = "own"
        if not (wall0 and rss0) and c in imported:
            wall0, rss0 = imported[c]["wall_s"], imported[c]["peak_rss_kb"]
            source = "imported"
        timed = bool(wall0) and bool(rss0)
        if not timed:
            no_f0.append(c)
        elif source == "imported":
            borrowed.append(c)
        below_w, below_rss, below_hosts = 0.0, 0, set()
        for n in common.NESTS:
            w, rss, nr, hosts, nar = chain(c, n, tsv[c][n])
            if nar:
                w += below_w
                rss = max(rss, below_rss)
                hosts |= below_hosts
            below_w, below_rss, below_hosts = w, rss, set(hosts)
            rows.append({
                "cell": c, "level": n, "band": tsv[c][n]["band"],
                "status": tsv[c][n]["status"], "rounds": nr, "narrowed": nar,
                "hosts": sorted(hosts),
                "f0_timing": source,
                "wall_ratio": (w / wall0) if (timed and nr) else None,
                "rss_ratio": (rss / rss0) if (timed and nr) else None,
                "size_ratio": len(common.points(tsv[c][n]["path"])) / n_f0 if n_f0 else None,
            })
        print(c, "done", flush=True)
    common.write("cost.json", {
        "rows": rows,
        "provenance": {
            "n_proteins": len(done),
            "members": done,
            "untimed_unconstrained_runs": no_f0,
            # Named, not counted: a ratio against a borrowed denominator is not the same
            # measurement as the rest and a reader has to be able to say which one it is.
            "imported_unconstrained_timings": borrowed,
            "finished_set_sha256": digest,
            "scripts": common.digest_of(__file__, common.__file__),
        },
    })


if __name__ == "__main__":
    sys.exit(main())
