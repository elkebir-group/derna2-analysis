"""The numbers the UTRs subsection of the paper quotes, as one text file.

    python story_stats.py BUNDLE E2_RESULTS HULLCHECK CAND CELL CAI_FLOOR OUT.txt

Arguments are those of fig_story.py, which draws the figure these numbers go with. Each line is
`key: value`, grouped by the paragraph that quotes it:

  concatenation  its designs are the hull vertices of its front, as in panels (a) and (b)
  mRNAfold       its designs are the ones its sweep returned
  DERNA2         the modular designs of its output, and the windows of its search by state
  instance       CELL at CAI_FLOOR: the designs of the three methods, the windows, and the bound
  instances      the instance's comparison over every cell
"""
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_story as F  # noqa: E402
from e3common import a0_path, rd, relocate  # noqa: E402


def pct(a, b):
    return f"{100 * a / b:.1f}%"


def main():
    bundle, e2, hullcheck, cand, cell, floor, out = sys.argv[1:8]
    floor = float(floor)
    F.CELL, F.THRESHOLD = cell, floor
    cells, head, masked, window, a4, hull, mf, refolded = F.load(bundle, e2, hullcheck)
    assert len(refolded) == len(cells)
    sets = {c: F.modular_sets(c, masked, window, hull, mf) for c in cells}
    combo, share = F.element_census(sorted(refolded), head, hull)
    split = F.pareto_split(cells, masked, window, sets)
    on = lambda p: [c for c in cells if c.endswith("__" + p)]  # noqa: E731
    lines = []

    def say(key, value):
        lines.append(f"{key}: {value}")

    say("cells", len(cells))
    say("proteins", len({c.split("__")[0] for c in cells}))
    say("residues min max", f"{min(len(head[c]['protein']['sequence']) for c in cells)} "
                            f"{max(len(head[c]['protein']['sequence']) for c in cells)}")

    # ---- concatenation
    total = sum(s[0] for s in share.values())
    nonmod = sum(s[1] for s in share.values())
    say("concatenation designs", total)
    say("concatenation non-modular designs", f"{nonmod} ({pct(nonmod, total)})")
    for p in F.PAIRS:
        say(f"concatenation modular share {p}", pct(share[p][0] - share[p][1], share[p][0]))
    classes = sum((combo[p] for p in F.PAIRS), Counter())
    say("non-modular designs with a pair at the Kozak sequence", pct(sum(v for k, v in classes.items() if "kozak" in k), nonmod))
    say("non-modular designs with a pair at a PRE", pct(sum(v for k, v in classes.items() if "pyr" in k), nonmod))
    say("concatenation modular designs", sum(split["hull"][p][0] for p in F.PAIRS))
    say("concatenation modular designs dominated", sum(split["hull"][p][2] for p in F.PAIRS))
    say("instances with a modular design: concatenation mRNAfold DERNA2",
        " ".join(str(sum(bool(sets[c][m]) for c in cells)) for m in ("hull", "mrnafold", "window")))
    for p in F.PAIRS:
        say(f"proteins with a modular design {p}: concatenation mRNAfold DERNA2",
            " ".join(str(sum(bool(sets[c][m]) for c in on(p))) for m in ("hull", "mrnafold", "window")))

    # ---- mRNAfold
    n_mrf = sum(len(mf[c]) for c in cells)
    m_mrf = sum(1 for c in cells for x in mf[c] if x["gap"] == 0)
    say("mRNAfold designs", n_mrf)
    say("mRNAfold modular designs", f"{m_mrf} ({pct(m_mrf, n_mrf)})")
    say("mRNAfold non-modular designs", f"{n_mrf - m_mrf} ({pct(n_mrf - m_mrf, n_mrf)})")

    # ---- DERNA2
    say("DERNA2 modular designs", sum(split["window"][p][0] for p in F.PAIRS))
    states = {p: sum((F.window_states(a4[c]) for c in on(p)), Counter()) for p in F.PAIRS}
    every = sum(states.values(), Counter())
    say("windows: full empty partial unstarted", f"{every['enum1']} {every['enum0']} {every['part']} {every['none']}")
    say("windows", sum(every.values()))
    for p in F.PAIRS:
        n = sum(states[p].values())
        say(f"windows {p}: all full", f"{n} {states[p]['enum1']} ({pct(states[p]['enum1'], n)})")
    hbb = on("hbb")
    say("hbb instances whose masked front has no modular design",
        sum(1 for c in hbb if not any(x["gap"] == 0 for x in masked[c])))
    say("hbb instances that are a single window", sum(1 for c in hbb if sum(F.window_states(a4[c]).values()) == 1))

    # ---- the instance
    say("instance", f"{cell} at CAI >= {floor}")
    say("instance concatenation hull vertices: all modular", f"{len(hull[cell])} {sum(1 for x in hull[cell] if x['gap'] == 0)}")
    say("instance mRNAfold designs: all modular", f"{len(mf[cell])} {sum(1 for x in mf[cell] if x['gap'] == 0)}")
    say("instance DERNA2 modular designs", len(sets[cell]["window"]))
    st = F.window_states(a4[cell])
    say("instance windows: all full empty partial unstarted", f"{sum(st.values())} {st['enum1']} {st['enum0']} {st['part']} {st['none']}")
    # In the result file `cai_log` is the total log adaptation over the protein: divide by the
    # residues to get the log CAI of a point.
    n_res = len(head[cell]["protein"]["sequence"])
    t = math.log(floor) * n_res
    walked = {}
    for s in a4[cell]["stretches"]:
        if s.get("walked"):
            k = (s["cai_log_lo"], s["cai_log_hi"])
            if k not in walked or s["round"] >= walked[k]["round"]:
                walked[k] = s
    inf = float("inf")
    for (lo, hi), s in sorted(walked.items(), key=lambda kv: (kv[0][0] is not None, kv[0][0] or 0.0)):
        if (-inf if lo is None else lo) <= t <= (inf if hi is None else hi):
            state = "full" if s["complete"] and s["points"] > 0 else "empty" if s["complete"] else "partial"
            say("instance window containing the cutoff: CAI range, state, points",
                f"{0.0 if lo is None else math.exp(lo / n_res):.4f} {1.0 if hi is None else math.exp(hi / n_res):.4f} {state} {s['points']}")
    mod = sorted(x["cai"] for x in masked[cell] if x["gap"] == 0)
    below = max(c for c in mod if c < floor)
    above = min(c for c in mod if c > floor)
    say("instance modular masked-front designs bracketing the cutoff", f"{below:.5f} {above:.5f}")
    say("instance masked-front designs between them: modular all",
        f"{sum(1 for c in mod if below < c < above)} {sum(1 for x in masked[cell] if below < x['cai'] < above)}")
    picks, _, _ = F.best_designs(bundle, e2, cand)
    for name, x in picks.items():
        say(f"instance best design {name}: CAI MFE CDS-UTR pairs",
            f"{x['cai']:.4f} {x['free_mfe']:.1f} {len(x.get('violations') or [])}")
    lowest = min(x["mfe"] for x in masked[cell] if x["cai"] >= floor)
    say("instance lowest UTR-masking front MFE at the cutoff", f"{lowest:.1f}")
    say("instance DERNA2 design MFE equals it", abs(picks["DERNA2"]["mfe_nopair"] - lowest) < 0.05)

    # ---- the comparison of the instance over every cell with a modular DERNA2 design at the cutoff
    having = cat_mod = mrf_mod = 0
    for c in cells:
        rows = list(rd(a0_path(bundle, c)))
        src = relocate(rows[0]["source"], bundle, e2)
        hv = {i for i, x in enumerate(r for r in rd(src) if r.get("record") == "point") if x.get("on_hull")}
        cc = [x for x in rows if x.get("record") == "point" and x["index"] in hv and x["cai"] >= floor]
        mm = [x for x in mf[c] if x["cai"] >= floor]
        dd = [x for x in sets[c]["window"] if x["cai"] >= floor]
        if not (cc and mm and dd):
            continue
        having += 1
        gap = {x["index"]: x["gap"] for x in hull[c]}        # the census gap, or the refold's where it has none
        cat_mod += gap[min(cc, key=lambda x: x["f0_mfe"])["index"]] == 0
        # Ranked on the free fold, the one energy every mRNAfold design has and the one the
        # figure plots this arm by: the census leaves the constrained energy null wherever the
        # free fold forms a forbidden pair, which is 2,572 of the 9,312 designs at this cutoff.
        mrf_mod += min(mm, key=lambda x: x["free"])["gap"] == 0
    say("instances with a modular DERNA2 design at the cutoff", f"{having} of {len(cells)}")
    say("  of these, the best design is modular: concatenation mRNAfold", f"{cat_mod} {mrf_mod}")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
