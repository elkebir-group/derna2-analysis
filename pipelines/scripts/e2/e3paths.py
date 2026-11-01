"""Where E2's composite-figure steps read the deposit and write their cache.

Two paths, and neither is a constant in a source file. They come from the environment rather
than from argv because these steps exec one another's prefixes -- composite_data.py runs the
head of draft_panel_a.py, dcai_at_mfe.py the head of pick_example3.py -- and an exec'd prefix
cannot be given its own arguments. A rule names both; there is no default, so a step run with
neither set fails at import rather than reading somebody's working directory.

E2_DEPOSIT is the experiment's results tree (`source("e2")`), E2_CACHE the directory the
steps write their JSON into and the figure reads it back from.
"""
import os

E = os.environ["E2_DEPOSIT"]
OUT = os.environ["E2_CACHE"]

# Panel (a)'s screen reaches outside E2's own deposit: the catalog is scored over the
# proteins E1 ran, and the motif census over an auxiliary tree of E1 instances screened for
# this figure. Both are named here rather than hard-coded, and both still have to be
# deposited before this chain is reproducible off this host.
E1_SWEEP = os.environ["E2_E1_SWEEP"]
E1_EXTRA = os.environ["E2_E1_EXTRA"]
# The motif patterns' display names, a static table that ships with the repository.
MOTIF_NAMES = os.environ["E2_MOTIF_NAMES"]
