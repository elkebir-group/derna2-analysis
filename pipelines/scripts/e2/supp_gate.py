"""The overlap gate: every pair of drawn labels that touch, in device pixels.

The composite plate (`fig_composite.py`) gates its own save this way and the supplement
plate is held to the same standard, because at 3 to 6 pt a collision is a few pixels and
survives a read-through -- (f)'s panel letter printed on (e)'s y label through two rounds of
looking at the figure, and was found by cropping the PNG rather than by reading it.

Kept here rather than imported from the composite: that script is a 1,500-line program that
builds its own figure on import, so importing it to borrow one function would run it. The
two copies are a duplicate of a test, not of a result, and either may tighten without the
other.

`owners` maps id(text) for figure-level texts to the axes each belongs to. A panel letter is
drawn on the figure but belongs to its panel; left unowned it is tested against the very box
it labels and reported as intruding on it.
"""


def _items(fig, owners):
    """(panel, group, text, bbox) for everything drawn that carries ink.

    `panel` is the axes a label may sit on without that being a collision; `group` is the
    set of labels laid out together, which are not compared with each other. The only
    groups are the two tick rows of an axes, because matplotlib places those and this is
    not the place to re-litigate its spacing, and whatever an artist's `gid` says it is
    laid out with -- which is how a caller declares that two lines are one block and their
    leading is deliberate. Everything else gets a group of its own, including two texts
    inside one axes: nothing positions those but the author, so a caption landing on a
    label it shares a panel with is exactly the kind of collision worth catching, and a
    panel letter is still compared with the title it stands beside.
    """
    items = []
    for ax in fig.axes:
        key = f"{id(ax)}"
        for i, t in enumerate([ax.xaxis.label, ax.yaxis.label, ax.title] + list(ax.texts)):
            if t.get_visible() and t.get_text().strip():
                items.append((key, t.get_gid() or f"{key}:t{i}",
                              t.get_text().strip()[:28], t.get_window_extent()))
        # A tick outside the view still has a label artist with a window extent and is not
        # drawn, and an axes with axis("off") keeps its tick labels and still reports them
        # visible: counting either invents collisions that are not on the page.
        # Tick labels carry their own owner per axis, so matplotlib's own layout of a tick
        # row is not re-litigated here, but they ARE compared with this axes' axis labels.
        for axis, lim, who in (() if not ax.axison else
                               ((ax.xaxis, ax.get_xlim(), "tx"),
                                (ax.yaxis, ax.get_ylim(), "ty"))):
            lo, hi = sorted(lim)
            for loc, lab in zip(axis.get_ticklocs(), axis.get_ticklabels(), strict=True):
                if lo - 1e-9 <= loc <= hi + 1e-9 and lab.get_visible() \
                        and lab.get_text().strip():
                    items.append((key, f"{key}:{who}", lab.get_text().strip()[:28],
                                  lab.get_window_extent()))
    for i, t in enumerate(fig.texts):
        if t.get_text().strip():
            items.append((owners.get(id(t), "fig"), f"fig:{i}",
                          t.get_text().strip()[:28], t.get_window_extent()))
    for lg in fig.legends:
        items.append(("fig", "legend", "<legend>", lg.get_window_extent()))
    return items


def _edges(fig, items, clear):
    """Anything drawn outside the page, or nearer an edge than the clearance.

    With no tight box an overrun is cut off silently, so the outside test has to be
    explicit. The margin is held to the same clearance as everything else: a label a
    pixel inside the canvas is not off it, and still reads as falling off the page.
    """
    w, h = fig.get_size_inches() * fig.dpi
    return [(text, "<outside the page margin>", round(bb.x0, 1), round(bb.y1, 1))
            for _p, _g, text, bb in items
            if bb.x0 < clear or bb.y0 < clear or bb.x1 > w - clear or bb.y1 > h - clear]


def _near(a, b, clear):
    """Overlap in both axes, counting a gap narrower than `clear` as a touch."""
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return (round(dx, 1), round(dy, 1)) if dx > -clear and dy > -clear else None


def _into_panels(fig, items, clear):
    """A label straying into a panel it does not belong to.

    Text against text cannot see this, and it is the class that bit both plates: a label can
    sit centred over a neighbour and read as that panel's caption without touching any other
    label. Clearance rather than penetration, for the same reason as below.
    """
    boxes = {f"{id(ax)}": ax.get_window_extent() for ax in fig.axes}
    bad = []
    for panel, _g, text, bb in items:
        for other, ob in boxes.items():
            hit = None if other == panel else _near(bb, ob, clear)
            if hit:
                bad.append((text, "<a panel it is not part of>", *hit))
    return bad


def _pairs(items, clear):
    """Labels that merely come close read as touching at this size."""
    bad = []
    for i, (_p, gi, ti, bi) in enumerate(items):
        for _q, gj, tj, bj in items[i + 1:]:
            hit = None if gi == gj else _near(bi, bj, clear)
            if hit:
                bad.append((ti, tj, *hit))
    return bad


def overlaps(fig, owners=None, clear_pt=2.9):
    """Every collision on the page, as (what, with what, dx, dy) in device pixels.

    The clearance is given in POINTS and converted here. A bbox is measured in device
    pixels, so a threshold written in pixels means a different distance on the page at
    every `figure.dpi`: the 4 px this was first written with is 2.9 pt at the matplotlib
    default of 100 and 1.4 pt at the 200 the figure now fixes, which would have relaxed
    the test by half without a line of it changing.
    """
    fig.canvas.draw()
    clear = clear_pt * fig.dpi / 72.0
    items = _items(fig, owners or {})
    return (_edges(fig, items, clear) + _into_panels(fig, items, clear)
            + _pairs(items, clear))


def check(fig, owners=None, clear_pt=2.9):
    """Report every collision and refuse to go on if there is one.

    An `assert`, so that `python -O` cannot drop the check: printed and not acted on, this
    would be a number beside a figure that shipped either way, which is not a check.
    """
    bad = overlaps(fig, owners, clear_pt)
    for ti, tj, dx, dy in bad:
        print(f"  OVERLAP  {ti!r} x {tj!r}  by {dx} x {dy} px")
    print(f"overlapping label pairs: {len(bad)}")
    assert not bad, f"{len(bad)} label collision(s); the figure was not written"
