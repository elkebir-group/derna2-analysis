"""A DOI in the README is either real or absent, never a placeholder.

This exists because the placeholder shipped once already. The sibling paper's analysis
repository, `elkebir-group/H2-opt-analysis`, is public and its README points a reader at
`https://doi.org/10.13012/B2IDB-XXXXXXX_V1`, which resolves to HTTP 404: the deposit's
identifier was written as a template, the repository went public, and nobody looked again.
A reader following it finds nothing, which is worse than a README that says the deposit is
not published yet.

The assertion is deliberately conditional rather than a demand for a DOI. Requiring one
would fail from now until the Illinois Data Bank record exists, which is a true statement
made at the wrong time and in the wrong place: the README's open items are where unfinished
work is recorded, not a red test. This fails only when a DOI-shaped string is present *and*
looks unfilled, so it is silent today and speaks on the day someone pastes a template in.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Any Illinois Data Bank or generic DOI, however it is written: bare, or as a doi.org URL.
DOI = re.compile(r"(?:doi\.org/|\bdoi:\s*)(10\.\d{4,9}/\S+)", re.I)

# What an unfilled identifier looks like. `X{3,}` catches B2IDB-XXXXXXX; the rest catch the
# other ways a template is spelled. A real IDB DOI is 10.13012/B2IDB-<7 digits>_V<n>.
PLACEHOLDER = re.compile(
    r"X{3,}|<[^>]*>|\bTBD\b|\bTODO\b|\bFIXME\b|\bPLACEHOLDER\b|0{5,}|\bNNNN+\b", re.I
)

DOCS = ["README.md", "CITATION.cff", "NOTICE"]


def test_no_readme_doi_is_a_placeholder():
    bad = []
    for name in DOCS:
        path = ROOT / name
        if not path.is_file():
            continue
        for doi in DOI.findall(path.read_text()):
            if PLACEHOLDER.search(doi):
                bad.append(f"{name}: {doi}")
    assert not bad, (
        "a DOI is written as an unfilled template; publish the deposit and paste the "
        f"reserved DOI, or say the deposit is unpublished and give no DOI: {bad}"
    )


def test_the_check_can_fail():
    # The mirror of the test above: an instrument that cannot fail is not a check, and this
    # one's whole value is that it fires on a string nobody has written yet.
    assert DOI.findall("see https://doi.org/10.13012/B2IDB-XXXXXXX_V1 for the data")
    assert PLACEHOLDER.search("10.13012/B2IDB-XXXXXXX_V1")
    assert PLACEHOLDER.search("10.13012/B2IDB-<reserved>_V1")
    # and that it passes a real one
    assert not PLACEHOLDER.search("10.13012/B2IDB-1234567_V1")
