# The pinned build for E3

    06605630   bin/derna_06605630   adopted 2026-10-21T22:20-0500
    549074e9   bin/derna_549074e9   previous pin, rounds before that time

Adopted on two independent grounds, both required.

CORRECTNESS. Fronts are bit-identical to the previous pin at F0, F1 and F4 on
P54852: 486/486, 479/479 and 507/507 points paired at exactly equal mfe, worst
cai_log deviation 0.000e+00, zero designs differing, and emission order
identical. Compared with a bijection at equal mfe rather than set equality on a
float-bearing key, because exact comparison on cai_log tests IEEE addition order
and a tolerance wide enough to absorb that also admits a real difference.

SPEED, measured at the width and on the socket the grid actually runs -- four
threads, node 1, cores 9,11,13,15, two runs per build alternated:

    arm     wall s   user s   occupancy   instructions
    pin_a    127.8    430.4       84.2%     2307.075 G
    pin_b    128.2    429.9       83.8%     2307.290 G
    new_a    118.8    414.2       87.1%     2256.597 G
    new_b    119.5    416.4       87.1%     2256.685 G

    wall 0.931x, user CPU 0.965x, instructions 0.9781x
    occupancy 84.0% -> 87.1%

The two runs per build are what licenses reading the 6.9%: a single wall run on
this host has a 2.8-3.6% spread, which once returned the wrong SIGN for a 3%
effect. Here each build's own two runs agree to 0.3% and 0.6%.

Why occupancy is reported beside instructions. The commit ladder between the two
pins is non-monotonic in occupancy -- 86.5, 87.5, 73.9, 62.2, 74.4 -- and the
commit that costs the most takes +7.4% instructions AND -13.4 pp occupancy, so
it is not a trade. An idle worker retires no instructions and burns no cycles, so
an instruction-only A/B reports that regression as an improvement. Two
consequences: a fix must be measured as a build rather than as a revert, and an
A/B between the two ends of a commit range cannot distinguish no improvement
from an improvement cancelled by a regression later in the range -- an endpoint
A/B over this range returned "do not move" for a range containing this 6.9%.

A front composed across rounds may be composed across builds. That is sound --
both builds compute the exact front, so inherited points and newly walked boxes
are each exact -- but a row's single build column would then name one build for
a front produced by two. Each round appends its own build to <label>.builds, and
that file is the authority on what produced a front.
