# The long band's supply, at 1000 aa and above

Assembled 2026-10-21 from two sources, and it is one protein short of comfortable.

    15 fastas   docs/reference/_cpp/derna-main/data/uniprotSeq  (>= 1000 aa)
     1 fasta    e2/out/hull/PRT-04.fasta                        (1053 aa)
    ----
    16 total, 1006-2828 aa

The mirror holds 102 proteins at 80-2876 aa and is restored by `just sync-cpp`;
its members are discovered by scanning the directory rather than enumerated, so
this count follows the corpus rather than freezing it.

TWO GAPS, both material.

PRT-05, PRT-06 and PRT-07 are not on disk. The plan counts four Cas orthologs at
1053-1629 aa and only PRT-04 is present, so three of the four are missing and the
band's longest Cas instances are absent. PRT-03 at 984 aa belongs to the middle
band by length and is not here.

ReverseCovid is a reversed synthetic sequence, not a protein. It is kept in the
corpus because the census over it is free and informative as a control, and it is
excluded from selection: a band characterizing the price of motif-freedom on real
proteins cannot reach its tenth instance by admitting a sequence that is not one.
Counting it is how "15 mirror proteins above 1000 aa" reaches fifteen; without it
the mirror supplies fourteen.

## Amended 2026-10-21T22:50, both gaps closed

PRT-05 (1082 aa), PRT-06 (1368 aa) and PRT-07 (1629 aa) were sourced from
`derna2-analysis/inputs/prt/`, so the supply is 19 instances at 1006-2828 aa and
the Cas ladder is whole. PRT-04's two copies were compared before either was
trusted and are byte-identical.

Eligibility was amended to `at least one entry per set` before any front ran.
The F4 tier is now 15 natural instances against ten wanted. PRT-05 and PRT-06
reach F4, so the tier holds Cas instances; PRT-04, O75151 and PRT-07 do not
spell RCCATGG and stop at F3, which means the band's LONGEST Cas ortholog
contributes three nests rather than four.

One limitation the selection cannot repair: the supply has no instance between
1871 and 2828 aa. Two log strata over that span are therefore empty and fall
back to the nearest unused qualifier, both of which sit near 1300 aa. The band
spans 1049-2828 aa at 2.70x, and its coverage inside that span is uneven by
supply rather than by choice.
