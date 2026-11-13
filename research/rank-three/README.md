# Rank-three research candidates

These unrefereed drafts study positive elementary shears on real 3-by-n matrices. They are not a proof of the full Grassmann convexity conjecture, nor are they independently peer-reviewed. All claims and the published inputs should be checked by independent experts.

## Quadratic bound for arbitrary finite words

For `q=n-3>=2`, the [quadratic proof candidate](independent_main_proof.md) argues that the number of strictly witnessed sign changes of `p123` is at most `2q^2-3q+4`; it gives 3 for `q=1`. This is a finite discrete statement, not the conjectured bound `3q`. Its all-source reduction, the planar rank formula, the pole insertion bounds, and the interface ledger are essential. The published finiteness paper gives the general rank-three bound `F < (q+1)^3`; subject to mathematical review, the candidate is quadratic in q.

The main draft links the supporting [planar insertion argument](insertion_lemma_independent.md), [elementary shear analysis](planar_shear_elementary.md), and [bridge to the published planar theorem](planar_source_bridge.md). The [chart audit](chart_redteam_notes.md) highlights a corrected nonconsecutive minor. The draft addresses strict sign changes; the distinct-zero and multiplicity questions for continuous curves require separate hypotheses and review.

## Exact maximum for reduced index words

The [reduced-word construction](reduced_sharp_all_n.md) proposes that, among Coxeter-reduced positive shear words, the largest attainable number of target flips is `ceil(11n/6)-5`. The upper bound uses the published Gaetz–Gao–Jiradilok–Nenashev–Postnikov sorting-network theorem; the new proof candidate constructs matrices and reduced words attaining it. Arbitrary words may backtrack, so this does not establish a general linear bound.

## Finite checks

From the repository root with Python 3 (standard library only):

```sh
python3 research/rank-three/reduced_sharp_all_n.py 18
python3 research/rank-three/insertion_exhaustive_check.py
```

The first constructs and checks exact integer reduced-word examples for every n from 4 through 18, including strictness of all vertex maximal minors. The second exhausts 50,360 oriented ray orders through seven rays and checks 347,736 deletions, including independent rank cross-checks. These are finite diagnostics; they do not prove the all-width claims. The scripts and notes originate in an internal Sprint 12 research checkpoint; the rank cross-check helper was extracted and its import path adjusted for this repository. Neither separate internal implementations nor replays constitute external peer review.
