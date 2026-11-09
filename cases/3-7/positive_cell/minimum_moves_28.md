# Exactly 28 moves suffice and are necessary for twelve target flips

A generic `(3,7)` trajectory with 12 changes of `p123` requires **at least 28 positive elementary moves**, and an exact 28-move trajectory attains that minimum. The lower bound is computer-assisted: it uses the separately certified `(3,7)` pointwise rank bound `r(M)<=8+m(M)`, the published rank-two theorem, the elementary zero-future-flips theorem, and the small exact over-approximation below. The witness `../optimal_witness.json` and its all-integer version `../optimal_witness_integer.json` are independently replayed in `the independent checks in the full archive`; each has 12 target flips and all 35 maximal minors nonzero at every recorded vertex.

## A 384-state obstruction

A 12-flip trajectory must start with `m=4` by `r(M)<=8+m(M)`. Normalize `p567>0`; its 35 minor signs are then alternating. In particular, the first five columns start with
\[
 \operatorname{sgn}p_{ijk}=(-1)^{i+j+k}\quad(1\le i<j<k\le5).
\]
Track only the ten maximal-minor signs of those first five columns. For each `s in {1,...,5}`, let `a<b<c<d` be the other four indices. Every real state satisfies the necessary sign relation from
\[
 p_{sab}p_{scd}-p_{sac}p_{sbd}+p_{sad}p_{sbc}=0,\tag{1}
\]
where each minor is oriented in the written order. Of the `2^10` sign assignments, exactly **384** satisfy all five sign conditions, including both overall orientations. They are an over-approximation; their individual or joint realizability is not assumed.

For a move `U_j`, `j=1,...,4`, an affected minor has
\[
 p_I(t)=p_I(0)+t\,p_{(I\setminus\{j\})\cup\{j+1\}}(0),
 \qquad j\in I,\ j+1\notin I.
\]
The graph permits the single sign flip whenever current signs oppose and the resulting sign state satisfies (1). This gives **1,152 move-labeled edges**. An index-5 move changes only column 5 among the first five, so it leaves all four minors on columns `1,2,3,4` fixed. We relax such a move to *any* jump between two of the 384 states having those four signs in common. Index-6 moves are invisible to the first five columns and may be omitted from this record.

To model at most **one actual index-1 move**, use three modes: before the index-1 event run, inside it, and after it. A sign event labeled 1 enters or stays inside the run; a later event labeled 2, 3 or 4, or a relaxed index-5 jump, ends the run and forbids future index-1 events. A zero-cost transition may also end or skip the run. This allows arbitrarily many consecutive index-1 sign events, so one actual index-1 shear is certainly covered even if it crosses several minor walls.

Starting at the alternating sign state, exact longest-path relaxation in this graph gives the following upper bounds, allowing the final state and mode to be arbitrary:

| Relaxed index-5 jumps | Maximum target changes, with at most one index-1 run |
|---:|---:|
| 0 | 6 |
| 1 | 8 |
| 2 | 10 |

Thus a 12-flip trajectory with at most two index-5 moves must use **at least two** index-1 moves. Equivalently, every 12-flip trajectory satisfies
\[
 \#U_1\ge2\quad\text{or}\quad\#U_5\ge3.\tag{2}
\]

`verify_minimum_28_independent.py` constructs the 384 states and all labeled edges independently, without importing the earlier graph or rank data. Its three-mode, three-jump-count automaton reaches 2,496 states and gives the displayed maxima. It also verifies one exact integer matrix at which all 12 pairs of simultaneously active minors for moves 1–4 have distinct roots. Each corresponding root-resultant polynomial is therefore nonzero. For any fixed finite seven-column move word, a prefix is invertible and its first-five-column projection is surjective as the initial `3`-by-`7` matrix varies. Thus an arbitrarily small generic perturbation separates all intermediate minor events while retaining the alternating initial signs, positive move parameters, generic recorded vertices, and every strict recorded target sign change. The relaxed graph covers the perturbed trajectory and proves the bound for the original one.

## Mandatory moves of each kind

There are at least **12 index-3 moves**, because no other generator changes `p123`. Between consecutive target-changing index-3 moves there must be an index-2 or index-4 move: without one, `p124` is fixed and each index-3 move changes `p123` in the same direction. Thus there are at least **11 index-2/index-4 separators**.

At least **one index-1 move** is needed. If `v1` were fixed, projecting `v2,...,v7` modulo `Rv1` would give a six-vector rank-two trajectory whose cyclic prerank is at most 8.

At least **one index-6 move** is needed. If none occurs, the first six columns form an autonomous `(3,6)` trajectory and `v6` is fixed. Choose a functional giving at most three sign variations on its initial six columns; such a functional exists by the elementary `m<=n-k` lemma. Append a fixed `v7` near `v6`, on the same side of that functional, and avoid the finitely many planes making a recorded triple singular. The padded initial matrix has `m<=3`; the certified `(3,7)` bound therefore gives at most 11 target changes, a contradiction.

At least **two index-5 moves** are needed. Without index 5 the first five columns form a `(3,5)` trajectory, which has at most six target changes by rank-two duality. With exactly one index-5 move, the portions before and after it have at most six each, so 12 would force exactly six in each. The first portion's final first-five matrix has zero future rank, since any further target flip would exceed six. The zero-future-flips theorem makes all its maximal minors the same sign. The index-5 move leaves `p123` and `p124` unchanged, so the second portion begins with these two minors of the same sign. It is therefore nonalternating, and the strict rank-two dual bound gives at most five further changes. The Gale-duality proof and its primary-source citation are in `minimum_moves_27.md`.

By (2), the indices 1 and 5 together occur at least four times: either `(at least 2, at least 2)` or `(at least 1, at least 3)`. All the listed index classes are disjoint, so the total move count is at least
\[
 12\ (j=3)+11\ (j=2\text{ or }4)+1\ (j=6)+4\ (j=1\text{ or }5)=28.
\]

## Reproduction

From the workspace root:

```sh
python3 grassmann_37_continuation/history_20261108/verify_minimum_28_independent.py
python3 grassmann_37_continuation/new_work/probe_two_j5_five_sign_graph.py
python3 grassmann_37_continuation/history_20261108/verify_dual_rank2_3x5.py
```

The first two check the same over-approximation with separately written state and closure code; the first imports no prior program. Their finite check is in addition to, not a replacement for, the independent `(3,7)` pointwise certificate and the mathematical coverage argument above.
