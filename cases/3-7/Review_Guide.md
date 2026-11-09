# Review guide: the (3,7) computer-assisted upper bound

**Status.** This is a candidate proof, not an independently peer-reviewed or
proof-assistant-verified theorem. The checker is exact and reads no input
certificate, but the mathematical state-coverage and perturbation arguments,
and the checker itself, require external review. The mathematical statement,
sharp examples, and general theorems are in `Proof_and_Sharp_Bounds.md`.

## What needs to be established

For a generic real `3 x 7` matrix, normalize `p567>0` and let `m` be the
minimum number of sign changes in a strict row-space sign vector. The goal is
a nonnegative integer `H` on a finite set of signs such that every actual
generic elementary-move event is represented and

* `H(state) <= 8 + m(state)` at every admitted state;
* `H(before) >= H(after) + 1` when `p123` changes sign; and
* `H(before) >= H(after)` at every other admitted event.

Telescoping these inequalities proves `r(M) <= 8+m(M)`. The separate
28-move witness proves the global lower bound 12. The source-only checker
constructs `H` and checks the three conditions directly; it does **not** need
the supplied older rank arrays or trust an optimality assertion about `H`.

## Source map

The minimal archive calls the checked-arithmetic single-file checker
`verify_upper_bound.cpp`; the full archive calls it
`grassmann_37_continuation/history_20261108/standalone_data_free_checked.cpp`. The full archive also
contains the original independently audited `standalone_data_free.cpp`.

| Source section | Mathematical obligation checked or implemented |
|---|---|
| `labels`, `bracket`, `relations` | Index the 35 Plücker coordinates and 15 disjoint-edge features; verify the 16 stated bracket identities exactly. |
| `direction`, `generate_keys` | Enumerate all sign words satisfying those necessary identities and possessing a consistent circular order of the ten affine edge rays; insert the two ordered axes in all allowed gaps and restore four column reorientations. |
| `check_affine_normalization` | Verify symbolically the exact determinant formulas that connect the 50 original signs to affine edge directions. |
| `symbolic_events`, `generic_separation` | Derive all 114 nonzero move derivatives as exact polynomials; check all inactive derivatives, the two G identities, and 1,070 same-generator event-separation resultants. |
| `gvalid` | Apply only necessary sign consequences of the two additional features, their defining differences of products, and exact Plücker identities. |
| `Variation` | Compute the minimum variation using the 35 four-column circuits and Gordan's alternative; independently compare the result with the sign-variation formula. |
| Graph and SCC code in `main` | Build a directed over-approximation, find its strongly connected components, and construct an integer candidate `H`. |
| Final parallel loop in `main` | Directly verify the bound and all move-labeled edge inequalities. This final loop is the actual numerical rank certificate. |

The checker rejects compilation with `-DNDEBUG`, since its exact assertions
are part of the proof check. The reference run needs roughly 2.4 GiB of RAM.
It uses exact integer operations; it performs no floating-point optimization.

## Coverage and events: the parts to scrutinize most

Write `v6=(0,1,0)`, `v7=(0,0,1)`, and `vi=xi(1,yi,zi)` for `i<=5`, using
the fixed sign of `p567` to have `x5>0`. After column reorientation, the
signs of the minors and the 15 `F` polynomials are precisely signs of
oriented brackets of the ten edge vectors of five planar points and the two
coordinate axes. Thus every real strict state has one of the cyclic edge
orders enumerated by `generate_keys`. The checker allows extra sign states;
the proof needs coverage, not realizability of every state.

The 15 disjoint-edge sign choices range over `2^15` possibilities; the
16 verified exterior identities leave 1,920. For each, the program tries all
`2^10` triangle sign choices and reconstructs/checks the full circular order.
It obtains 4,560 angular types. Twenty gaps for axis 6, eleven allowed gaps
for axis 7, and 16 reorientations produce 16,051,200 distinct base keys.
The additional G tests produce 23,246,848 augmented states. Independent code
in the full archive reconstructs the same base states by a different route.

Every tracked feature is affine in any one shear parameter. The checker
derives the actual derivative signs rather than reading a move table. It
adds an event edge when a feature and its derivative have opposite signs and
the resulting sign vector satisfies the necessary filters. The resulting
graph may contain paths that no real matrix can follow, which is safe for an
upper bound. It checks 230,629,888 **move-labeled** edges; 144,788,896
distinct feature-flip edges remain after move labels are merged.

For any fixed finite real trajectory, perturb the starting matrix slightly
while keeping the positive shear parameters. Every prefix is an invertible
linear map of the starting matrix, so finitely many proper algebraic sets
can be avoided simultaneously. The exact witness built into the checker
shows that no two active tracked polynomials for one generator have an
identically coincident affine root; endpoint feature values can also be
made nonzero. The perturbation preserves every strict sampled target sign
and all initial maximal-minor signs, hence preserves `m`. This is how the
event-graph bound applies to arbitrary generic starts and finite sequences.

## Independent checks and limits

The minimal `run.sh` checks the data-free graph, exact rational trajectories,
and a rational construction from a simple moment-curve start. Optional SymPy
scripts check the positive-cell chart in two different symbolic directions.
For a lower-memory check of the upper bound, `run_low_memory.sh` reads two
bundled rank arrays and verifies the bound and every transition inequality
directly; correctness of the arrays' generation is not assumed.
The full archive adds independent direction enumeration, direct rank checks,
variation/circuit checks, a generic-event audit, a second rank construction,
independent replays of the exact 28-move witness, and source-free rank-one and pencil lower-bound
proofs. The full script and retained logs identify each check. A combined
address/undefined-behavior sanitizer run was also completed.

The full archive additionally contains a 41-sign potential with 4,960,256
projected states. A direct checker tests it on all full states and all
230,629,888 labeled edges, including the pointwise bound. It is an additional
representation, not a proof that the hidden 52-sign enumeration can be skipped.

The 28-move lower bound for a 12-crossing example uses the separately
established `8+m` certificate, the published rank-two theorem, and an
independently checked 384-state five-column sign automaton. The 28-move
witness attains that bound. The all-dimensional positive-cell transfer and rank-one criterion do not imply
the conjectured upper bound for other `(k,n)`.
