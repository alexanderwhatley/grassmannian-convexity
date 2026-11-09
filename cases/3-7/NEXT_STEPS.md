# Grassmann convexity (3,7): continuation notes

The main report describes a **candidate computer-assisted proof** of
`r(3,7)=12` and the pointwise bound `r(M)<=8+m(M)`. It has not been
independently peer reviewed or formalized in a proof assistant. In particular,
the finite-state coverage, the treatment of sign-wall degeneracies, and the
source code remain review obligations.

## Established within the candidate proof

* The data-free checker constructs 23,246,848 sign states, checks
  230,629,888 move-labeled edge inequalities exactly, and finds a nonnegative
  rank bounded by `8+m`. The source with checked signed arithmetic has
  SHA-256 `ad5acbcba71e601f16eb628982988088ae953d797a83794304ece4170ab36946`.
* A second potential depends on 41 signs and has 4,960,256 projected states;
  a direct check covers all original 52-sign states and move-labeled edges.
  Its current construction still enumerates the hidden full states.
* Exact rational and all-integer 28-move paths have 12 target changes; the
  integer example is printed in the report. The simple alternating moment-curve
  matrix also has an exact 28-move path. All have independent exact replays.
  Thus the upper bound is sharp.
* Every generic `m=4` start has a 12-change path of at most 40 moves: a
  12-shear positive-cell transfer followed by the 28-move witness. Hence,
  subject to the checked upper bound, `r(M)=12` exactly for `m(M)=4`.
  The exact rational constructor handles any rational generic `m=4` input.
* Every 12-change path needs at least 28 moves, and the 28-move witness
  attains it. The 384-state obstruction for the lower bound has multiple
  independent implementations. Its reduction to the full problem uses the
  certified `8+m` bound and the rank-two theorem.
* In all dimensions, alternating/maximum-variation starts have the same
  future rank, which is the global rank; the positive-cell chart and this
  transfer have an elementary proof. Also `r(M)<=1` exactly when the
  negative maximal-minor indices form a lower ideal, and ordered pencils of
  minors give a general sign-only lower bound. None of these general results
  establishes the Grassmann convexity upper bound for other dimensions.

## Outstanding mathematical directions

1. **Independent proof audit.** Reimplement the 50-sign coverage and the two
   G filters from the mathematical definitions, review generic-event
   separation, and independently verify the rank inequalities. The source-only
   and lower-memory checkers both avoid trusting precomputed graph paths, but
   their coverage argument is still a mathematical dependency. See
   `Review_Guide.md` and the independent audits in the full archive.
2. **Understand optimal 28-move paths.** Every such path needs exactly twelve
   U3 target flips, eleven U2/U4 separators, two U1 moves, two U5 moves and
   one U6 move; each of U2 and U4 appears 4 to 7 times. Both U5 moves must
   flip p345, and the unique U6 flips p456. Their order is U5, U6, U5.
   A 3x6 sign relaxation excludes the alternative
   one-U1/three-U5 allocation. Characterize the realizable words and
   investigate whether the 12-move preparation can be shortened uniformly. The lower-bound checker is
   `verify_minimum_28_independent.py`.
3. **Is the graph rank actual?** For realizable strict states with `H<=3`,
   the rank-one/terminal theorems and an exhaustive tally imply `r=H`.
   No general equality is asserted for `H>=4`. A concrete unresolved `H=4`
   matrix, its exact two-flip path, rank-two upper bound 4, and carefully
   limited failed searches are in `sample289_open_rank_case.md`: its true
   rank is presently only known to be 2 or 4.
4. **Human-scale upper bound.** The `m=0` case follows from the published
   rank-two theorem; a suffix-halfspace potential handles another class.
   Coarse projected graphs that retain only 35 Plücker signs admit abstract
   rank 24, while some tested 38–41-sign models still admit rank 16. The
   successful 41-sign potential is smaller but has no known human-scale
   formula. Those failed models are over-approximations, not actual
   counterexamples; see
   `human_potential_obstructions.md` and the projection logs. For the 3-by-5
   sign graph, an exact check shows its rank is `V` or `V+2`, where `V` is
   the variation of `(p123,p124,p134,p234,p345)`. A nonenumerative proof
   of this formula, or an analogue for seven columns, is a concrete target.

## Reproduction

The minimal archive's `./run.sh` builds the checked source and verifies the
integer witnesses, the constructor, the small move-count obstruction, and
the full upper certificate. `./run_low_memory.sh` verifies the upper bound
using two bundled compressed rank arrays and about 565 MiB on the reference
machine; it checks every inequality directly. The complete archive also
contains the earlier rank construction, independent audits, experiments,
`./run_all.sh`, and `./run_supplement.sh` for the 41-sign and 28-move
refinements. Assertions must remain enabled. Symbolic chart audits in the
complete archive require SymPy; the minimal core needs only Python's
standard library and a GCC/Clang-compatible C++17 compiler.
