# Further restrictions on an optimal 28-move trajectory

The audited lower bound 28 and the exact 28-move witness in
`../optimal_witness.json` establish the
minimum number of positive shears needed for 12 target changes. The
3-by-6 certificate in `../Proof_and_Sharp_Bounds.md` shows every such
optimal trajectory has
\[
 (N_1,N_3,N_5,N_6)=(2,12,2,1),\qquad N_2+N_4=11.\tag{1}
\]
The additional restrictions below hold for **every** generic optimal
trajectory, not merely for the saved witness. The `p345` statement uses a
separate small finite sign check; the remaining statements are elementary
consequences of (1) and the global 12-change upper bound.

Normalize `p567>0`. A 12-change start has maximum minimum variation and is
therefore alternating, so in particular `p123>0`, `p124<0`, `p345>0` and
`p456<0` initially. The invariant last minor stays positive. After the
12th target change no future target change is possible by the global upper
bound. The terminal dominance lemma then forces **all** maximal minors
positive at the endpoint.

## The 11 separators

Every one of the 12 `U3` moves flips `p123`. During `U3`,
`p123 -> p123+t p124`, while `p124` stays fixed. Its sign must reverse
between consecutive target changes. Only `U2` and `U4` change `p124`, and
it changes affinely during each such move. There are exactly 11 of them;
therefore one lies between each consecutive pair of target-changing moves,
and each flips `p124`. There are no `U2` or `U4` moves before the first or
after the last target change.

Write the 11 separator labels as a word in `2` and `4`. Two adjacent `2`
separators must have a `U1` strictly between them: their `p124` derivatives
are `p134`, whose sign must change, and only `U1` or `U4` can change it.
Similarly, two adjacent `4` separators must have a `U5` strictly between
them: their derivatives are `p125`, which only `U2` or `U5` can change.
The intervals for distinct adjacent pairs are disjoint. Since there are two
`U1` and two `U5` moves, the separator word has at most two adjacent `22`
pairs and at most two adjacent `44` pairs. In particular,
\[
 \boxed{4\le N_2\le7,\qquad 4\le N_4\le7.}\tag{2}
\]
If `N2=7,N4=4`, the separator word starts and ends with `2`, has no `44`,
and has exactly two `22` adjacencies; both `U1` moves lie in those two
intervals. If `N2=4,N4=7`, the same statement holds with `2,4` and
`U1,U5` interchanged.

## The final-column and p345 crossings

The unique `U6` must flip `p456` from negative to positive. This minor starts
negative, is positive at the all-positive endpoint, and no generator other
than `U6` changes it.

The two `U5` moves **strictly surround** that unique `U6`:
\[
 \boxed{U_5\quad\cdots\quad U_6\quad\cdots\quad U_5.}\tag{3}
\]
This ordering has a noncomputational reduction to the existing `8+m` upper
bound. Any generic six-column trajectory can be padded by a seventh
column close to its fixed sixth column so that its starting minimum
variation is at most 3; a generic perturbation preserves any finite list of
strict target signs. Consequently a 12-change trajectory must contain a
`U6` move, as proved in `minimum_moves_27.md`.

If the unique `U6` came after both `U5` moves, delete it: the first five
columns, and hence every recorded target sign, stay exactly the same,
because the changed sixth column is never used again. If it came before
both `U5` moves, replace the **initial** sixth column by its value just
after `U6` and then delete `U6`. Earlier moves at indices 1–4 cannot use
the sixth column; after the deleted move every subsequent `U5` sees the
same sixth column as before. Thus the entire first-five-column trajectory
and all recorded target signs are again identical. In either case, a small
generic perturbation of the modified six-column input and move times avoids
any newly singular recorded minors while retaining the 12 strict target
changes. This would give a 12-change path with no `U6`, a contradiction.
The same argument shows more generally that if a 12-change trajectory has
exactly one `U6`, at least one `U5` must occur on each side of it.

There is also at least one `U2` strictly between the two `U1` moves. If not,
`U1` commutes with every intervening generator (all are at least `U3`), so
the two positive `U1` shears can be merged into one. This reordering
preserves all target changes because `U1` itself fixes `p123`. Perturb any
new intermediate degeneracies without altering the strict target signs; the
result would be a 27-move 12-change path, contradicting the sharp lower
bound.

**Both `U5` moves must flip `p345`**: the first changes it from positive to
negative, and the second from negative to positive. To see this, `p345`
is unchanged by every generator except `U5`, and varies affinely during a
single `U5` move. With two such moves and positive initial and final signs,
either both flip it or neither does. The latter possibility is excluded by
the following finite sign-prerank statement, which in fact allows unlimited
`U1`, `U2`, `U3` and `U4` moves and does not assume an all-positive endpoint.

For the first five columns, impose the five necessary three-term
Grassmann–Plücker sign relations. The resulting 384 sign states have 1,152
move-labeled edges for `U1,...,U4`. Their future target-sign-change
prerank `R` has maximum 6. In fact, the exact finite graph satisfies
\[
 R\in\{V,V+2\},\qquad
 V=\operatorname{var}(p_{123},p_{124},p_{134},p_{234},p_{345}).\tag{4}
\]
Each value occurs in 192 of the 384 states. This is a **finite sign-graph
identity**; a nonenumerative proof remains open. The independently checked
rank-two cyclic-word formula yields the same identity. Among any two
accepted states that agree on the four minors supported on columns
`1,...,4` **and** on `p345`, the prerank can therefore increase by at most 2. A `U5` move preserves the four first-four minors; if
it also preserves the sign of `p345`, it is contained in such an arbitrary
jump. Therefore at most two such `U5` jumps allow at most `6+2+2=10`
target changes, not 12. `U6` moves are invisible to the first five columns.

`verify_3x5_group_rank_jump.py` reconstructs all 384 states and their 1,152
edges, verifies every prerank inequality, and checks all 4,608 ordered
pairs in the 32 equal-sign groups. An independent rank-two cyclic-word
verification is `verify_rank2_cyclic_jump.py`: after normalizing `p345>0`,
Gale duality gives 192 cyclic words in 16 fixed-four-sign groups, and in
each group their cyclic preranks range over exactly two values differing by
2. This is a finite rank-two verification, not a claimed nonenumerative
formula. Its 16 group ranges are recorded in `verify_rank2_cyclic_jump.log`.
The direct 384-state verifier does not rely on the rank-two theorem.

## Reproduction

```sh
python3 verify_3x5_group_rank_jump.py
python3 verify_rank2_cyclic_jump.py
python3 ../verify_compact_witness.py ../optimal_witness.json
```

The corresponding outputs are saved as `.log` files. For the exact witness,
`../verify_compact_witness.py` directly replays all 28 rational shears and
checks all 35 minors and 17 auxiliary features at all 29 vertices. The
first-five sign graph is an overapproximation, and simultaneous 3-by-5
minor events can be separated; in fact their pairwise root resultants factor
into products of two nonzero invariant minors, so such ties do not occur
inside a shear from a strict 3-by-5 start.
