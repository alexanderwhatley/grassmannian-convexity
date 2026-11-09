# An elementary rank-one theorem and exactness through H=3

This note proves a source-free result in every dimension. Its final numerical
application for `(3,7)` additionally depends on the independently checked
finite-graph upper bound described below.

## The theorem

Let `M` be a real generic `k`-by-`n` matrix, `0<k<n`, with positive last
maximal minor. Write `p_I` for its maximal minors, indexed by increasing
`k`-subsets, and order those subsets componentwise:
$$
 I=(i_1,\ldots,i_k)\le J=(j_1,\ldots,j_k)
 \quad\Longleftrightarrow\quad i_a\le j_a\ \text{ for all }a.
$$
A positive move is right multiplication by `U_j(t)=Id+tE_{j+1,j}`, `t>0`.
Let `r(M)` be the supremum of the numbers of sign changes of `p_[k]` along
finite positive-move sequences with generic recorded vertices.

**Theorem.** The following are equivalent:

1. `r(M)<=1`.
2. The negative `p_I` form a lower order ideal in the above order; equivalently,
   there are no `I<J` with `p_I>0>p_J`.

If these conditions hold, `r(M)=0` when `p_[k]>0` and `r(M)=1` when
`p_[k]<0`. If they fail, `r(M)>=2` for a positive target and `r(M)>=3` for
a negative target. No computer check or external positivity theorem is used
in this all-dimensional result. The last maximal minor remains positive and
unchanged under every allowed move. A negative last minor can be normalized
by one orientation-reversing row operation, which does not change the number
of target sign changes.

### 1. The ideal condition is preserved

A cover `A < B` changes one index `b` in `A` to `b+1`, provided that
`b+1` is not already in `A`. Covers generate the whole subset order. Let
$A^j=(A\setminus\{j\})\cup\{j+1\}$ when `j in A` and `j+1 not in A`;
call such a subset *active*. On applying `U_j(t)`,
$$
 p_A(t)=\begin{cases}
p_A+t\,p_{A^j},& A\text{ active},\\
p_A,& A\text{ inactive}.
\end{cases}\tag{1}
$$
Fix a starting generic matrix satisfying the ideal condition, a cover
`A<B`, and `t>=0` for which the endpoint is generic. We show that an
inversion `p_A(t)>0>p_B(t)` is impossible.

* If neither member is active, the two signs do not change.
* If `A` is active and `B` is inactive, necessarily `B=A^j`. If `p_B<0`
  then the ideal condition gives `p_A<0` and (1) keeps `p_A(t)<0`;
  otherwise `p_B(t)>0`.
* If only `B` is active, then `p_A(t)=p_A`. If this is positive, both
  `p_B` and `p_{B^j}` are positive by `A<B<B^j`, so `p_B(t)>0`.
* Suppose both are active. Set
  $a=p_A,\ b=p_B,\ a'=p_{A^j},\ b'=p_{B^j}$.
  Since `A^j<B^j`, if `b>0` then `b'>0` and `p_B(t)>0`. If `b<0`,
  the ideal condition gives `a<0`. If `a'<0`, then `p_A(t)<0`. The
  remaining case has `a'>0`, hence also `b'>0`.

  Write `S=A∩B` with `j` removed. The cover changes an index `b0` to
  `b0+1`, and both members being active forces the pairs `{j,j+1}` and
  `{b0,b0+1}` to be disjoint. Let `x<x+1<y<y+1` be these four indices,
  put `L=S∪{x,x+1}` and `H=S∪{y,y+1}`, and sort all subsets. The short
  Plücker relation is
  $$
     a b'-b a'=p_Lp_H.\tag{2}
  $$
  Here `L<=A` and `A^j<=H`; hence `p_L<0<p_H`, so the right side is
  strictly negative. At the zero `t_A=-a/a'>0` of `p_A(t)`,
  $$
     p_B(t_A)=\frac{b a'-a b'}{a'}>0.
  $$
  As `b'>0`, `p_B(t)>0` whenever `p_A(t)>=0`. This excludes an
  inversion.

Thus the ideal condition holds after each positive shear at every generic
endpoint, even if intermediate wall events coincide. Moreover, whenever
`p_A>0`, its possible coefficient `p_{A^j}` is also positive; (1) shows
that this minor stays positive after a shear. By iteration, once the target
becomes positive it stays positive at all later generic vertices. Therefore
an ideal-sign start has at most one future target sign change.

### 2. A positive word with all allowed minors positive

Put `r=n-k` and
$$
 L_0=\bigl(U_1(1)U_2(1)\cdots U_{n-1}(1)\bigr)^{kr}.
$$
For any two `k`-subsets `K,T`, its `k`-by-`k` minor satisfies
$$
 \Delta_{K,T}(L_0)>0\quad\text{if }K\ge T,
 \qquad
 \Delta_{K,T}(L_0)=0\quad\text{otherwise}.\tag{3}
$$
The zero assertion follows from lower triangularity. For positivity, move
occupied indices of `K` to `T` one adjacent step at a time without collisions;
at most `kr` steps are needed, and each can be selected from a separate full
sweep. Every elementary factor has nonnegative minors, so Cauchy–Binet has
only nonnegative contributions and this positive history cannot cancel.
The diagonal minors `Delta_{T,T}(L_0)` equal 1.

Given any chosen subset `I`, set `d_i=t` on `I` and `d_i=1` elsewhere,
`D_t=diag(d_1,...,d_n)`, and conjugate the word:
$$
 L(t)=D_t^{-1}L_0D_t.
$$
Every conjugated elementary factor is a positive shear. Cauchy–Binet yields
$$
 p_T(M L(t))=d_T\sum_{|K|=k}
      p_K(M)d_K^{-1}\Delta_{K,T}(L_0),\qquad d_K=\prod_{i\in K}d_i.\tag{4}
$$
For `T=[k]`, (3) makes every coefficient positive and `K=I` is the unique
term of maximal order in `t^{-1}`. Hence its sign dominates at the endpoint
for all sufficiently small positive `t`. This proves the terminal lemma:
**the target can be made to have the sign of any chosen initial maximal
minor**. In particular, a negative target can be made positive using the
last minor; a positive target with some negative minor can be made negative
and subsequently positive, giving at least two changes. Conversely, all
positive minors remain positive by (1).

All recorded vertices can be made generic while preserving the desired strict
endpoint signs. Perturb the finitely many positive parameters of `L(t)`:
for any prescribed vertex and minor, its value is a nonzero polynomial in
the parameters because at the all-zero parameter vector it equals the
corresponding nonzero starting minor. The union of finitely many zero sets has
empty interior. Repeat independently for each concatenated word. This argument
also applies to any finite collection of initially nonzero polynomial features.

### 3. Any inversion forces three changes from a negative target

If the negative minors do not form an ideal, a chain of covers contains a
cover `I<J` with `p_I>0>p_J`. Write
$$
 S=I\cap J,\qquad I=S\cup\{i\},\qquad J=S\cup\{j\},\quad i<j.
$$
Use the word in (3), with a **different** diagonal scaling:
$$
 d_s=t^3\ (s\in S),\qquad d_i=t^2,\qquad d_j=t,\qquad
 d_\ell=1\ (\ell\notin I\cup J).\tag{5}
$$
Among all `k`-subsets, `I` is the unique maximum-weight subset and `J` is
the unique second maximum for the exponent `w(K)` defined by `d_K=t^{w(K)}`.
Indeed their weights are `3(k-1)+2` and `3(k-1)+1`; any subset missing
one of the elements of `S` has weight at most `3(k-2)+2+1`, and any subset
containing `S` and neither `i` nor `j` has smaller weight. This also covers
`k=1`, where `S` is empty.

Apply (4) with (5). At `T=[k]`, the `I` term dominates and is positive.
At `T=J`, the `I` term vanishes because `I` is not componentwise at least
`J`, while `Delta_{J,J}(L_0)=1`; the `J` term is therefore the unique
leading term and is negative. For small positive `t`, the endpoint has
positive target and negative `p_J`, simultaneously. A generic word with the
same strict endpoint signs follows by the perturbation argument above.

If the initial target is negative, reaching this endpoint gives at least one
change. The terminal lemma then makes the target negative using the still
negative `p_J`, and positive again using the invariant positive last minor.
This gives at least three changes. If the initial target is positive, the
terminal lemma already gives at least two. Together with §1, the theorem
follows.

## Application to the certified `(3,7)` graph

Let `H` be the independently audited optimal rank of the accepted 52-sign
formal graph. The graph certificate, together with its real-trajectory
coverage, gives `r(M)<=H(M)` for a generic real starting matrix with all 52
features strict. The following **formal-state** tally is obtained by scanning
all accepted augmented states, their H bytes, target signs, and all 60 covers
of the 3-subset order:

| H | Target sign | Inversion `p_I>0>p_J` for some cover `I<J`? | Formal states |
|---:|:---:|:---:|---:|
| 0 | − | Yes | 25,519 |
| 0 | + | Yes | 25,553 |
| 0 | + | No | 36 |
| 1 | − | No | 7,770 |
| 2 | + | Yes | 83,845 |
| 3 | − | Yes | 631,987 |

The first two rows cannot come from actual strict matrices: the terminal
lemma gives at least one change, contrary to `H=0`. For each of the other
rows, the theorem supplies a lower bound equal to H. Consequently, on every
real generic 52-feature-strict state,
$$
 \boxed{H(M)\le3\quad\Longrightarrow\quad r(M)=H(M).}\tag{6}
$$
This does not assert that every formal state is realizable, nor does it prove
exactness for `H>=4`. The theorem itself is valid without the graph and with
only the 35 maximal minors strict. It is not necessary to assume target
parity from the formal graph; the target sign is checked in the table.

## Reproduction and scope

From this directory:

```sh
c++ -O3 -std=c++17 -pthread check_rank3_sign_ideal.cpp -o check_rank3_sign_ideal
./check_rank3_sign_ideal
python3 check_ideal_shear_lemma.py
```

The first program independently reconstructs the accepted states from the
50-sign generator and algebraic two-G filter, reads
`../grassmann_37_continuation/new_work/augmented_optimal_rank.bin`, and checks
all 60 covers and target signs. It reports counts for all H, including those
above 3. The second program checks 7,637 combinatorial cover/shear cases,
1,204 exact short-Plücker identities, and 4,030 strict shear intervals from
640 integer sign-ideal starts; these finite tests corroborate but do not
replace the all-dimensional proof.

The graph upper-bound dependency is independently checked in
`verify_optimal_rank_bound_independent.cpp/.log`; the reduction from real
generic trajectories to its finite graph is explained in
`audit_50_and_g_coverage.md`. The all-dimensional theorem, §§1–3, does not
rely on either file.
