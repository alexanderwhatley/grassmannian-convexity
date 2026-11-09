# A general pencil-variation lower bound

Let `M` be a generic real `k`-by-`n` matrix, with positive last maximal
minor, and let `r(M)` be its future target-sign-change rank under positive
shears, as in `rank_one_ideal_theorem.md`. For a fixed `(k-1)`-subset `S`,
list the complement as `i_1<...<i_q`, `q=n-k+1`, and put
\(I_h=S\cup\{i_h\}\), sorted. Let `var` count sign changes in a list of
nonzero real numbers. Then
\[
 \boxed{\displaystyle r(M)\ge
 \operatorname{var}\bigl(p_{[k]},p_{I_1},\ldots,p_{I_q},p_{[n-k+1,n]}\bigr).}\tag{1}
\]
The result holds for every choice of `S`; maximize (1) over them. The final
minor is the positive invariant last minor. For `k=1`, this recovers the
usual sign-variation lower bound for a row vector. The proof is source-free
and makes no claim that an abstract graph path is realizable.

## Proof

Let `L_0` be the positive word from §2 of `rank_one_ideal_theorem.md`, with
`Delta_{K,T}(L_0)>0` for exactly the componentwise inequalities `K>=T`.
Given a generic current matrix and `h=1,...,q`, we construct a finite
positive word with the following strict endpoint properties:

* the target sign is the current sign of `p_{I_h}`;
* for every `l>=h`, the sign of `p_{I_l}` is unchanged.

Choose integer weights on the column labels by
\[
 w_s=W\quad(s\in S),\qquad
 w_{i_l}=\begin{cases}0,&l<h,\\q-l+1,&l\ge h,\end{cases}
 \qquad W>q(q+1)/2.\tag{2}
\]
Set `d_i=t^{w_i}`, `D_t=diag(d_i)`, and conjugate the word to
`D_t^{-1}L_0D_t`. It is again a positive word. In Cauchy–Binet, the
relative size of the term indexed by `K` is controlled by
`t^{-sum_{i in K} w_i}`. Since `W` exceeds the sum of all weights outside
`S`, every `k`-subset containing `S` has larger weight than every one that
does not. Of the subsets `S∪{i_l}`, the unique largest is `I_h`; thus the
`I_h` term dominates the target and has strictly positive coefficient.

For `T=I_l` with `l>=h`, a subset `I_j=S∪{i_j}` satisfies `I_j>=I_l`
exactly when `j>=l`. Every subset `K>=I_l` not containing `S` has smaller
weight than `I_l`, while among those containing `S`, the unique largest
allowed weight is `I_l` itself. Since `Delta_{I_l,I_l}(L_0)=1`, its initial
minor is the unique leading term of the transformed `p_{I_l}`. Thus all
properties hold for sufficiently small positive `t`. Perturb the finite
word parameters slightly to make all recorded maximal minors nonzero; each
relevant vertex-minor polynomial equals a nonzero starting minor when all
parameters vanish, and the prescribed strict endpoint signs survive.

Starting from `M`, perform this construction for `h=1,2,...,q`. By induction,
when stage `h` begins, each later pencil minor still has its **original**
sign; at its endpoint, the target has the original sign of `p_{I_h}`.
Consequently the successive stages give at least the sign changes in
`(p_[k], p_{I_1},...,p_{I_q})`. Finally use the invariant positive last minor
to make the target positive, contributing any final sign change in (1).
The intermediate word can contain additional flips, which only strengthen
the lower bound.

## An exact integer example

For the matrix with columns
\[
 (0,0,1),(5,25,1),(-3,-7,1),(-2,-5,1),(-7,0,1),(5,-24,1),(4,-10,1),
\]
all maximal minors are nonzero, the target is `40`, and the last minor is
`144`. For `S={2,6}` in one-based indexing, the five pencil minors in
increasing omitted index are
\[
 -245,\quad392,\quad343,\quad588,\quad-49.
\]
Thus the displayed seven-sign list is `+ − + + + − +`, with four changes,
so (1) proves at least four future target changes. The independent formal
52-sign graph gives upper bound `H=4`; therefore its true rank is exactly 4.
The same conclusion has a separate explicit 13-move rational witness in
`projected_rank8_actual_rank4.md`.

## Exhaustive `(3,7)` tally

The following table counts *formal augmented states*. `P` is the maximum of
(1) over all `C(7,2)=21` pencils. Any formal state with `P>H` is unrealizable,
because (1) is a real-matrix lower bound and the certified graph gives upper
bound H. If `P=H`, every actual generic strict matrix in that state has true
rank H. Other entries leave a gap; no assumption of realizability is made.

| H | P values and formal-state counts |
|---:|:---|
| 0 | 0: 36; 2: 222; 3: 3,389; 4: 16,956; 5: 22,130; 6: 8,375 |
| 1 | 1: 7,770 |
| 2 | 2: 83,845 |
| 3 | 3: 631,763; 5: 224 |
| 4 | 2: 122,729; 4: 2,318,427; 6: 48 |
| 5 | 3: 1,657,164; 5: 3,529,923 |
| 6 | 2: 6,682; 4: 5,170,834; 6: 1,421,258 |
| 7 | 3: 164,447; 5: 4,998,915 |
| 8 | 4: 903,488; 6: 1,503,764 |
| 9 | 5: 605,474 |
| 10 | 4: 19,107; 6: 47,641 |
| 11 | 5: 2,225 |
| 12 | 4: 12 |

The graph upper bound and strict-real-trajectory coverage are the independently
checked dependencies described in `audit_50_and_g_coverage.md` and
`verify_optimal_rank_bound_independent.cpp/.log`. The pencil inequality itself
needs neither of them.

From this directory, reproduce the tally with:

```sh
c++ -O3 -std=c++17 -pthread analyze_pencil_lower_bound.cpp -o analyze_pencil_lower_bound
./analyze_pencil_lower_bound
```

The program reconstructs the accepted 50+2 states and reads the saved
`augmented_optimal_rank.bin` byte certificate. Its exact output is retained in
`analyze_pencil_lower_bound.log`.
