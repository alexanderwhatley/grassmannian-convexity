# Positive rectangular charts and maximal future rank in every dimension

This note records a generalization of the positive-cell transfer and an independent audit of the backward-extension argument. The rectangular chart, transfer, and backward extension are proved here without external positivity theorems; the published lower bound is used only for the stated numerical consequence. No claim of novelty is made.

Fix integers `0 < k < n`, put `r=n-k` and `d=kr`, and write
\[
 U_i(t)=I+tE_{i+1,i},\qquad E=[0_{k\times r}\mid I_k].
\]
Multiplication on the right by `U_i(t)` with `t>0` is the allowed move: it adds `t` times column `i+1` to column `i`. A matrix is *generic* if every maximal minor is nonzero. Its future rank is the supremum of the numbers of sign changes of the minor in columns `1,...,k`, over all finite positive-move sequences whose recorded matrices are generic. The published finiteness theorem shows this supremum is finite; finiteness is not needed for the equalities below.

A generic matrix `M` is *alternating* if after normalizing its last maximal minor `p_{r+1,...,n}` to be positive,
\[
 \operatorname{sgn} p_I(M)=(-1)^{\sum I+\sum_{j=r+1}^n j}.
\]
Equivalently, its row space becomes Plücker positive after multiplying column `i` by `(-1)^i`. These are precisely the starts with maximum possible minimum sign variation `n-k`, by the elementary lemma below.

### Elementary characterization of the alternating chamber

For a generic `M`, let `m(M)` be the minimum number of sign changes in the sequence `(ell(v_1),...,ell(v_n))` over linear functionals `ell` that are nonzero on every column. Then `m(M)<=n-k`, and equality holds exactly for the alternating matrices defined above.

For any `(k-1)`-subset `J`, the columns indexed by `J` are linearly independent. Choose a nonzero functional `ell_J` annihilating them. Its evaluations on the other `n-k+1` columns are nonzero. Perturb `ell_J` so the zero entries acquire chosen small signs without introducing further sign changes: this can be done independently because the `J` columns are independent. Hence `m(M)` is at most the variation of the shorter sequence, which is at most `n-k`.

If `m(M)=n-k`, every such shorter sequence must alternate. Write its indices as `i_1<...<i_{n-k+1}`. Up to an overall sign,
\[
 \ell_J(v_{i_s})=(-1)^{|\{j\in J:j>i_s\}|}p_{J\cup\{i_s\}}.
\]
Since `i_s+|{j in J:j>i_s}|` has the same parity as `s+k-1`, the numbers `(-1)^{sum I}p_I` have one common sign among all `k`-subsets `I` containing `J`. The graph of `k`-subsets joined when they share `k-1` indices is connected, so the same is true for every `I`: the minors alternate.

Conversely, suppose the minors alternate and that a strict functional has a sign sequence `sigma_i` with at most `n-k-1` changes. The sequence `tau_i=(-1)^i sigma_i` has at least `k` changes, so it has an alternating subsequence of length `k+1` on indices `i_1<...<i_{k+1}`. In the standard cofactor dependence
\[
 \sum_{a=1}^{k+1}(-1)^a p_{I\setminus\{i_a\}}v_{i_a}=0,
 \qquad I=\{i_1,...,i_{k+1}\},
\]
the coefficient signs are proportional to `(-1)^{a+i_a}`. The chosen alternating subsequence makes `sigma_{i_a}` proportional to exactly the same signs. Applying the functional would therefore give a sum of nonzero terms all with the same sign, a contradiction. Thus `m(M)=n-k`.

**Theorem A.** For every alternating generic start `M`, its future rank equals the supremum over *all* generic starts. Equivalently, any finite generic trajectory from any generic start can be matched, with at least as many target sign changes, by a finite generic trajectory from the prescribed alternating `M`. More precisely, between any two alternating starts there is a `d`-move preparation sending the first to a positive column scaling and left `GL_k` transform of the second; every maximal minor keeps its strict initial sign throughout this preparation.

In particular, the known existence of a generic sequence with `k(n-k)` target sign changes implies that **every** alternating generic start has at least `k(n-k)` future sign changes. No upper bound `k(n-k)` is asserted here for general `(k,n)`; it is exactly the open Grassmann convexity assertion. For `(3,7)`, the separately certified upper bound 12 and exact 28-move witness imply that every alternating start has rank exactly 12, attained in at most `12+28=40` moves. The rational witness is `../optimal_witness.json`, with an all-integer version `../optimal_witness_integer.json`; independent exact replays are `the independent checks in the full archive`.

## 1. An elementary positive rectangular chart

Let the ordered rectangular word be
\[
\begin{split}
 \mathbf i_{k,n}={}&(r,r-1,\ldots,1;\ r+1,r,\ldots,2;\\
 &\ldots;\ r+k-1,r+k-2,\ldots,k),
\end{split}\tag{1}
\]
with `d=kr` entries. Number the blocks by `i=1,...,k`, and the positions in each block by `j=1,...,r`, so the generator at position `(i,j)` is `U_{r+i-j}(t_{i,j})`. Let `T(t)` be their ordered product, and let `Gr_{>0}(k,n)` consist of row spaces admitting a representative whose maximal minors are all strictly positive.

**Theorem B.** The map
\[
 \Phi:(\mathbb R_{>0})^d\longrightarrow Gr_{>0}(k,n),\qquad
 t\longmapsto\operatorname{row}(E T(t)),\tag{2}
\]
is a real-analytic diffeomorphism. The proof below is elementary and does not require a positive-flag completion or a general reduced-word factorization theorem.

For `0<=h<=k`, `0<=ell<=r`, put
\[
 C_{h,\ell}=\Delta_{[1,h]\cup[h+\ell+1,k+\ell]}(E T(t)),\tag{3}
\]
using the empty set when the upper endpoint is smaller than the lower endpoint. In particular, `C_{k,ell}=C_{h,0}=Delta_[k]` for all allowed indices. We claim that
\[
 C_{h,\ell}=\prod_{i\le h,\,j\le r}t_{i,j}
             \prod_{i>h,\,j\le r-\ell}t_{i,j}.\tag{4}
\]
Thus every parameter has the rational inverse
\[
 t_{h+1,r-\ell+1}=
 \frac{C_{h,\ell-1} C_{h+1,\ell}}
      {C_{h,\ell} C_{h+1,\ell-1}},
 \qquad 0\le h<k,\quad 1\le\ell\le r.\tag{5}
\]
These are ratios of homogeneous Plücker coordinates, so they do not depend on the chosen representative of the row plane.

### Proof of (4) by forced particles

Expand `Delta_I(E T)` in the `k`-th exterior power, starting from occupied positions `r+1,...,r+k`. Traversing a factor `U_j(t)` from left to right either leaves the occupied positions unchanged or, if `j+1` is occupied and `j` is empty, moves that particle from `j+1` to `j`, with weight `t`. No particles pass each other. All coefficients are nonnegative: in increasing exterior order an adjacent move introduces no sign. Hence the minor is a sum of nonnegative monomial weights of such histories.

Label particles from left to right by `1,...,k`. Particle `i` starts at `r+i`, and cannot make its first move before block `i`: its first required generator is `r+i-1`, while all earlier blocks have smaller maximum generators. In block `i`, generators decrease consecutively from `r+i-1` to `i`.

For the index set in (3), the final position of particle `i` is `i` for `i<=h` and `i+ell` for `i>h`. If `i<=h`, its last required move is `i+1 -> i`; generator `i` never appears after block `i`. Since the particle cannot begin before block `i`, all `r` of its moves are forced to occur in block `i`.

Put `q=r-ell`. If `q=0`, every particle `i>h` ends at its starting position and never moves. Suppose `q>0`. Particle `k` cannot move before block `k`, so it makes all `q` moves in block `k`. For these descending moves to succeed, particle `k-1` must already be at position at most `k-1+ell` at the start of block `k`: if it occupied any position particle `k` must traverse, its freeing generator would appear later in that same decreasing block. When `k-1>h`, its final position is exactly `k-1+ell`, and positions only decrease, so it is already at its final position. It cannot have begun before block `k-1`; therefore it makes its `q` moves in block `k-1`. Repeat backward for every particle with index greater than `h`. The earlier particles are already forced by the preceding paragraph. This proves there is exactly one history, selecting every generator in blocks `i<=h` and the first `q` generators in blocks `i>h`, giving (4).

The edge cases are included: when `h=k` every block is forced in full; when `ell=r` the later particles never move; and when `h=0,ell=r` the sole history is empty and gives the last minor `1`.

For *any* increasing `I=(I_1,...,I_k)`, one legal history moves particle `i` in its own block from `r+i` to `I_i`, selecting its first `r+i-I_i` generators. This is legal because `i<=I_i<=r+i` and particles `1,...,i-1` have already stopped at `I_1<...<I_{i-1}<I_i`. Its weight is strictly positive; all other histories have nonnegative weight. Consequently **all** maximal minors of `E T(t)` are positive for every `t>0`.

### Surjectivity and smooth inverse

Formula (5) defines a rational map `Psi` from the positive Grassmannian to the positive parameter orthant. Formula (4) gives `Psi(Phi(t))=t` exactly. Work in the affine Grassmannian chart `Delta_[k]!=0`, normalized as `[I_k | B]`; it has dimension `d`. Differentiating that identity shows that the square Jacobian of `Phi` is invertible everywhere for `t>0`, so its image is a nonempty open set of this chart.

Both `Phi` and `Psi` are rational in these affine coordinates. Since `Phi(Psi(X))=X` on that nonempty open set, each coordinate of their difference is the zero rational function on the entire affine chart. For a Plücker-positive `X`, every denominator in (5) is positive and all parameters `Psi(X)` are positive; the output of `Phi` is defined there and has a positive first maximal minor. The rational identity therefore holds at **every** positive `X`. This proves surjectivity and global injectivity, and (5) is a smooth inverse. Only surjectivity and the local Jacobian are needed in the rank argument below.

There is also a direct reconstruction argument, useful as an independent check. In `[I_k | B]`, put `q=min(k-h,ell)`. The coordinate `C_{h,ell}` is, up to a fixed nonzero sign, the contiguous `q`-by-`q` minor of `B` using rows `h+1,...,h+q` and columns `ell-q+1,...,ell`. Process `h` downward and `ell` upward. Expansion in the top-right entry `B_{h+1,ell}` has nonzero cofactor (a smaller minor of `B`, hence up to sign a nonzero maximal minor of the original plane); all other entries in that submatrix are already known. Thus the normalized `C` values recover `B` uniquely. Formula (5) and the boundary `C_{h,0}=C_{k,ell}=Delta_[k]` recursively determine all normalized `C` values from `t`, giving the same surjectivity proof without the rational-identity principle.

The sources [Bloch–Karp](https://arxiv.org/abs/2206.05806), Theorem 1.1 and Remark 2.6, and [Berenstein–Fomin–Zelevinsky](https://math.hse.ru/data/2020/02/04/1571787948/BFZ_TP.pdf), Proposition 1.7, provide an alternative derivation via positive flags and reduced-word factorization; that derivation is preserved in `general_positive_chart_and_global_rank.md`. They are not inputs to the proof above.

## 2. Transfer between alternating starts

Let `M,N` be alternating, and put `D=diag((-1)^1,...,(-1)^n)`. The row spaces `X=row(MD)` and `Z=row(ND)` are positive. Let
\[
 \Lambda_\epsilon=\operatorname{diag}(\epsilon,\ldots,\epsilon,1,\ldots,1),
\]
with `r` copies of `epsilon` followed by `k` copies of `1`, and put `Y_epsilon=Z Lambda_epsilon`. Its last maximal minor is nonzero, so it tends to `row(E)` in the trailing-minor chart as `epsilon -> 0+`.

By Theorem B choose `a^0>0` with `X=row(E T(a^0))`. Apply the implicit function theorem to
\[
 (Y,a)\longmapsto\operatorname{row}(Y T(a))
\]
in the trailing-minor chart near `(row(E),a^0)`. For sufficiently small positive `epsilon` there is `a(epsilon)>0` with
\[
 X=\operatorname{row}(Y_\epsilon T(a(\epsilon))).\tag{4}
\]
Starting at `M`, apply the `d` positive shears of `D T(a(epsilon))^{-1} D`, in reverse word order. This sends `M` to `A N Lambda_epsilon` for some `A in GL_k`: use `D U_i(-t)D=U_i(t)` and the commutativity of the two diagonal matrices.

Every preparation matrix is generic, with the same maximal-minor signs as `M`. Indeed, after alternating reorientation, every intermediate row space is `Y_epsilon` multiplied by a prefix of the positive word `T`; during a partial reverse shear the remaining parameter is still nonnegative. Multiplication of a Plücker-positive plane by nonnegative lower elementary shears preserves strict positivity: Cauchy–Binet has a positive original-minor term and no negative terms. Therefore the preparation has no target sign changes.

A given subsequent move `U_i(t)` from `N` is reproduced from `A N Lambda_epsilon` by `U_i(t lambda_i/lambda_{i+1})`, where `lambda_i>0` are the diagonal entries of `Lambda_epsilon`. Thus every subsequent matrix is the corresponding original one multiplied by `A` on the left and `Lambda_epsilon` on the right. All maximal minors remain nonzero, and target sign changes are identical. This proves the claimed preparation and equality of future ranks for all alternating starts, with no perturbation of the preparation needed.

## 3. Backward extension from an arbitrary start

Let `M` be any generic `k`-by-`n` matrix and normalize its last maximal minor positive. Put
\[
 L_0=(U_1(1)U_2(1)\cdots U_{n-1}(1))^d,\qquad
 H_\epsilon=\operatorname{diag}(1,\ldots,1,\epsilon,\ldots,\epsilon),
\]
with `r` ones and `k` copies of `epsilon`, and define
\[
 L_\epsilon=H_\epsilon^{-1}L_0H_\epsilon,\qquad
 M_0=M L_\epsilon^{-1}.
\]
`L_epsilon` is a word of `d(n-1)` positive elementary shears. For sufficiently small `epsilon>0`, `M_0` is alternating.

To prove this, let `K_*={r+1,...,n}`, let `h_i` be the diagonal entries of `H_epsilon`, and apply Cauchy–Binet:
\[
 p_I(M_0)=\Bigl(\prod_{i\in I}h_i\Bigr)
 \sum_{|K|=k}p_K(M)\Bigl(\prod_{j\in K}h_j^{-1}\Bigr)
 \det((L_0^{-1})_{K,I}).\tag{5}
\]
In the sum, `K=K_*` is the unique term with `epsilon^{-k}`. The prefactor in (5) is positive and may itself depend on `I`; it does not affect the dominant sign. The coefficient of that unique term is nonzero, with sign
\[
 \operatorname{sgn}p_{K_*}(M)\,(-1)^{\sum K_*+\sum I}
 =(-1)^{\sum K_*+\sum I}.
\]
For the nonvanishing, `D L_0^{-1}D` is the reversed word of positive elementary shears, so all its minors are nonnegative. Its minor with row set `K_*` and column set `I` is strictly positive: move the `k` occupied indices of `I` one step at a time to `K_*`, rightmost first. At most `kr=d` adjacent moves are needed, and the successive moves can be chosen as a subword of the `d` complete sweeps. Every term in the exterior-power expansion is nonnegative, so this one positive term suffices. Conjugating by `D` contributes the displayed checkerboard sign. Thus (5) proves strict alternating signs for all `I` when `epsilon` is small enough.

The unperturbed word need not have generic intermediate vertices. Vary its `d(n-1)` positive parameters slightly, and define the initial matrix by `M_0(t)=M L(t)^{-1}` so the final endpoint remains **exactly** `M`. Alternating signs of the first vertex are stable. For any intermediate vertex and any maximal minor, its value is a polynomial in the word parameters; at the all-zero parameter vector it equals the corresponding nonzero minor of `M`. It is therefore a nonzero polynomial, so a sufficiently small perturbation avoids all finitely many zero sets while preserving positivity of the parameters. We have obtained a finite generic positive-move sequence from an alternating `M_0` to the original `M`.

Given any finite generic trajectory starting at `M`, prepend that sequence. Its number of sign changes cannot decrease. By §2 every prescribed alternating start has the same future rank as `M_0`; this establishes Theorem A. The argument also covers a hypothetical infinite supremum, with finite trajectories treated one at a time.

## 4. The case `k=1`

No step above requires `k>=2`. Here `r=n-1`, the rectangular word is simply `(n-1,n-2,...,1)`, and the chart is explicit. If its ordered positive parameters are `(a_1,...,a_{n-1})`, then
\[
 E T(a)=(a_1a_2\cdots a_{n-1},\ a_1a_2\cdots a_{n-2},\ \ldots,\ a_1,\ 1).
\]
Thus for a positive row `(x_1,...,x_n)` normalized by `x_n=1`, the inverse is `a_j=x_{n-j}/x_{n-j+1}`. In particular the Jacobian is nonzero on the positive orthant. The parabolic word has only its first block, and the second `1`-by-`1` block is the identity; the preceding proof applies without alteration. The backward extension uses `d=n-1` full sweeps, and all of its sign arguments remain valid.

There is also an elementary independent check of the rank conclusion: the number of sign changes in the `n` scalar columns is nonincreasing under a positive move, and a sign change in the first column decreases that number by one. Hence every start has rank at most `n-1`. From an alternating start, the word
\[
 (1;\ 2,1;\ 3,2,1;\ \ldots;\ n-1,n-2,\ldots,1)
\]
with each positive parameter large enough to flip the updated scalar makes the first entry change sign once in each block, for `n-1` changes total. All vertices are generic. Thus Theorem A gives the correct exact rank `n-1` in this boundary case.

## 5. Scope and checks

The main result here is structural; it does **not** replace the independent finite-graph verification of the `(3,7)` upper bound, and it does not establish that every abstract graph edge path is realizable. The chart and global-rank equality are proved here without external total-positivity inputs.

For `(3,7)`, `verify_positive_cell_jacobian.py` independently computes the 12-variable chart Jacobian exactly. `trajectory_search_20261108/verify_positive_cell_factorization.py` independently verifies an explicit rational inverse on all positive Plücker data. `python3 grassmann_37_continuation/history_20261108/verify_general_rectangular_chart.py` checks 435 reduced-word pairs for `n<=30`, the exact braid identity, 21 rational chart Jacobians, 240 positive chart minors, 494 backward dominant coefficients, and epsilon dominance on 21 generic matrices. Its output is saved in `verify_general_rectangular_chart.log`. These finite checks illustrate the general proof; they do not replace it.

Additional independent checks: `python3 grassmann_37_continuation/history_20261108/verify_elementary_positive_chart.py` verifies 2,821 unique forced border paths for `n<=14`, legal positive paths to 32,738 minors, cofactor reconstruction in 28 generic matrices, and the exact rational inverse for `n<=9`. Output: `verify_elementary_positive_chart.log`. These are finite checks of the stated conventions, not proofs in all dimensions.
