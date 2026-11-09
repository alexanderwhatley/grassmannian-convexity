# Positive rectangular charts and maximal future rank in every dimension

This note records a generalization of the positive-cell transfer and an independent audit of the backward-extension argument. The positive Grassmannian facts quoted in §1 are published inputs. The remaining deductions, including the uniform-rank statement, are given in full. No claim of novelty is made.

Fix integers `0 < k < n`, put `r=n-k` and `d=kr`, and write
\[
 U_i(t)=I+tE_{i+1,i},\qquad E=[0_{k\times r}\mid I_k].
\]
Multiplication on the right by `U_i(t)` with `t>0` is the allowed move: it adds `t` times column `i+1` to column `i`. A matrix is *generic* if every maximal minor is nonzero. Its future rank is the supremum of the numbers of sign changes of the minor in columns `1,...,k`, over all finite positive-move sequences whose recorded matrices are generic. The published finiteness theorem shows this supremum is finite; finiteness is not needed for the equalities below.

A generic matrix `M` is *alternating* if after normalizing its last maximal minor `p_{r+1,...,n}` to be positive,
\[
 \operatorname{sgn} p_I(M)=(-1)^{\sum I+\sum_{j=r+1}^n j}.
\]
Equivalently, its row space becomes Plücker positive after multiplying column `i` by `(-1)^i`. By the generic Karp sign-variation formula, these are precisely the starts with minimum sign variation `n-k` (see `universal_first_wall_all_k_n.md`).

**Theorem A.** For every alternating generic start `M`, its future rank equals the supremum over *all* generic starts. Equivalently, any finite generic trajectory from any generic start can be matched, with at least as many target sign changes, by a finite generic trajectory from the prescribed alternating `M`. More precisely, between any two alternating starts there is a `d`-move preparation sending the first to a positive column scaling and left `GL_k` transform of the second; every maximal minor keeps its strict initial sign throughout this preparation.

In particular, the known existence of a generic sequence with `k(n-k)` target sign changes implies that **every** alternating generic start has at least `k(n-k)` future sign changes. No upper bound `k(n-k)` is asserted here for general `(k,n)`; it is exactly the open Grassmann convexity assertion. For `(3,7)`, the separately certified upper bound 12 and exact 207-move witness imply that every alternating start has rank exactly 12, attained in at most `12+207=219` moves.

## 1. A positive rectangular chart

Let the ordered rectangular word be
\[
\begin{split}
 \mathbf i_{k,n}={}&(r,r-1,\ldots,1;\ r+1,r,\ldots,2;\\
 &\ldots;\ r+k-1,r+k-2,\ldots,k),
\end{split}\tag{1}
\]
with `d` entries, and define `T(a)` as the product, in this order, of the corresponding `U_i(a_j)`. Let `Gr_{>0}(k,n)` consist of row spaces admitting a representative whose maximal minors are all strictly positive.

**Theorem B.** The map
\[
 \Phi:(\mathbb R_{>0})^d\longrightarrow Gr_{>0}(k,n),\qquad
 a\longmapsto\operatorname{row}(E T(a)),\tag{2}
\]
is a real-analytic diffeomorphism. In particular it is surjective, and its differential is invertible at every positive parameter vector.

We use these precise external facts:

* Bloch–Karp, [On two notions of total positivity for partial flag varieties](https://arxiv.org/abs/2206.05806), Theorem 1.1 applied to a one-element set and Remark 2.6: every Plücker-positive `k`-plane is spanned by the first `k` columns of a lower unipotent matrix whose minors not forced to vanish by triangularity are strictly positive.
* Berenstein–Fomin–Zelevinsky, [Parametrizations of canonical bases and totally positive matrices](https://math.hse.ru/data/2020/02/04/1571787948/BFZ_TP.pdf), Proposition 1.7: for any reduced word for the longest permutation, the corresponding positive elementary product is a bijection onto the totally positive unipotent group. We apply the transpose form to lower unipotents; reversal of a reduced longest word is again a reduced longest word.

For clarity, the second fact also yields a smooth inverse on the positive group, not merely a set-theoretic bijection. One elementary way to see this is to use the reduced word `(1; 2,1; ...; n-1,...,1)`: the entries of the last row of its positive product are successive products of the parameters in the final block, so those parameters are smooth ratios of positive last-row entries. Remove that block and proceed inductively. Any two reduced words are connected by commuting and braid moves. The braid reparameterization for lower shears is
\[
 U_i(a)U_{i+1}(b)U_i(c)
 =U_{i+1}\!\left(\frac{bc}{a+c}\right)
  U_i(a+c)
  U_{i+1}\!\left(\frac{ab}{a+c}\right),\tag{3}
\]
and the analogous relation with the indices interchanged is its smooth positive inverse. Thus each change of reduced-word parameters is a positive smooth diffeomorphism.

### Proof of surjectivity and the Jacobian claim

Let `R` be the `n`-by-`n` reversal matrix. Reversing all columns of a positive row plane changes all its maximal minors by the same sign `(-1)^{k(k-1)/2}`. Given positive row plane `X`, apply the first quoted fact to the column plane `(X R)^T`: it is spanned by the first `k` columns of a lower unipotent matrix `C` that is positive relative to the lower-unipotent group. Then
\[
 L=R C^T R
\]
is also lower unipotent and positive relative to that group, and `row(E L)=X`.

Let `w_P` reverse the first `r` positions and separately reverse the last `k` positions. Choose any reduced word for `w_P`, using only indices unequal to `r`. Its length is `binom(r,2)+binom(k,2)`. The word (1) is reduced of length `kr`, and appending it to a reduced word for `w_P` gives a reduced word for the longest permutation. Indeed (1) takes the identity list to `(r+1,...,n,1,...,r)`; the concatenated word takes it to `(n,...,1)` and its length is `binom(n,2)`.

By the second quoted fact, factor the above `L=P(b) T(a)` with all parameters positive. Since `P(b)` is block diagonal for the split `r|k`, its bottom `k` rows are an invertible row transformation of `E`. Hence
\[
 X=\operatorname{row}(E P(b)T(a))=\operatorname{row}(E T(a)).
\]
This proves surjectivity. Conversely, for any `a>0`, choose any `b>0`; `P(b)T(a)` is positive relative to the lower-unipotent group. Its bottom `k` by arbitrary `k` minors are strictly positive, and those minors agree with the maximal minors of `E T(a)` up to the fixed positive determinant of the bottom block of `P(b)`. Thus the entire image of (2) is positive.

To prove the Jacobian claim, write a generic lower unipotent matrix in blocks as
\[
 L=\begin{bmatrix}L_{11}&0\\L_{21}&L_{22}\end{bmatrix}.
\]
The coordinates of `row(E L)` in the trailing-minor chart are `L22^{-1} L21`. The projection `L -> row(E L)` is therefore a submersion: varying `L21` spans all `d` coordinates. The full positive reduced-word factorization `(b,a) -> P(b)T(a)` is a diffeomorphism onto an open set of lower unipotents. Its composite with the projection is `Phi(a)`, independent of `b`; consequently the `d` derivatives in `a` have full rank `d` at every positive `a`.

### Proof of global injectivity

Suppose `Phi(a)=Phi(a')`. Since `T(a)` and `T(a')` are lower unipotent, equality of their bottom row spaces implies
\[
 C_0=T(a)T(a')^{-1}
\]
is block diagonal for the split `r|k`: its bottom-left block vanishes and its top-right block vanishes by triangularity.

Choose a fixed block-diagonal lower unipotent `P_*` whose blocks are positive relative to their unipotent groups. Put `D_t=diag(t,t^2,...,t^n)` and `P_t=D_t P_* D_t^{-1}`. Then `P_t` is block positive, and so is `P_t C_0` for all sufficiently large `t`. Indeed,
\[
 P_t C_0=D_t\bigl[P_*(D_t^{-1}C_0D_t)\bigr]D_t^{-1},
\]
where the middle conjugate of `C_0` tends to the identity and strict block positivity is open. Factor both `P_t` and `P_t C_0` along the chosen reduced parabolic word with positive parameters. We obtain two positive factorizations of the same matrix,
\[
 P_t T(a)=(P_t C_0)T(a'),
\]
along the same reduced longest word. Uniqueness in the second quoted fact gives `a=a'`. Combined with surjectivity and the everywhere-invertible differential, the inverse function theorem gives the globally smooth inverse in Theorem B.

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

The main result here is structural; it does **not** replace the independent finite-graph verification of the `(3,7)` upper bound, and it does not establish that every abstract graph edge path is realizable. The two primary-source positivity results are essential external inputs as presented.

For `(3,7)`, `verify_positive_cell_jacobian.py` independently computes the 12-variable chart Jacobian exactly. `trajectory_search_20261108/verify_positive_cell_factorization.py` independently verifies an explicit rational inverse on all positive Plücker data. `python3 grassmann_37_continuation/history_20261108/verify_general_rectangular_chart.py` checks 435 reduced-word pairs for `n<=30`, the exact braid identity, 21 rational chart Jacobians, 240 positive chart minors, 494 backward dominant coefficients, and epsilon dominance on 21 generic matrices. Its output is saved in `verify_general_rectangular_chart.log`. These finite checks illustrate the general proof; they do not replace it.
