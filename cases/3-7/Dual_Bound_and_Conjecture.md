# A dual pointwise bound and a general question

**Status:** This note extends the existing computer-assisted candidate; it
is not a separate external audit of the state coverage or perturbation
argument. `verify_dual_bound.cpp` reuses the graph construction in
`verify_upper_bound.cpp`. The general inequality below is a conjecture.

## The dual convention

Let $M$ have full row rank $k$, let $B$ be any full-rank matrix with
row space $\ker M$, and let $J$ reverse the order of the $n$ coordinates.
Set $D=\operatorname{diag}((-1)^1,\ldots,(-1)^n)$ and define the
**positive-move dual** $M^\vee=BJD$. Write $m(M)$ for the minimum sign
variation of a row-space vector with no zero coordinates, and put
$m^\vee(M)=m(M^\vee)$.

The reversal and alternating signs matter. If $M$ is multiplied on the
right by $U_j(t)=I+tE_{j+1,j}$, then its orthogonal complement is
multiplied by $U_j(t)^{-T}$. Conjugating that action by $JD$ gives the
positive move $U_{n-j}(t)$ on $M^\vee$. Complementary-minor duality maps
the target minor of $M$ to the target minor of $M^\vee$ up to a fixed
nonzero factor. This provides a bijection of positive move sequences and
preserves their target sign-change counts. In particular, the claimed
$(3,7)$ upper bound implies the $(4,7)$ upper bound, subject to the same
review of the $(3,7)$ argument. Together with the published cases,
including the reported computer verification of $(3,6)$, this would cover
all $0<k<n\le7$; the current work remains a candidate proof.

The *unreoriented* kernel cannot be substituted for $M^\vee$ in the
pointwise formula. For the included 28-move witness
`simple_moment_28.json`, the columns are
$v_i=(-1)^i(1,i,i^2)$, and there are 12 target changes. The ordinary
kernel of $M$ contains $(1,6,15,20,15,6,1)$: its product with each row of
$M$ vanishes by finite differences. Thus the minimum variation in the
ordinary kernel is zero, and the putative bound $r(M)\le9+0$ is false.
The signed-reversed convention above is the one used below.

## The finite check for (3,7)

For every one of the 23,246,848 admitted 52-sign states, the extended
checker reconstructs the original rank $H$ and directly checks

$$
0\le H\le\min\{8+m,\ 9+m^\vee\}.
$$

It also checks the original transition inequalities and computes the
maximum $H$ for each occupied $(m,m^\vee)$ cell. The following are maxima
in the **abstract over-approximating graph**; they do not alone show that
each maximum is attained by a real matrix or trajectory.

| $m$ | $m^\vee$ | States | Maximum $H$ |
|---:|---:|---:|---:|
| 0 | 0 | 5,863,236 | 8 |
| 0 | 1 | 2,105,499 | 8 |
| 0 | 2 | 22,369 | 8 |
| 1 | 0 | 8,125,049 | 9 |
| 1 | 1 | 4,185,923 | 9 |
| 1 | 2 | 86,863 | 9 |
| 2 | 0 | 1,253,230 | 9 |
| 2 | 1 | 1,471,294 | 10 |
| 2 | 2 | 69,350 | 10 |
| 3 | 0 | 14,229 | 9 |
| 3 | 1 | 30,387 | 10 |
| 3 | 2 | 19,383 | 11 |
| 4 | 3 | 36 | 12 |

All 13 maxima equal the smaller of $8+m$ and $9+m^\vee$.
The computation uses Karp's sign-variation formula on complementary minors
to calculate $m^\vee$. `check_dual_convention.py` separately compares that
formula with direct exact nullspace and circuit calculations for 30 generic
matrices and checks the ordinary-kernel counterexample. Those samples are
diagnostics, not a proof of state coverage.

Reproduce the full table, with a C++17 compiler and approximately the same
memory requirements as the main verifier:

```sh
c++ -std=c++17 -O3 -pthread verify_dual_bound.cpp -o /tmp/verify_dual_bound
/tmp/verify_dual_bound
```

Compile with assertions enabled. The optional convention diagnostic requires
Python 3 and SymPy:

```sh
python3 check_dual_convention.py
```

## General conjecture

For generic full-rank $k\times n$ matrices, the two pointwise estimates
suggest the following question:

$$
r(M)\stackrel{?}{\le}
\min\left\{(k-1)(n-k)+m(M),\ k(n-k-1)+m^\vee(M)\right\}.
$$

The second term is the first applied to the positive-move dual. The
$(3,7)$ finite check above supports this formula, subject to the coverage
argument; it does not prove the conjecture in general. The abstract graph
maxima do not prove that every parameter cell is sharp for real matrices.
A proof of the first term for all $k=3$ would give the conjectured global
bound $3(n-3)$ because $m(M)\le n-3$.
