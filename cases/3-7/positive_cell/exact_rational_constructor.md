# Exact terminating construction from any rational m=4 start

`exact_ift_construct_rational.py` takes any generic rational 3-by-7 matrix
with $m=4$, together with an exact $N$-move 12-crossing witness, and
constructs $N+12$ strictly positive rational moves with exactly 12 target
sign changes. It uses Python's standard library and exact rational arithmetic
only. It does not use floating-point optimization, the finite upper-bound
rank, or any precomputed sign-graph data.

Example, from this directory with an independently checked witness:

```sh
python3 exact_ift_construct_rational.py matrix_simple.json \
  ../optimal_witness.json result.json
python3 verify_constructed_12.py result.json
```

The input is either a JSON list of seven 3-component columns or an object
with an `initial_columns` or `columns` field. Integers and rational strings
are accepted. If $p_{567}<0$, the constructor temporarily reflects the
first row; it writes the original input and the same move parameters, since
invertible row operations do not affect target sign changes. It rejects a
nongeneric or nonalternating input. The witness may have any positive move
count: the program checks its initial alternating signs, positivity of every
parameter, strictness of all 35 minors at each vertex, and exactly 12 target
changes. The earlier 207-move witness remains valid input. Assertions must
remain enabled; the program explicitly refuses Python `-O`.
The independently checked 28-move witness gives at most **40 moves** for
every rational generic (m=4) input. The same existence proof applies to
real inputs; rationality is needed only for this exact-arithmetic program.

## Exact cubic construction

Use the notation and factorization in
`positive_cell_universal_reachability.md`: $D=\operatorname{diag}((-1)^i)$,
$E=[0\mid I_3]$, and $T(t)$ is the 12-shear product with word
$(4,3,2,1,5,4,3,2,6,5,4,3)$. Normalize $MD=[X\mid I_3]$ and
$WD=[A_W\mid I_3]$, with all ordered maximal minors positive. The
explicit positive Plücker-ratio formulas give the rational $t^0$ with
$ET(t^0)=[X\mid I_3]$.

For $Y_\varepsilon=[\varepsilon A_W\mid I_3]$, write
$Z(\varepsilon,t)=Y_\varepsilon T(t)=[Z_L\mid Z_R]$, and define the
12 rational-polynomial equations
$$
G(\varepsilon,t)=Z_L-Z_R X=0.\tag{1}
$$
At $(0,t^0)$, $G=0$, and its $t$-Jacobian $J_0$ is invertible by
the positive-cell factorization. Each equation is separately affine in each
$t_i$, so the program computes each column of $J_0$ by the **exact**
difference $G(0,t^0+e_i)-G(0,t^0)$. It inverts $J_0$ over the rationals
and successively solves for rational coefficients $t^{(1)},t^{(2)},t^{(3)}$
such that, for
$$
\widetilde t(\varepsilon)=t^0+t^{(1)}\varepsilon+
 t^{(2)}\varepsilon^2+t^{(3)}\varepsilon^3,
$$
the coefficients through degree three of $G(\varepsilon,\widetilde t)$
vanish. The implicit-function theorem then gives
$\widetilde t(\varepsilon)-t(\varepsilon)=O(\varepsilon^4)$, where
$t(\varepsilon)>0$ is the exact positive solution for sufficiently small
$\varepsilon>0$.

The program tries $\varepsilon=2^{-1},2^{-2},2^{-4},2^{-8},\ldots$.
For each trial it rejects nonpositive parameters, applies the 12 reversed
conjugate shears, and checks that all 35 maximal-minor signs at every prefix
vertex remain the original alternating signs. It then appends the fixed
$N$-move witness, multiplying the parameter of each generator-4 move by
$\varepsilon$, and compares all 35 minor signs at every vertex with the
exact witness template. A trial is returned only after all signs are strict,
the prefix has no target change, and the witness contributes exactly 12.

This search terminates. For the ideal solution, each prefix is
$Y_\varepsilon$ times a prefix of the positive product $T(t(\varepsilon))$.
The Cauchy–Binet summand indexed by the same column set gives
$\Delta_I\ge c_I\varepsilon^{|I\cap[4]|}\ge c_I\varepsilon^3>0$
for $0<\varepsilon\le1$; all other summands are nonnegative. The cubic
parameter error changes each prefix minor by $O(\varepsilon^4)$, so the
prefix signs eventually agree. At the endpoint, dividing the first four
columns by $\varepsilon$ turns the error into $O(\varepsilon^3)$;
the resulting matrix tends to the fixed witness up to an invertible row
operation. Every one of the finitely many strict witness sign conditions
therefore eventually agrees as well. Thus some rational trial succeeds.
This is an existence/termination argument, not an explicit uniform upper
bound on the required exponent or on running time for arbitrary inputs.

## Bit sizes and checks

If every rational Taylor coefficient has numerator and denominator at most
$B$ bits and a successful trial is $\varepsilon=2^{-p}$, a common
denominator bound gives each cubic prefix parameter numerator/denominator
at most $4B+3p+2$ bits. If the witness parameters have at most $B_W$
bits, scaling generator 4 adds at most $p$ denominator bits. Intermediate
rational sizes and the successful $p$ depend on the input; no uniform
polynomial-time complexity claim is made.

The simple moment input and 20 strictly convex, non-conic integer inputs
were constructed and independently replayed. Their sampled curvature ratios
range from about 2.4 to 1,000,000; all succeeded at $p=64$, and the largest
output parameter numerator/denominator was 588 bits. One additional test used
an exact rational row transformation and positive rational column scalings;
another used a row reflection with $p_{567}<0$. These finite tests check the
implementation only; universal termination follows from the argument above.
`verify_constructed_12.py` independently evolves all 35 maximal minors by
exterior-power updates and checks the 12 changes and strictness at all
recorded vertices. The numerical bit-size/timing samples above used the
earlier 207-move witness; the constructor also accepts shorter ones.
