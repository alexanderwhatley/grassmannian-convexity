# Universal 12-crossing attainment from the alternating chamber

This is a proposed elementary reachability proof, supported by the exact
symbolic checks in `verify_positive_cell_factorization.py`. It should be
reviewed independently, especially the inverse formula and conjugation.
It implies that the ratio ordering previously considered as a possible
additional condition is automatic and is sufficient, because it holds
throughout the alternating Plücker chamber.

Let $D=\operatorname{diag}((-1)^1,\ldots,(-1)^7)$, let
$E=\operatorname{rowspan}[0\mid I_3]$, and put
$L_j(t)=I+tE_{j+1,j}$. Let
$$
 (j_1,\ldots,j_{12})=(4,3,2,1,5,4,3,2,6,5,4,3),\qquad
 T(t)=L_{j_1}(t_1)\cdots L_{j_{12}}(t_{12}).
$$
If $X\in\operatorname{Gr}(3,7)$ has all ordered Plücker coordinates
positive, scale its Plücker coordinates so that $\Delta_{567}=1$, and
set
$$
\begin{array}{ll}
t_1=\Delta_{167}\Delta_{456}/\Delta_{156},&
t_2=\Delta_{156}\Delta_{345}/(\Delta_{145}\Delta_{456}),\\
t_3=\Delta_{145}\Delta_{234}/(\Delta_{134}\Delta_{345}),&
t_4=\Delta_{134}/\Delta_{234},\\
t_5=\Delta_{127}\Delta_{156}/(\Delta_{126}\Delta_{167}),&
t_6=\Delta_{126}\Delta_{145}/(\Delta_{125}\Delta_{156}),\\
t_7=\Delta_{125}\Delta_{134}/(\Delta_{124}\Delta_{145}),&
t_8=\Delta_{124}/\Delta_{134},\\
t_9=\Delta_{126}/\Delta_{127},&t_{10}=\Delta_{125}/\Delta_{126},\\
t_{11}=\Delta_{124}/\Delta_{125},&t_{12}=\Delta_{123}/\Delta_{124}.
\end{array}\tag{1}
$$
All parameters are strictly positive. Direct multiplication and row
normalization verify
$$
 X=E\,T(t).\tag{2}
$$
In detail, `verify_positive_cell_factorization.py` starts with a symbolic
matrix $[A\mid I_3]$, substitutes (1) into the product $T$, and verifies
all 12 rational identities for the entries of the last-three-column
normalization. Conversely, the relevant 12 minors of $ET(t)$ are monomials
in $t$, giving the formulas (1) back. Consequently the map
$t\mapsto ET(t)$ has invertible differential at every positive parameter
point used below. Independently, `verify_positive_cell_jacobian.py` computes
its Jacobian determinant in the same chart to be
$$
t_1^5t_2^4t_3^3t_4^2t_5^4t_6^3t_7^2t_8t_9^3t_{10}^2t_{11}>0.
$$

For completeness, the alternating chamber is exactly the class $m=4$.
First, $m\le4$: take a functional vanishing on columns 1 and 2, then
perturb it so its first two signs equal its third; only four adjacent changes
remain possible. For any generic rank-three matrix $N$, all its ordered
maximal minors have one sign if and only if every generic row $\omega N$
has at most two sign changes. If all minors have one sign, an alternating
signs subsequence on four columns contradicts their cofactor dependence. If
the minors do not all have one sign, two differently signed adjacent triples
lie in some four-set $a<b<c<d$. The image of
$\omega\mapsto(\omega N_a,\omega N_b,\omega N_c,\omega N_d)$ is the
hyperplane orthogonal to
$(\Delta_{bcd},-\Delta_{acd},\Delta_{abd},-\Delta_{abc})$. Because its
coefficients, multiplied by $(+,-,+,-)$, have both signs, that hyperplane
contains a vector with signs $(+,-,+,-)$; a small perturbation avoids
zeros at the remaining columns. Thus some row has at least three changes.
Finally $\operatorname{var}(\omega MD)=6-\operatorname{var}(\omega M)$,
so $m=4$ exactly when $MD$ has same-sign maximal minors.

Now let $M$ be any generic matrix in the alternating chamber
$\operatorname{sgn}p_{ijk}(M)=(-1)^{i+j+k}$ after arranging
$p_{567}(M)>0$. Thus $X=\operatorname{rowspan}(MD)$ has all ordered
Plücker coordinates positive. Let $W$ be one fixed generic matrix with a
known positive 12-crossing path, also in that chamber, and set
$$
 \Lambda_\varepsilon=\operatorname{diag}(\varepsilon,\varepsilon,
 \varepsilon,\varepsilon,1,1,1),\qquad
 Y_\varepsilon=\operatorname{rowspan}(WD\Lambda_\varepsilon).
$$
After normalizing the last three columns, $Y_\varepsilon\to E$ as
$\varepsilon\downarrow0$. Choose the positive $t^0$ in (1) with
$ET(t^0)=X$. The map
$(\varepsilon,t)\mapsto\operatorname{rowspan}(Y_\varepsilon T(t))$
is smooth in the last-three-column chart near $(0,t^0)$, and its
$t$-derivative there is invertible. The implicit-function theorem gives
positive $t(\varepsilon)$ for all sufficiently small $\varepsilon>0$
with
$$
 Y_\varepsilon T(t(\varepsilon))=X.
$$
Multiplying by $T^{-1}$ and then by $D$ gives
$$
 \operatorname{rowspan}\bigl(M\,D T(t(\varepsilon))^{-1}D\bigr)
 =\operatorname{rowspan}(W\Lambda_\varepsilon).
$$
Since $DE_{j+1,j}D=-E_{j+1,j}$, the left factor
$D T^{-1}D$ is the reversed word of **12 positive elementary shears**.
Thus every alternating-chamber $M$ can be moved to a positive column
scaling of $W$, up to an invertible row operation.

The prefix itself remains strictly in the alternating chamber. Choose the
row transformation $A$ such that
$$
 MD=A(WD\Lambda_\varepsilon)L_{j_1}(t_1)\cdots L_{j_{12}}(t_{12}).
$$
It has positive determinant: the $567$ minor is positive on both sides
and is unchanged by these lower shears. After reversing the last $k$
factors and conjugating by $D$, the reoriented matrix is
$$
 A(WD\Lambda_\varepsilon)L_{j_1}(t_1)\cdots L_{j_{12-k}}(t_{12-k}).
$$
The ordered maximal minors of $WD\Lambda_\varepsilon$ are all positive,
and a positive shear replaces each affected minor by the sum of itself and
$t$ times another positive ordered minor. Therefore all ordered minors
stay strictly positive, and the original prefix vertices all have the
alternating signs. In particular no prefix target change occurs.

A witness move $L_j(s)$ from $W$ is transferred to $W\Lambda$ by
using the positive parameter $s\lambda_j/\lambda_{j+1}$. Row operations
and positive column scalings preserve target sign changes, so the 12
witness crossings follow. Thus the 207-move witness produces a total of
**219 positive moves and exactly 12 target changes**. If $M$ and $W$
are rational, positive rational approximations of the 12 prefix parameters
preserve all the finitely many strict vertex signs, so all 219 moves can
be chosen rational.

Therefore every generic $M$ in the alternating chamber has at least 12
target sign changes. Together with the separate global upper certificate,
this gives $r(M)=12$ for every such $M$. It is not merely a claim for
one representative in each abstract sign chamber.
