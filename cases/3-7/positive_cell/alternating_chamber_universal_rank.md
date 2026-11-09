# Every alternating (3,7) start has the same future rank

This is an elementary structural lower-bound argument, conditional only on
the explicit positive-cell chart identities below. It is independent of the
finite sign-graph upper bound. Combined with that bound and the exact
207-move witness, it shows that **every** generic alternating \((3,7)\)
start has future rank exactly 12.

## Statement

Call a generic \(3\times7\) matrix alternating if, after normalizing
\(p_{567}>0\), every ordered maximal minor has sign
\((-1)^{i+j+k}\). Equivalently, after multiplying column \(i\) by
\((-1)^i\), its row space is in the positive Grassmannian. These are the
\(m=4\) starts in the accompanying notes.

**Theorem.** All alternating starts have the same future rank, and that
common value equals the unrestricted \(r(3,7)\). Moreover, from any such
start one can reproduce any chosen finite target-sign sequence of any other
alternating start after 12 positive preparation moves, up to an arbitrarily
small perturbation. In particular, using the saved 207-move witness gives
12 target changes in at most 219 positive moves from every alternating
start. All sampled matrices can be made to have every maximal minor nonzero.

## The positive 12-parameter chart

Let \(U_j(t)=I+tE_{j+1,j}\), so right multiplication adds \(t\) times
column \(j+1\) to column \(j\). Put \(E=[0_{3\times4}\mid I_3]\) and
\[
T(a)=U_4(a_1)U_3(a_2)U_2(a_3)U_1(a_4)
U_5(a_5)U_4(a_6)U_3(a_7)U_2(a_8)
U_6(a_9)U_5(a_{10})U_4(a_{11})U_3(a_{12}).
\]
The row spaces \(\operatorname{row}(ET(a))\), with all \(a_j>0\), are
exactly \(Gr_{>0}(3,7)\), and the chart is a diffeomorphism.

For completeness, normalize any positive representative to \(X=[A\mid I_3]\)
and write its minors as \(p_{ijk}\), so \(p_{567}=1\). Its inverse chart is
\[
\begin{array}{llll}
a_1=\frac{p_{167}p_{456}}{p_{156}},&
a_2=\frac{p_{156}p_{345}}{p_{145}p_{456}},&
a_3=\frac{p_{145}p_{234}}{p_{134}p_{345}},&
a_4=\frac{p_{134}}{p_{234}},\\[3pt]
a_5=\frac{p_{127}p_{156}}{p_{126}p_{167}},&
a_6=\frac{p_{126}p_{145}}{p_{125}p_{156}},&
a_7=\frac{p_{125}p_{134}}{p_{124}p_{145}},&
a_8=\frac{p_{124}}{p_{134}},\\[3pt]
a_9=\frac{p_{126}}{p_{127}},&
a_{10}=\frac{p_{125}}{p_{126}},&
a_{11}=\frac{p_{124}}{p_{125}},&
a_{12}=\frac{p_{123}}{p_{124}}.
\end{array}\tag{1}
\]
All are positive. `audit_positive_cell_elimination.py` independently
verifies (1) by direct inverse-word column elimination and checks the
forward composition. It also checks that the 35 minors of \(ET(a)\) are
positive-coefficient polynomials in the positive parameters. Thus the
inverse exists smoothly and the parameter Jacobian is invertible.

## Transfer between arbitrary alternating starts

Let \(M,N\) be alternating and let
\(D=\operatorname{diag}((-1)^1,\ldots,(-1)^7)\). Then \(X=MD\) and
\(Y=ND\) have positive row spaces. Let
\(\Delta_\epsilon=\operatorname{diag}(\epsilon,\epsilon,\epsilon,\epsilon,1,1,1)\).
Because the last three columns of \(Y\) are independent,
\(\operatorname{row}(Y\Delta_\epsilon)\to\operatorname{row}(E)\) as
\(\epsilon\downarrow0\).

Write \(\operatorname{row}(X)=\operatorname{row}(ET(a^0))\) with all
\(a^0_j>0\). The map taking \((Y',a)\) to
\(\operatorname{row}(Y'T(a))\), in the standard trailing-minor chart, is
smooth near \((E,a^0)\). Its derivative in \(a\) is invertible by (1).
The implicit function theorem therefore gives \(a(\epsilon)>0\) such that
\[
\operatorname{row}(X)
=\operatorname{row}(Y\Delta_\epsilon T(a(\epsilon)))
\]
for all sufficiently small positive \(\epsilon\). Thus for some
\(A\in GL_3\),
\[
M=A N\Delta_\epsilon T(-a(\epsilon)),
\]
where \(T(-a)\) denotes the same word with all parameters negated; this uses
\(DU_j(t)D=U_j(-t)\). Right multiplication by \(T(-a)^{-1}\) is the
**reversed** word, all with positive parameters, and sends \(M\) exactly
to \(AN\Delta_\epsilon\).

For a positive shear \(U_j(t)\) in a fixed trajectory from \(N\), start
from \(AN\Delta_\epsilon\) and use the conjugated shear
\(U_j(t d_j/d_{j+1})\), where \(d_i\) are the diagonal entries of
\(\Delta_\epsilon\). The resulting trajectory is exactly the original
trajectory multiplied on the left by \(A\) and on the right by
\(\Delta_\epsilon\), so all target sign changes agree.

Some exact preparation vertices might lie on a minor wall. Perturb the 12
positive preparation parameters slightly. Each preparation-vertex minor is
a nonzero polynomial in these parameters: at the all-zero parameter vector
it equals the corresponding nonzero minor of \(M\). Hence finitely many
minor walls can be avoided, while the finitely many strict signs of the
chosen target trajectory persist. This proves the rank comparison in both
directions; exact row-space equality after the perturbation is unnecessary.

## Why the common rank is global

Take any generic start \(M\) with \(p_{567}>0\). Let
\(L_0=(U_1(1)\cdots U_6(1))^{12}\) and
\(H_\epsilon=\operatorname{diag}(1,1,1,1,\epsilon,\epsilon,\epsilon)\).
The matrix \(L_\epsilon=H_\epsilon^{-1}L_0H_\epsilon\) is a finite
product of positive shears. Put \(M_0=M L_\epsilon^{-1}\), so the same word
moves \(M_0\) exactly to \(M\).

For any ordered triples \(K\ge I\) componentwise, both \(L_0\) and
\(D L_0^{-1}D\) have strictly positive minors \(\Delta_{K,I}\). One way to
see this is to move the three occupied indices of \(I\) one step at a time
to those of \(K\), rightmost first. This takes at most 12 adjacent particle
moves; twelve full sweeps contain that word as a subword, and every
contribution to the relevant exterior-power coefficient is nonnegative.

Cauchy–Binet expands each \(p_I(M_0)\) as a sum over triples \(K\). The
term \(K=567\) is the unique term with three factors \(\epsilon^{-1}\)
and therefore dominates all other terms as \(\epsilon\to0^+\). Its sign
is
\[
\operatorname{sgn}p_{567}(M)\,(-1)^{5+6+7+\sum_{i\in I}i}
=(-1)^{\sum_{i\in I}i}.
\]
For small enough \(\epsilon\), all 35 minors of \(M_0\) have alternating
signs. Perturb the positive word parameters, redefine its starting matrix
as \(M\) times the exact inverse word, and avoid the finitely many
intermediate minor walls. The target endpoint remains exactly \(M\). Thus
every finite trajectory from any generic \(M\) can be prefixed by one from
an alternating start. The common alternating-start rank is consequently
the unrestricted supremum.

## Reproduction and scope

```sh
python3 audit_positive_cell_elimination.py
python3 audit_positive_cell_transfer_numeric.py
```

The first check uses SymPy for exact symbolic identities. The second is a
nonproof diagnostic: it solves one local chart equation numerically, rounds
the 12 positive parameters to rationals, and then exactly replays all 219
positive shears from a different alternating start. It checks all 35 minors
at every vertex and obtains 12 target changes. The fixed 207-move witness and
its independent rational audit are described in
`maximal_closed_certificate.md`.

This argument does not assert a uniform lower bound of 12 for starts with
\(m<4\), and does not require every path in the finite sign graph to be
realizable.
