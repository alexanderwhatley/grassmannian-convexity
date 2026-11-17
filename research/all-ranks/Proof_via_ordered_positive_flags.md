# Grassmann convexity from ordered-positive-flag transversality

## Abstract and scope

Assuming Proposition 4.1 of Qiyue Tang's preprint [T], we prove the sharp Grassmann convexity bound for every $0<k<n$. The proof applies to arbitrary real $k\times n$ starting matrices in the discrete formulation and to arbitrary full-row-rank starting matrices in the strictly positive continuous formulation. In the latter setting it counts zeros with multiplicity, including even-multiplicity zeros.

The mathematical deduction in this note is conditional on [T, Proposition 4.1]; it is not a statement quoted from [T]. Sections 2–6 give the deduction and the discrete proof in full, apart from that explicitly stated input. The final corollary uses the published equivalence with the disconjugate-ODE formulation from [SSS21]. Sharpness is an established result of [SSS21, SSS23]. No independent proof of Tang's proposition or of the published equivalence is claimed here.

## 1. Statements

Fix integers $n\ge2$ and $1\le k<n$, and write $[k]=\{1,\ldots,k\}$. Let $E_{ij}$ denote the matrix with one in position $(i,j)$ and zeros elsewhere.

**Theorem 1 (continuous form, assuming [T, Proposition 4.1]).** Let $J\subseteq\mathbb R$ be a nondegenerate interval, let $V\in\mathbb R^{k\times n}$ have rank $k$, and let $A:J\to\operatorname{Lo}_n^1$ be a smooth curve of lower unitriangular matrices satisfying

$$
A'(t)=A(t)B(t),\qquad
B(t)=\sum_{j=1}^{n-1}\beta_j(t)E_{j+1,j},\qquad
\beta_j(t)>0.
\tag{1}
$$

Set

$$
f(t)=\det\bigl((VA(t))_{[k],[k]}\bigr)
    =\det\bigl((VA(t))_{[:,\,1:k]}\bigr).
\tag{2}
$$

Then every zero of $f$ has finite order and the sum of the orders of all zeros on $J$ is at most $k(n-k)$. At a finite endpoint, the order is understood one-sided; smooth positive coefficients and the flow can be extended slightly beyond that endpoint.

**Theorem 2 (discrete form, assuming [T, Proposition 4.1]).** Start with any real $k\times n$ matrix $V$ and apply an arbitrary finite sequence of moves

$$
v_j\longmapsto v_j+u v_{j+1},\qquad 1\le j<n,\quad u>0,
\tag{3}
$$

with the other columns fixed. At each sampled matrix record the sign of its first $k$-column determinant, delete all zero entries, and count sign changes in the remaining sequence. The number of changes is at most $k(n-k)$. In particular, the $3\times n$ bound is $3(n-3)$.

The discrete theorem includes non-generic and rank-deficient starting matrices. The continuous theorem requires rank $k$: if $\operatorname{rank}V<k$, its determinant is identically zero. Strict positivity in (1) is essential to the asserted finite zero count; weakly nonnegative controls can leave the determinant identically zero on an interval.

**Sharpness.** The bound $k(n-k)$ is attained for every $0<k<n$. For flows, [SSS21, Theorem 2] gives examples with exactly this many simple zeros of the relevant lower-left minor. For discrete moves, [SSS23, Introduction and §3] records examples attaining $k(n-k)$ sign changes.

## 2. The external transversality input

We first specify the convention from [T, Definition 1.1]. A real complete flag $F_\bullet$ is represented by an invertible row frame $G$: $F_j$ is spanned by the first $j$ rows. An increasing tuple $F_\bullet(t_1),\ldots,F_\bullet(t_m)$ is *ordered positive* if it has row representatives $G_i$ for which every left quotient $G_{i+1}G_i^{-1}$ is upper unitriangular and all its minors not forced to vanish by upper triangularity are positive.

We use the following consequence of [T, Proposition 4.1 and Corollary 4.2]:

> **Ordered-positive incidence bound.** Let $F_\bullet^1,\ldots,F_\bullet^m$ be ordered positive real complete flags in $\mathbb R^n$, and let $W$ be a fixed real $p$-plane. For each $i$, take the actual open Schubert cell of $W$ with respect to $F_\bullet^i$, and denote its codimension by $c_i$. Then
>
> $$
> \sum_{i=1}^m c_i\le p(n-p).
> \tag{4}
> $$

Indeed, [T, Proposition 4.1] asserts simultaneous transversality of those actual cells at their common point $W$, for arbitrary codimension sums. The independent conormal spaces embed in the cotangent space to $\operatorname{Gr}(p,n)$, which has dimension $p(n-p)$.

## 3. Ordered positivity of the flow flags

For $s<t$ let

$$
T(s,t)=A(s)^{-1}A(t).
$$

Its lower unitriangularity follows from the equation. More precisely, for every $q$ with $1\le q\le n-1$ and every increasing $q$-tuples $I=(i_1,\ldots,i_q)$ and $J=(j_1,\ldots,j_q)$,

$$
\det T(s,t)_{I,J}>0\quad\text{if }i_a\ge j_a\text{ for all }a;
\tag{5}
$$

the other minors are forced to be zero by triangularity.

For completeness, here is a proof of (5). In the standard ordered wedge basis of $\bigwedge^q\mathbb R^n$, the additive compound $B^{[q]}(t)$ of $B(t)$ has nonnegative entries: $E_{j+1,j}$ replaces $j$ by $j+1$ in a basis index set when $j$ is present and $j+1$ is absent, with coefficient $+1$, and otherwise acts by zero. The matrix of $\bigwedge^q T(s,t)$ satisfies

$$
\frac{d}{dt}\bigl(\bigwedge\nolimits^q T(s,t)\bigr)
 =\bigl(\bigwedge\nolimits^q T(s,t)\bigr)B^{[q]}(t),
 \qquad \bigwedge\nolimits^q T(s,s)=I.
$$

Its time-ordered integral expansion has only nonnegative terms. If $I\ge J$ componentwise, there is a path from $J$ to $I$ by increasing indices one step at a time without collisions: move the rightmost index to its target first, then the next, and so on. For a nonempty path, the corresponding iterated integral is strictly positive because every $\beta_j$ is strictly positive on $(s,t)$; a path of length zero contributes the identity term. It follows that the $(I,J)$ entry, which is the minor in (5), is positive.

Identify column vectors with their transposed rows and put $G(t)=A(t)^T$. Then $G(t)$ represents the flag whose $j$-plane is

$$
F_j(t)=\operatorname{span}\{A(t)e_1,\ldots,A(t)e_j\}.
$$

Furthermore,

$$
G(t)G(s)^{-1}=\bigl(A(s)^{-1}A(t)\bigr)^T=T(s,t)^T.
\tag{6}
$$

By (5), this quotient is strictly totally positive upper unitriangular in the convention of [T]. Thus every finite increasing tuple of flags $F_\bullet(t)$ is ordered positive.

## 4. A local multiplicity–codimension lemma

Fix $t_0\in J$ and set $M=VA(t_0)$. Since $V$ has rank $k$ and $A(t_0)$ is invertible, $M$ has rank $k$. Define its *greedy basis indices* $I=(i_1<\cdots<i_k)$ by taking $i_a$ to be the least index for which the first $i_a$ columns have rank $a$. Put

$$
d(I)=\sum_{a=1}^k(i_a-a).
\tag{7}
$$

Let $\Delta_I(M)$ be the determinant of the columns indexed by $I$.

**Lemma 3.** At $t_0$,

$$
f(t_0+h)=c_I\,\Delta_I(M)h^{d(I)}+O(h^{d(I)+1}),
\qquad c_I>0.
\tag{8}
$$

Consequently $f$ has order $d(I)$ at $t_0$ if $d(I)>0$, and is nonzero there if $d(I)=0$.

*Proof.* Set $T(h)=A(t_0)^{-1}A(t_0+h)$ and $b_j=\beta_j(t_0)>0$. Define $g_1=1$ and $g_r=\prod_{\ell=1}^{r-1}b_\ell$. Repeated integration of the differential equation gives, for $i\ge j$,

$$
T(h)_{ij}
 =\frac{g_i}{g_j}\frac{h^{i-j}}{(i-j)!}
   +O(h^{i-j+1}),
\tag{9}
$$

while $T(h)_{ij}=0$ for $i<j$. The formula also holds one-sided at an endpoint after positive extension of the controls.

For an increasing $k$-tuple $J=(j_1<\cdots<j_k)$, (9) and the determinant expansion imply

$$
\det T(h)_{J,[k]}
 =c_J h^{d(J)}+O(h^{d(J)+1}),
\tag{10}
$$

where

$$
c_J=
\frac{\prod_{a=1}^k g_{j_a}}{\prod_{b=1}^k g_b}
\frac{\prod_{a<b}(j_b-j_a)}
     {\prod_{a=1}^k(j_a-1)!}>0.
\tag{11}
$$

For the last identity, interpret $1/m!=0$ for negative integers $m$ and observe that

$$
\det_{a,b}\!\left(\frac1{(j_a-b)!}\right)
=\frac{\det_{a,b}\bigl((j_a-1)_{b-1}\bigr)}
       {\prod_a(j_a-1)!}
=\frac{\prod_{a<b}(j_b-j_a)}{\prod_a(j_a-1)!}.
\tag{12}
$$

Here $(x)_r=x(x-1)\cdots(x-r+1)$ is a monic polynomial of degree $r$.

Cauchy–Binet gives

$$
f(t_0+h)=\sum_{|J|=k}\Delta_J(M)\det T(h)_{J,[k]}.
\tag{13}
$$

By construction $\Delta_I(M)\ne0$. If $\Delta_J(M)\ne0$, then the first $j_a$ columns of $M$ have rank at least $a$; hence $j_a\ge i_a$ for every $a$. Therefore $d(J)\ge d(I)$, with equality only for $J=I$. In (13), the unique lowest-order nonzero term is the one indexed by $I$. Equations (10)–(11) prove (8). $\square$

**Lemma 4.** Let $W=\ker V$, identified with a row subspace by transpose. The actual Schubert cell containing $W$ relative to $F_\bullet(t_0)$ has codimension $d(I)$.

*Proof.* Let $p=n-k$ and write the complement of $I$ in $[n]$ as $H=(h_1<\cdots<h_p)$. For every $j$,

$$
\dim(W\cap F_j(t_0))
 =j-\operatorname{rank}M_{[:,\,1:j]}
 =\#\{b:h_b\le j\}.
\tag{14}
$$

Thus $H$ is the pivot set of the actual Schubert cell. In a basis adapted to the flag, rightmost-pivot echelon coordinates for that cell have $\sum_{b=1}^p(h_b-b)$ free entries. Its codimension in $\operatorname{Gr}(p,n)$ is therefore

$$
kp-\sum_{b=1}^p(h_b-b)
 =\sum_{a=1}^k(i_a-a)=d(I).
\tag{15}
$$

The equality follows because $I$ and $H$ partition $[n]$. $\square$

## 5. Proof of the continuous bound

Let $t_1<\cdots<t_m$ be any finite collection of distinct zeros of $f$. By §3, the flags at these times are ordered positive. Let $c_i$ be the codimension of the actual cell of the fixed plane $W=\ker V$ at $t_i$. Lemmas 3 and 4 imply

$$
c_i=\operatorname{ord}_{t_i}f\ge1.
$$

The input (4), with $p=n-k$, consequently gives

$$
\sum_{i=1}^{m}\operatorname{ord}_{t_i}f
 =\sum_{i=1}^{m}c_i
 \le k(n-k).
\tag{16}
$$

In particular, there cannot be more than $k(n-k)$ distinct zeros, since any larger finite set would violate (16). All zeros have finite order by Lemma 3, so applying (16) to the finite set of all zeros proves Theorem 1. $\square$

## 6. From flows to arbitrary finite words

Consider a word of $m$ moves (3), with parameters $u_s>0$ and indices $j_s$. The case $m=0$ is immediate. On each interval $(s-1,s)$ choose a nonnegative smooth bump $\phi_s$ of integral $u_s$, with compact support inside that interval. Define the nonnegative controls

$$
\beta_j^0(t)=\sum_{\{s:\,j_s=j\}}\phi_s(t).
$$

Their supports ensure that, starting from $A_0(0)=I$, the solution of (1) satisfies

$$
A_0(s)=\prod_{r=1}^{s}(I+u_rE_{j_r+1,j_r}),
\qquad s=0,\ldots,m.
\tag{17}
$$

Consequently $VA_0(s)$ is exactly the matrix after $s$ moves.

For $\varepsilon>0$ put $\beta_j^\varepsilon=\beta_j^0+\varepsilon$ for every $j$, and let $A_\varepsilon$ start at $I$. These controls are smooth and strictly positive. Continuous dependence of linear ODE solutions on their coefficients gives

$$
A_\varepsilon(s)\longrightarrow A_0(s)
\quad\text{for all }s=0,\ldots,m
\quad\text{as }\varepsilon\downarrow0.
\tag{18}
$$

The sample sequence is finite, so for all sufficiently small $\varepsilon$, every originally nonzero sampled leading minor retains its sign. Each sign change between consecutive nonzero samples forces, by the intermediate value theorem, at least one zero of $t\mapsto\det((VA_\varepsilon(t))_{[:,1:k]})$ in the corresponding interval. These intervals have disjoint interiors. If $\operatorname{rank}V=k$, Theorem 1 bounds their number by $k(n-k)$. If $\operatorname{rank}V<k$, every sampled leading minor is zero, and there are no sign changes. This proves Theorem 2, including the deletion-of-zeros convention. $\square$

## 7. Grassmann-convexity and differential equations

For $V=[\,0_{k\times(n-k)}\ \ I_k\,]$, the function in Theorem 1 is the lower-left $k\times k$ minor of $A(t)$. This is the flag-convex formulation in [SSS21, §1]. That paper states its equivalence, with zeros counted with multiplicity, to the Wronskian formulation of [SSS21, Conjecture 1.2]; see also its discussion of Grassmann curves in §2. Applying the published equivalence gives the following corollary, with the regularity and multiplicity conventions of [SSS21].

**Corollary 5 (assuming [T, Proposition 4.1] and the cited equivalence).** For an order-$n$ real linear ODE with continuous coefficients that is disconjugate on an interval, the Wronskian of any $k$ linearly independent solutions has at most $k(n-k)$ real zeros on that interval, counted with multiplicity.

Theorem 2 also proves the finite-word formulation discussed in [SSS23, Introduction and §2]. The flow proof above directly accounts for multiplicity; it does not infer multiplicity merely by counting sign changes.

## References

- **[T]** Qiyue Tang, *The Shapiro–Shapiro conjecture for ordered positive flags*, Hexagon preprint **hexagon:2611.00409v1** (2026), Definition 1.1, Proposition 4.1 and Corollary 4.2. [Preprint](https://hexagonmath.org/pdf/2611.00409v1).
- **[SSS21]** Nicolau Saldanha, Boris Shapiro and Michael Shapiro, *Grassmann convexity and multiplicative Sturm theory, revisited*, Moscow Mathematical Journal **21** (2021), 613–637, especially §1–2 and Theorem 2. [DOI: 10.17323/1609-4514-2021-21-3-613-637](https://doi.org/10.17323/1609-4514-2021-21-3-613-637); [preprint](https://arxiv.org/abs/1902.09741).
- **[SSS23]** Nicolau Saldanha, Boris Shapiro and Michael Shapiro, *Finiteness of rank for Grassmann convexity*, Comptes Rendus Mathématique **361** (2023), 445–451, especially Introduction and Lemma 2. [DOI: 10.5802/crmath.343](https://doi.org/10.5802/crmath.343).

