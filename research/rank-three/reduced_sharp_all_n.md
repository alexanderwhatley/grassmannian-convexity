# Sharp target-flip bound for reduced shear words: proof candidate

> **Research status:** The all-width construction is an unrefereed proof candidate. Its sorting-network input is published, and finite integer replays are diagnostic rather than a proof for all widths.

## Statement and scope

Let `n=q+3>=4`. Consider physical shears `U_j(t): v_j -> v_j+t v_(j+1)`
with positive times. Call the word *Coxeter-reduced* when its generator
indices, interpreted as adjacent transpositions of `1,...,n`, form a
reduced word. For all real `3 x n` source matrices and all such words,

    F <= M(3,n) = ceil(11n/6)-5 = ceil((11q+3)/6).          (1)

Here `F` counts strictly witnessed changes of sign of `det(v1,v2,v3)`.
For every `n>=4`, equality is attained by a fully generic rational source
and a positive word that is reduced for the longest permutation `w0`.
Moreover, **all times in the word can be taken equal**. A construction
of the source and whole word is given below.

The upper bound is a direct application of the exact sorting-network
theorem of Gaetz, Gao, Jiradilok, Nenashev and Postnikov, *Repeatable
patterns and the maximum multiplicity of a generator in a reduced word*,
arXiv:2204.03033v2, July 2024, Theorem 3.23:
https://arxiv.org/abs/2204.03033 . Every target flip uses a `U3`, and any
reduced word extends to a reduced word for `w0`. This is the **corrected
2024 formula with a ceiling**; an earlier conference version states a
different rounding.

The construction and the flag-minor lemma below establish the matching
physical realization. This is a strict class of words and says nothing
by itself about the unrestricted upper bound `F<=3q`. In particular,
for `q>=2`, the maximum in (1) is strictly less than `3q`, whereas the
known unrestricted lower construction attains `3q`. Thus one cannot
replace every unrestricted history by a single reduced history while
preserving its target flip count, even allowing the source to change.

## A flag-minor dominance lemma

The following holds for arbitrary rank `k`. Let

    w = s_(i1) ... s_(im)

be reduced and `A_w(t)=(I+t E_(i1+1,i1)) ... (I+t E_(im+1,im))`.
Write `W` for the permutation obtained by applying these transpositions to
the initially ordered list, and `W[k]` for the **sorted set** of its first
`k` values. If `V` is a real `k x n` matrix with `Delta_(W[k])(V)!=0`, then
as `t -> +infinity`,

    sign Delta_[k](V A_w(t)) = sign Delta_(W[k])(V).         (2)

The same single sufficiently large `t` works simultaneously for every
prefix of any fixed reduced word whose corresponding source minors are
nonzero. Here and below all row and column minors use increasing indices.

**Proof.** Each elementary factor is totally nonnegative. In the
Cauchy--Binet expansion, a nonzero elementary `k`-minor either leaves a
`k`-set unchanged or replaces `j` by `j+1` when `j` is present and `j+1`
is absent. Call this latter greedy operation `G_j`; it is monotone for
componentwise order on sorted `k`-sets. Processing the factors from right
to left, it follows inductively that every reachable row set lies
componentwise below the greedy reachable set, which itself has positive
minor.

For a reduced suffix `v`, prefixing `s_j` increases Coxeter length, hence
`v^(-1)(j)<v^(-1)(j+1)`. Consequently, if `j+1` belongs to `v([k])`, then
`j` belongs too. Thus `G_j(v([k]))=s_j v([k])`. Induction through the
reversed word identifies the greedy set as `W[k]`. This proves that the
unique coordinatewise maximal `I` with `Delta_(I,[k])(A_w(1))>0` is
`I=W[k]`.

Let `D(t)=diag(t,t^2,...,t^n)`. Since
`A_w(t)=D(t) A_w(1) D(t)^(-1)`, every flag minor is a nonnegative
constant times `t^(sum(I)-k(k+1)/2)`. A strictly smaller componentwise set
has strictly smaller index sum. The `W[k]` term in

    Delta_[k](V A_w(t)) = sum_I Delta_I(V) Delta_(I,[k])(A_w(t))

therefore has strictly largest degree and a positive coefficient. This
proves (2). There are finitely many prefixes, so one t works for all. ∎

If the source is fully generic, every maximal minor of every word vertex
is a polynomial in `t` with its nonzero source minor as constant term.
There are only finitely many forbidden roots. Thus t can be chosen to be
a positive integer, with every vertex fully generic. At `U3`, the target
determinant is affine in the shear parameter, so opposite endpoint signs
give exactly one crossing on that segment; every other generator leaves
the target unchanged.

## Optimal paths and their alternating source minors

The paper's Section 3.5 gives an eleven-edge path `P` of triples,

    123, 124, 125, 145, 245, 345, 456, 457, 567, 578, 678, 789.

It repeats with an index shift by 6. Its six terminal paths `P_r` are:

| r | triples |
|---|---|
| 4 | 123, 124, 134, 234 |
| 5 | 123, 124, 125, 145, 245, 345 |
| 6 | 123, 124, 125, 145, 245, 345, 456 |
| 7 | 123, 124, 125, 145, 245, 345, 456, 457, 567 |
| 8 | 123, 124, 125, 145, 245, 345, 456, 457, 567, 578, 678 |
| 9 | 123, 124, 125, 145, 245, 345, 456, 457, 567, 578, 579, 589, 789 |

For a given `n`, write `n=6m+r`, `4<=r<=9`. Concatenate `m` copies of
`P` shifted by `0,6,...,6(m-1)`, followed by `P_r` shifted by `6m`;
identify adjacent matching endpoint triples. As proved in the cited
paper, this is a monotone weakly separated path with exactly
`ceil(11n/6)-5` edges.

We exhibit a matrix whose minors indexed by this particular path have
alternating signs. Start with six column vectors

    b1=(1,0,0),   b2=(0,1,0),  b3=(0,0,1),
    b4=(-1,-1,-1), b5=(-1,0,1), b6=(-2,-1,-1),

and `D=diag(1,-1,1)`, which has determinant `-1`. For each block, set

    b_(6h+i) = D^h b_i,        1<=i<=6.                       (3)

If `r=9`, replace the very last column, `b_n`, by
`D^m (0,1,-1)`. For the displayed path `P`, the twelve minors using
the unperturbed block `b1,...,b9` are, in order,

    1, -1, 1, -1, 2, -1, 1, -1, 1, -1, 1, -1.             (4)

The `P4` values are `1,-1,1,-1`; the paths `P5,...,P8` are prefixes
of (4). For `P9`, the final-column replacement instead gives

    1, -1, 1, -1, 2, -1, 1, -1, 1, -1, 1, -1, 1.        (5)

Each full `P` block has eleven edges, while applying `D` to all three
columns negates the determinant. Therefore the signs alternate across
every repeated block and across every terminal path, in all six cases.
The `r=9` change affects only the terminal block, not earlier repeats.

To make **all** maximal minors nonzero without changing those signs, let

    v_i(R) = R b_i + (1,i,i^2).                              (6)

For each fixed triple, its determinant is a polynomial in `R` whose
constant term is the positive Vandermonde determinant, so it is not the
zero polynomial. On the displayed path its leading term is the nonzero
determinant of the corresponding `b_i`. Hence every sufficiently large
integer `R`, outside a finite set of roots, gives a fully generic integer
source and retains the alternating path signs. The script finds such an R
by starting at `100n^2` and incrementing it when necessary.

## Explicit reduced full-word lift

The paper's Proposition 2.3 gives a reduced-word lift of every such path.
For reproducibility, here is its deterministic form. Start with the list
`[1,2,...,n]`. At each path transition `A -> B`, write
`x` for the element removed and `y` for the element added. If their
current 1-based positions are `a<=3<b`, append the transpositions

    s_a s_(a+1) ... s_2,   s_(b-1) s_(b-2) ... s_4,   s_3,

omitting an empty product. Apply each swap to the current list as it is
appended. The first group moves `x` to position 3 and the second moves
`y` to position 4. The cited proposition proves every swap is ascending,
so the resulting word is reduced. The `s3` swaps occur exactly once per
path transition and its successive first-three sets are the listed
triples.

At the end the first three entries are the three largest labels. Repeatedly
swap any ascending adjacent pair wholly within positions `1,2,3` or wholly
within `4,...,n`. These additional swaps are ascending, use no `s3`, and
finish at `w0`; the complete word is reduced and has `n(n-1)/2` factors.
Choose `R` by (6), and then choose one large enough common integer time `t`
by the lemma. The determinant signs alternate at all successive `s3`
positions and remain constant at all other positions. This gives exactly
`ceil(11n/6)-5` strict target crossings, proving the matching lower bound.

## Independent finite checks

Run:

    python3 research/rank-three/reduced_sharp_all_n.py 18

The checker independently builds the source and the complete `w0` word,
verifies every ascending swap, every source and vertex maximal minor, and
every sign change, using integer arithmetic and one common time. It is
diagnostic; the path theorem and all-n argument above provide the proof.

