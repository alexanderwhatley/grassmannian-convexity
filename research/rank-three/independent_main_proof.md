# Rank-three quadratic bound: proof candidate

> **Research status:** This is an unrefereed proof candidate. Its mathematical assumptions and its transfer to continuous formulations need independent expert review.

**Status.** Unrefereed independent audit draft. The planar rank formula and
its monotonicity/payment for positive adjacent shears appear in
Saldanha–Shapiro–Shapiro, *Grassmann convexity and multiplicative Sturm
theory, revisited*, Definition 3.13, Lemma 3.16 and Proposition 3.17,
DOI 10.17323/1609-4514-2021-21-3-613-637. Companion notes give independent
elementary derivations of the needed monotonicity and insertion bounds.
No graph enumeration or SMT result is used. The proof below is conditional
on the correctness of those arguments; it is not the conjectured linear
bound and has not been externally reviewed.

## Precise discrete statement

Let n=q+3, q>=2. Start with arbitrary real 3-by-n columns and apply a finite
word of positive elementary shears U_j(t): v_j -> v_j+t v_(j+1), t>0.
At a finite ordered list of strict target brackets, let F count sign
changes of p123=det(v1,v2,v3), with disjoint brackets counted separately.
Then the argument below gives

    F <= 2q^2 - 3q + 4.                                       (T)

For q=1 the same argument gives F<=3. The conventional discrete sequence
of strictly full-minor-generic vertices is a special case. Brackets may be
selected inside a shear; split the shear at those times. The claimed bound
does not count tangential zeroes of even multiplicity. A transfer to smooth
strictly convex curves' *distinct* zeros uses the separate published
discrete/continuous Lemma 2 in DOI 10.5802/crmath.343.

## Reduction to strict recorded vertices

Keep all positive shear parameters fixed. The prefix at any specified
vertex acts on the initial matrix by right multiplication by an invertible
unipotent matrix. Hence each recorded 3-by-3 minor is a nonzero polynomial
in the initial 3n entries: under the invertible map its value can be made
nonzero by choosing a generic matrix at that vertex. A finite union of
these zero sets has empty interior. Perturb the initial matrix arbitrarily
slightly so every recorded maximal minor is nonzero while preserving every
strict target-bracket sign. Therefore it is enough to prove (T) for strict
vertices. In particular the final three initial columns are independent.

For strict vertices, fix one covector omega positive on those final three
columns. Such covectors form a nonempty open cone. At each recorded vertex
let L_ab=span(v_a,v_b). Strict 3-minors imply that distinct unordered pairs
give distinct two-planes. Require omega to be nonzero on each recorded
column and on each intersection line L_ab∩L_cd. These are finitely many
proper hyperplanes, so a choice exists in the open cone. The final three
scalar values remain positive under all later positive shears.

The latter hyperplane exclusion does more than make recorded charts strict.
On U_j only v_j moves, inside the fixed starting-vertex plane L_(j,j+1).
Every scalar pole is isolated; at that pole its nonzero column spans
L_(j,j+1)∩ker omega. After deleting the singular projected ray, the
retained limiting ray is this line; the remaining rays use distinct
starting-vertex pair planes. Thus their retained limiting list is strict.
The j=3 chart switch has retained planes L12,L34,L45,... on both sides;
the j=2 B chart is regular. The j=1 comparison below is direct.

On a pole-free piece of one physical shear all but at most one selected
projected ray have fixed, pairwise distinct directions. This makes each
projected collision isolated. The exceptional moving ray cannot vanish:
in chart A use the nonzero invariant `p_(j-1,j,j+1)` for U_j, j>=2;
in chart B use `p123` for U1 and U2, `p134` for U3, and
`p_(j-1,j,j+1)` for U_j, j>=4. Each is a recorded maximal minor and is
invariant under its stated shear. More precisely, for U3 in chart B,
`det(r3,z3)=p134/a1`; this case needs a nonconsecutive minor. In chart A
U1 merely rescales a nonzero first ray. At a U3 target zero the nonzero
moving column lies on L12∩L34, so the chosen omega is nonzero on it. Thus
no scalar pole coincides with a target zero. Other U_j preserve p123.

## The planar rank and elementary insertion lemma

For a strict list of m>=2 oriented plane rays, take the m-1 signed turns
of magnitude in (0,pi), group consecutive turns of equal sign into r runs,
and let c be the sum over runs of the number of complete pi-turns. The
published rank is rho=r-1+2c. It is nonnegative and at most 2(m-2), is
unchanged by positive ray scalings, orientation reversal or reversal of the
list, and is nonincreasing under a positive adjacent shear. A sign flip of
the first 2-by-2 minor pays at least one rank unit. Reversing the list
also gives monotonicity for adding a positive multiple of a predecessor.
An independent run-length proof of the monotonicity is given in
[planar_shear_elementary.md](planar_shear_elementary.md); [planar_source_bridge.md](planar_source_bridge.md) spells out the
normalization needed to invoke the published result directly.

**Insertion lemma.** From a strict list of at least three rays, deleting
one ray lowers rho by an integer in [0,4].
If the ray is second or second-to-last the drop is at most 3; at an
endpoint it is at most 2. Reflecting the first ray changes rho by at most 1.

The all-length elementary proof is given in the companion
[insertion_lemma_independent.md](insertion_lemma_independent.md). In particular, the reflection claim
follows immediately from endpoint insertion: both rho values are in
`[rho(tail),rho(tail)+2]`, while reflecting the first ray reverses the
first turn and hence reverses their rank parities. For two rays both ranks
are zero. No physical-path assumption enters this planar lemma.

## Charts and payment away from poles

In coordinates omega=(first coordinate), write v_i=(a_i,b_i). Set

    z_i=b_(i+1)-(a_(i+1)/a_i)b_i,
    r_i=b_i-(a_i/a1)b1,
    A=(z1,z2,z3,...,z_(n-1)),
    B=(r2,r3,z3,...,z_(n-1)).

Use A when a1*a3>0 and B when a1*a3<0. In either chart the determinant
of the first pair equals p123/a1, up to the fixed orientation convention.
For U_j on an interval where a_j keeps sign put c=a_j'/a_j>0. In A,
U1 rescales z1 positively; for j>=2,

    z_(j-1)'=c z_(j-1)+t z_j,      z_j'=z_j/c.

In B, U2 sends r2 -> r2+t r3; U3 and subsequent generators have the
same formula with r3 in the z2 position; and U1 has

    r2'=(a1/a1')r2,   r3'=r3-(t*a3/a1')r2.

The U1 coefficient is positive in B. Hence rho never increases between
scalar poles. A target flip away from poles flips the first planar minor
while a1 keeps sign and costs at least one rank unit.

## Rank jumps at scalar poles

Only indices j<=q can have poles, since the last three scalar coordinates
stay positive. All comparisons below are one-sided strict markers.

* j>=4: delete z_j on each side. The remaining oriented lists have the
  same strict limit. By the insertion lemma both ranks lie in a common
  interval of length 4. Thus the upward jump is <=4.
* j=3: A and B switch. After deleting z3, the first rays are identically
  z1=r2; z2 and r3 both tend to the same pole column v3; later rays are
  unchanged. The strict common-list argument again gives jump <=4.
* j=2 in B: the chart is regular and its limit is strict, so the jump is
  zero. In A, delete z2. Both ranks are in a common interval of length 3.
  The first turn sign equals sign(p123/a1), unchanged by U2, and the last
  turn sign is unchanged. Since q>=2 there is a separate last turn. Rank
  parity is the parity of changes in the turn-sign row, so the jump is even
  and therefore <=2.
* j=1: put c=a3/a2 and compare equal physical-time offsets about the
  pole. The exact identity is r2^+=-r2^- and r3=z2+c r2, with z2 fixed.
  If c>0 the chart goes B -> A: reflect the first ray and then add c times
  that reflected predecessor to the second. If c<0 the chart goes A -> B:
  add (-c) times the predecessor to the second and then reflect the first.
  The positive predecessor addition cannot increase rho and the reflection
  costs at most one; the upward jump is <=1. p123 is unchanged by U1.

## Interface potential and conclusion

Let e_j=[sign(a_j)!=sign(a_(j+1))]. Since a_(q+1),a_(q+2),a_(q+3) are
positive, e_(q+1)=e_(q+2)=0. For q>=2 set

    H=e1+e2+sum_(j=3)^q (4j-7)e_j.

At a j-pole e_j changes 1 -> 0 and only e_(j-1) can toggle. The drop in H
is 1 for j1; at j2 it is 0 if e1=0 (exactly the regular B case) and 2 if
e1=1 (the A case); at j>=3 it is at least 4. These drops pay the upper
rank jumps above. Summing over all scalar-pole and strict target intervals,

    F + rho_end + H_end <= rho_start + H_start.

There are n-1=q+2 projected rays, so rho_start<=2q. Also
`H_start <= 2q^2-5q+4`. Discarding nonnegative endpoint terms proves (T).
For q=1, only a1 can have a pole, H=e1<=1 and rho<=2, yielding F<=3.

## Audit limits

This argument is deliberately stated in the finite discrete setting and
imports a published planar theorem. Its mathematical merit and novelty
require independent expert review. It does not establish F<=3q, says
nothing about direct multiplicity counts, and does not use any conjectural
paid-history graph. Exact finite computations are diagnostics only.
