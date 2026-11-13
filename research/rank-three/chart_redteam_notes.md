# Red-team notes: chart nonvanishing and simultaneous events

The candidate proof is an unrefereed argument, not an externally verified
theorem. This note makes explicit a corrected dependency in
`independent_main_proof.md` and checks the nearby event logic.

## Which maximal minor prevents a projected ray from vanishing?

Write the fixed nonzero covector as the first coordinate and use
`v_i=(a_i,b_i)` with `b_i` in R^2. The charts are

    z_i=b_(i+1)-(a_(i+1)/a_i)b_i,
    r_i=b_i-(a_i/a1)b1,
    A=(z1,z2,z3,...),  B=(r2,r3,z3,...).

Put `p_abc=det(v_a,v_b,v_c)`. Direct expansion gives, whenever the
denominators are nonzero,

    det(z_(j-1), z_j) = p_(j-1,j,j+1) / a_(j-1),
    det(r2,r3)         = p123 / a1,
    det(r3,z3)         = p134 / a1.

The source is strict, so all minors on the right are initially nonzero.
Each minor below is invariant under the shear being considered:

| Chart and generator | Moving direction | Nonzero invariant |
| --- | --- | --- |
| A, U1 | none: z1 only scales | p123, ensuring z1 is nonzero |
| A, U_j, j>=2 | z_(j-1) toward z_j | p_(j-1,j,j+1) |
| B, U1 | r3 toward r2; r2 scales | p123 |
| B, U2 | r2 toward r3 | p123 |
| B, U3 | r3 toward z3 | **p134** |
| B, U_j, j>=4 | z_(j-1) toward z_j | p_(j-1,j,j+1) |

For example, under U3, `v3 -> v3+t v4`, so
`det(v1,v3+t v4,v4)=det(v1,v3,v4)`. The relevant B-chart determinant is
`det(r3,z3)=p134/a1`, not a determinant controlled solely by p234.
Thus the earlier blanket phrase "consecutive invariant minor" was
incorrect for this case; strictness of every source maximal minor supplies
the required nonconsecutive minor. These determinants prove the moving ray
and its donor remain nonzero and independent on each pole-free piece.
The remaining fixed-direction rays are nonzero by source strictness and
their own invariant pair planes.

## Simultaneous walls and endpoints

Before analyzing a finite word, split shears at all marker times, make all
recorded maximal minors strict by an arbitrarily small perturbation of
the initial matrix, and choose one covector omega which is nonzero on
every recorded column and on each line `L_ab ∩ L_cd` for distinct source
pair planes. Strict 3-minors make all pair planes distinct. These are
finitely many hyperplane exclusions; their union cannot exhaust the open
cone positive on the final three initial columns.

During one elementary shear only one scalar `a_j` moves. All others are
constant and nonzero, and `a_j` is nonzero at both recorded endpoints.
Therefore no two scalar poles occur simultaneously and a scalar pole is
never a recorded endpoint. The scalar is affine and has at most one pole
in that shear. The final three scalar coordinates remain positive along
the whole word.

On each pole-free piece, all rays except at most one have fixed projective
directions, represented by distinct source pair planes. The fixed-plane
sets are, explicitly:

* A: all `L_(i,i+1)` except the plane of the moving ray;
* B under U1 or U3: `L12,L34,L45,...`;
* B under U2: `L13,L34,L45,...`;
* B under U_j, j>=4: `L12,L13,L34,L45,...` with the moving ray's plane
  omitted.

The covector exclusion prevents any pair of those fixed directions from
coinciding. A nonzero moving ray varies by a positive adjacent or
predecessor shear, so each collision with a fixed direction is isolated;
there cannot be two simultaneous collisions with different fixed rays.
For U3 the target minor is affine, `p123(t)=p123(0)+t p124`, and
`p124` is a nonzero source minor. Its target zero is simple, and its
column lies on `L12 ∩ L34`; the covector exclusion makes its scalar
nonzero. It therefore cannot coincide with a scalar pole.

At a U2 pole in chart B, the chart is regular and its limiting first ray
is the line `L23 ∩ ker(omega)`; all other retained pair planes are
distinct from L23, so the limit is strict. For a U2 pole in chart A,
delete z2; the retained list is strict and the first and last turn signs
are unchanged. For q>=2 its last turn is disjoint from the changed
first two rays, including the borderline four-ray case q=2. Rank parity
therefore forces an even jump, and the insertion interval of length 3
gives jump at most 2.

At a U3 pole, delete z3. The A and B retained lists have the same limits:
`z1=r2`, `z2 -> b3`, `r3 -> b3`, and later z-rays are unchanged. Their
limiting pair planes are `L12,L34,L45,...`, so the retained list is
strict. At a U_j pole with j>=4, deleting z_j leaves a common strict
limit in which the preceding retained ray lies in
`L_(j,j+1) ∩ ker(omega)`. In both cases the insertion interval bounds
the upward jump by 4.

At a U1 pole, compare symmetric physical-time offsets about the zero of
a1. With `c=a3/a2` and z2 fixed,

    r2(+)=-r2(-),    r3=z2+c r2.

If c>0 the chart comparison is B-to-A: reflect the first ray, then add
`c` times it to the second. If c<0 it is A-to-B: first add `-c` times
the first ray to the second, then reflect the first. Reflection costs at
most one rank unit; the positive predecessor shear is nonincreasing after
reversing the ray order. The comparisons use strict markers and need no
continuity of a singular chart at the pole.

## Remaining scope

The planar prerank theorem and its elementary insertion proof remain
independent substantive ingredients. The finite-word and smooth-curve
multiplicity bridges use different hypotheses and are documented
separately. None of the finite diagnostics alone proves the candidate.
