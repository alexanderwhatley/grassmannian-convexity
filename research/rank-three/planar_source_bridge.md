# Exact bridge to the published planar theorem

The only external planar input used by `independent_main_proof.md` is the
rank monotonicity in Proposition 3.17 of N. Saldanha, B. Shapiro and
M. Shapiro, *Grassmann convexity and multiplicative Sturm theory, revisited*,
DOI 10.17323/1609-4514-2021-21-3-613-637,
https://arxiv.org/pdf/1902.09741.  Definition 3.13 gives its rank as
`2 Cont + s - 1`; in angular terms this is exactly the rank in the audit
draft: each maximal same-sign sequence of shortest turns is one of the
paper's maximal monotone arcs, and its content counts complete half-turns.
Lemma 3.16 gives the parity and Proposition 3.17 the monotonicity of each
generic crossing and strict decrease when the first pair's sign changes.

The paper uses the normalized space `Lo^1_(2 x m)`, with its last two
columns forming the standard basis, and hence positive last determinant.
For an arbitrary strict ray list, choose a fixed invertible plane map taking
its last two columns to that basis.  An invertible linear map acts on the
circle of rays by a homeomorphism commuting with the antipodal map. It
either preserves or reverses every turn sign. On a monotone arc it
preserves the number of complete half-turns, so the stated rank is unchanged.
The fixed map also commutes with each operation of adding a positive
multiple of the next column. The last two columns are either fixed or the
penultimate is sheared by the last, so their determinant stays nonzero;
for a single-shear comparison we normalize at the starting endpoint.
For the paper's literal normalized cross-section, a harmless continuously
varying invertible row map restores the standard basis; it changes neither
the cyclic word nor the positive-shear path's admissible crossing direction.

There is a minor textual point: Definition 3.9 as printed says `k > 2`,
while Remark 3.15 and the proof of Proposition 3.17 discuss motion of the
first ray toward the second. The endpoint motion can be proved directly,
without relying on how that index restriction is read. Hold the tail
`v2,...,vm` fixed and move `v1` toward `v2` by `v1+t v2`, t increasing.
The shortest signed turn from `v1` to `v2` has a fixed sign and decreasing
absolute length `d in (0,1)`, measured in half-turns. If the following
turn has opposite sign, the rank contribution of the first turn over the
tail is constantly 1. If it has the same sign and its initial run has
length B, the contribution is

    2 (floor(B+d)-floor(B)),

which cannot increase as d decreases. When the tail has only one ray,
both ranks are 0. The first determinant cannot change sign under this
motion. Thus the endpoint generator has exactly the needed properties.

For moving `v_i` toward its successor for `i>=2`, the motion is an
admissible generic wall crossing in the sense used in the proof of
Proposition 3.17; that proof explicitly treats `i=2` crossing ray 1 and
its antipode. If a continuous interval contains several walls, apply the
inequality one wall at a time. In the application the single-covector
choice makes each wall isolated, and multiple fixed rays cannot coincide.

Reversing the column order reverses all turns and their order, preserving
the rank. Consequently adding a positive multiple of a predecessor is
also rank nonincreasing. This is how the `U1` chart transition is handled;
it is not a separate assertion from the paper about its preferred matrix
normalization.
