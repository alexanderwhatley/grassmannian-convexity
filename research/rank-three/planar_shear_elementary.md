# Direct planar-rank monotonicity under adjacent shears

This independent argument derives the planar input used in the quadratic
audit from the run/carry formula. It is written out to expose the local
wall cases and remains an unrefereed proof candidate. The rank formula
agrees with Definition 3.13 of Saldanha–Shapiro–Shapiro,
https://arxiv.org/pdf/1902.09741; their Proposition 3.17 supplies the
published general monotonicity result, and the normalization bridge is
recorded in `planar_source_bridge.md`.

For a strict list of m>=2 oriented rays, measure shortest turns in units
of pi. A maximal run of equal turn signs has positive length L. Write
`C(L)=1+2 floor(L)`; then `rho+1` is the sum of C over all runs. Any run
length is nonintegral because its endpoint rays are linearly independent.

Consider a positive shear `v_i -> v_i+t v_(i+1)` with the other rays
fixed. Assume first that all fixed rays are pairwise independent and that
the endpoints of the comparison are strict. The moving ray follows the
shortest angular arc toward its successor, without reaching it. Reverse
the orientation of the plane if necessary, so the signed turn

    b = angle(v_i,v_(i+1))

stays in `(0,1)` and decreases. For i>1 write

    a = angle(v_(i-1),v_i) in (-1,1).

Between its crossings of 0 and the endpoints ±1, the quantity a
increases at exactly the rate at which b decreases. Thus a+b is constant.

If a>0, the two turns are in the same positive run. Its length is
constant, as are all other runs, so rho is constant on that piece. If
a<0, the negative run ending with a and the positive run starting with b
both decrease in length; their floor contributions cannot increase, and
all other runs are unchanged. Therefore rho cannot increase between the
two types of local wall `a=0` and `a=±1`. This statement also includes
crossings of other fixed-ray directions: such a crossing merely changes a
run floor, in the nonincreasing direction just described.

At a wall, all relevant fixed-ray directions are distinct. Consequently
the limiting run lengths used below are nonintegral, as are sums that
would otherwise identify their fixed endpoint rays. Let B>0 be the length
of the positive run beginning with b at the wall, and let A>0 be the
length of the immediately preceding run, if there is one. Put
`x=frac(A)`, `y=frac(B)`.

## Wall a=0

Immediately before the wall, a is a small negative turn; immediately
afterward it is positive.

* If the preceding run is negative, the small negative turn joins it on
  the left side and disappears on the right; the two limiting rank sums
  are equal.
* If the preceding run is positive, the local contributions to `rho+1`
  change from `C(A)+C(0+)+C(B)` to `C(A+B)`, where `C(0+)=1`. The drop is

      2 - 2 floor(x+y),

  hence 0 or 2.
* If no preceding turn exists (i=2), the drop is exactly 1.

## Wall a=±1

Immediately before the wall a tends to +1; immediately afterward its
principal value tends to -1. The positive b-turn stays positive.

* If the preceding run is negative, the local contributions change from
  `C(A)+C(1+B)` to `C(A+1)+C(B)`, which are equal.
* If the preceding run is positive, they change from `C(A+1+B)` to
  `C(A)+C(1-)+C(B)`, where `C(1-)=1`. The drop is

      2 floor(x+y),

  hence 0 or 2.
* If no preceding turn exists (i=2), the drop is exactly 1.

Thus every adjacent shear is rank nonincreasing. For i=1 there is no a:
the first run simply shortens as b decreases, so this endpoint generator
is nonincreasing as well. Its first pair determinant never changes sign.
A first-pair determinant sign change under any other generator can occur
only for i=2, and it is exactly one of the two walls above, with rank
drop 1 at that wall. This proves the needed payment.

Reversing the ray order negates and reverses the row of turns, preserving
the run lengths and rho. It converts adding a positive multiple of a
predecessor into an adjacent successor shear, and therefore gives the
same monotonicity for positive predecessor additions.

## Genericity and limits

For one abstract elementary shear starting at a strict list, all fixed
rays are automatically pairwise independent, and the moving ray remains
independent of its donor. A collision of the moving ray with a fixed ray
is a zero of an affine determinant; two such collisions cannot be
simultaneous, since that would make two fixed rays collinear. Thus strict
start and end configurations already suffice for this single-shear
argument; no separate path genericization is needed. A finite word with
non-strict recorded endpoints can instead be perturbed at its original
source while preserving finitely many strict target-bracket signs, as
explained in the main proof. No rank is assigned here at a non-strict
endpoint.
At a physical scalar pole the rank is compared by a separate insertion
argument; this note does not assert monotonicity through a chart pole.
