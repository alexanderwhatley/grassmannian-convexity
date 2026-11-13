# Elementary proof of planar one-ray insertion bounds

Angles are measured in units of pi. Let C(x)=1+2 floor(x) for x>0.
For a strict ordered list of at least two rays, rho+1 is the sum of C over
the lengths of its maximal constant-sign shortest-turn runs. No run length
ending at a listed ray is an integer, because any two listed rays are
linearly independent. Starting from a list of at least two rays, we prove
that inserting one ray increases rho by
between 0 and 4; the bound is 3 when inserted second or second-to-last,
and 2 when inserted at an endpoint.

At an endpoint the new turn has length <1. It either adds a singleton run
(cost 1) or extends the neighboring run (0 or 2); this gives [0,2].

For an interior insertion, reverse all turn signs if necessary so the old
shortest turn has length d in (0,1). Let the two new signed turns be a,b
in (-1,1). Their sum is congruent to d modulo 2.

## Same direction

If a,b>0, necessarily a+b=d. The original angular arc is merely subdivided
and rho is unchanged.

## A spike

Suppose a>0>b. Necessarily a=d+u, b=-u for 0<u<1-d. Let A>=0 be the
length of the positive run attached just before d, and let B>0 be the
length of the run immediately after the old turn, if there is one.

If that following run is negative, the cost change is

    2([frac(A+d)+u>=1] + [frac(B)+u>=1]),

which is 0, 2 or 4. If it is positive, one positive run splits into three
and the change is

    2+2(floor(A+d+u)+floor(B)-floor(A+d+B)),

again 0, 2 or 4: subtracting the integer part of B leaves the difference
of two floors with arguments differing by u-frac(B) in (-1,1). If there is
no following turn, the change is

    1+2(floor(A+d+u)-floor(A+d)),

which is 1 or 3. If the inserted ray is second, A=0 and d+u<1; the first
two cases are then at most 2, and the missing-following case at most 3.
The case a<0<b is the same calculation after reversing the ray list. Thus
the spike contributes at most 4 in general and at most 3 next to either
endpoint.

## The long opposite arc

If a,b<0, the two new turns form one negative run of total length 2-d.
Let A,B>0 be the lengths of the immediately neighboring old runs, when
present, and write x=frac(A), y=frac(B). According to their signs relative
to the old positive d-turn, direct subtraction of C costs gives

    ++ : 4-2 floor(x+y+d),
    +- : 2(1 + [y>=d] - [x+d>=1]),
    -+ : 2(1 + [x>=d] - [y+d>=1]),
    -- : 2 floor(x+y+1-d).

Each expression is 0, 2 or 4. If just one neighboring run is absent,
directly applying C instead gives 1 or 3. If both are absent, it gives 2.
This proves both the general and the near-endpoint bounds for this case.

These three possibilities exhaust the signs of a,b, so every insertion
has nonnegative cost at most 4. Removing the inserted ray gives the stated
deletion bound. The first-ray reflection bound follows from the endpoint
bound and parity, as explained in `independent_main_proof.md`. ∎

## Independent finite diagnostics

The companion `insertion_exhaustive_check.py` enumerates oriented ray
orders via two antipodal circular copies, without using coordinate-angle
sorting. Through seven rays it checks every deletion, both boundary bounds
and first-ray reflection. Finite enumeration is diagnostic; the arbitrary
length proof is the case split above.
