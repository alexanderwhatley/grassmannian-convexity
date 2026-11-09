# Grassmann convexity conjecture

This repository collects research on the Grassmann convexity conjecture. Its first case study is a **candidate computer-assisted proof for (k,n) = (3,7)**. The argument and code have not been independently peer reviewed or formally verified. In particular, the mathematical coverage of the finite graph, the perturbation argument, and the checker itself need external scrutiny.

## The problem

Consider an nth-order real linear ordinary differential equation that is *disconjugate* on an interval: every nonzero solution has at most n−1 zeros there, counting multiplicities. The Grassmann convexity conjecture, formulated by Boris and Michael Shapiro, asserts that the Wronskian of any k linearly independent solutions has at most **k(n−k)** zeros on that interval, counted with multiplicity. The number k(n−k) is also the dimension of the Grassmannian Gr(k,n).

An equivalent discrete formulation tracks the signs of maximal minors of a real k × n matrix under positive elementary column moves. In the (3,7) case, write its columns as v₁,…,v₇ and pᵢⱼₖ = det(vᵢ,vⱼ,vₖ). A move replaces vⱼ by vⱼ + t vⱼ₊₁, with t > 0. The question is whether p₁₂₃ can change sign more than 12 times along any finite sequence of such moves. See the finiteness paper below for the precise formulations and equivalence.

## Prior work

- Saldanha, Boris Shapiro, and Michael Shapiro proved the conjecture for k = 2 and, by duality, k = n−2. See [*Grassmann convexity and multiplicative Sturm theory, revisited*](https://arxiv.org/abs/1902.09741).
- Their [*Finiteness of rank for Grassmann convexity*](https://arxiv.org/abs/2110.07389) gives an equivalent discrete formulation, the general lower bound k(n−k), and a general explicit upper bound. It also reports a computer-assisted verification of (3,6). The general conjecture remains the broader problem addressed by this repository.
- Steven Karp's [*Sign variation, the Grassmannian, and total positivity*](https://arxiv.org/abs/1503.05622) supplies sign-variation results used in the (3,7) analysis.

These references establish context; they should not be read as independent verification of the new certificate in this repository. No claim of priority or exhaustive literature review is made here.

## The (3,7) candidate

For a generic 3 × 7 matrix M, let m(M) be the minimum number of sign changes in the sequence (ω(v₁),…,ω(v₇)) over linear functionals ω that are nonzero on every column. The proposed certificate establishes the stronger upper bound **r(M) ≤ 8 + m(M)**, where r(M) is the maximum number of p₁₂₃ sign changes reachable from M. Since m(M) ≤ 4, this would give the upper bound 12. Exact rational examples give 12 changes, matching the known lower bound.

The proof uses a finite over-approximation of real trajectories. A state records 35 maximal-minor signs and 17 auxiliary signs. The checker reconstructs 23,246,848 states and directly checks 230,629,888 move-labeled transition inequalities using exact arithmetic. The reduction from real trajectories to this graph is a mathematical part of the argument, not a consequence of the computation alone. The detailed report explains the assumptions, checks, witnesses, and remaining review obligations.

## Files and reproduction

- [Proof and sharp bounds](cases/3-7/Proof_and_Sharp_Bounds.md): the candidate argument, definitions, witnesses, and further claims.
- [Review guide](cases/3-7/Review_Guide.md): a map of the mathematical and computational obligations.
- [Certificate instructions](cases/3-7/README.md): requirements, detailed checks, and lower-memory alternative.
- [Source-only upper-bound checker](cases/3-7/verify_upper_bound.cpp), [witnesses](cases/3-7/optimal_witness_integer.json), and supporting scripts are included in `cases/3-7/`.
- [Complete audit archive](cases/3-7/complete-audit.zip): additional source, data, examples, and historical verification logs.

To run the main package, from the repository root:

```sh
cd cases/3-7
./run.sh
```

The script needs Python 3.10+, a GCC/Clang-compatible C++17 compiler with pthreads, and substantial memory (about 2.4 GiB on the reference run). It compiles with assertions enabled. `./run_low_memory.sh` checks the upper bound with bundled rank data; see the certificate instructions for details. The archived logs describe earlier runs and are not a substitute for independent review.
