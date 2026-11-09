# Data-free certificate for Grassmann convexity (3,7)

**Status:** candidate computer-assisted proof, not independently peer
reviewed. `verify_upper_bound.cpp` reconstructs its entire finite sign graph
and its candidate rank from scratch. It reads no certificate, event table,
or other input data file. It directly checks every rank inequality and
the sharper bound $r(M)\le8+m(M)$, where $m$ is the minimum sign
variation of a nonzero-evaluation row-space vector.

Run:

```sh
./run.sh
```

Requirements: a GCC/Clang-compatible C++17 compiler with pthreads,
Python 3.10 or newer, and sufficient RAM.
On the reference machine, the graph check took about 35 seconds and peaked
at approximately 2.4 GiB; other machines may take longer. Compile without
`-DNDEBUG`: the source deliberately requires assertions. The included
historical sanitizer log records output attributed to an address/undefined-
behavior run; the log alone does not establish its provenance or flags.

For a machine with less memory, `./run_low_memory.sh` checks the same upper
bound using two bundled compressed rank arrays. This version reconstructs
the states and verifies every inequality directly, so it does not require
the arrays' generation to be trusted. It checks the upper bound only; run
the Python witness checks separately for the lower bound. On the reference
machine the original version of this check used about 565 MiB. It requires
`gzip` in addition to the compiler.

The checker derives all 114 nonzero move derivatives with exact symbolic
polynomials, verifies generic event separation, constructs 23,246,848
states and 144,788,896 distinct directed sign-flip edges, and checks
230,629,888 move-labeled transitions. It computes a nonnegative integer
rank bounded by $8+m$, dropping at each target sign change. The theorem
still relies on the mathematical coverage and perturbation arguments in
the accompanying proof report, as well as the correctness of the program.

`verify_sharpness.py` independently replays `sharp_trajectory.json` using
exact rational/integer arithmetic. The trajectory has 207 positive moves,
12 target changes, and all 52 tracked features nonzero at every vertex.
The suffixes beginning after moves 0, 15, 36, 57, and 89 start with
$m=4,3,2,1,0$ and respectively have 12, 11, 10, 9, and 8 target changes.
Thus, subject to independent review of the upper certificate, the bound is
sharp in each class and gives $r(3,7)=12$.

`compact_witness.json` is a shorter 105-move rational witness. It has 12
target changes, all 52 tracked features strict at each vertex, and its
suffixes are sharp in all five variation classes. `verify_compact_witness.py`
replays it by evolving integer column representatives and taking exact
determinants. The still shorter `short_witness.json` has an integer starting
matrix, 29 integer move sizes, and 12
changes. The optimal `optimal_witness.json` has 28 positive rational moves,
and `optimal_witness_integer.json` is an all-integer version. These and the
12-move preparation give a 40-move path from every generic start with $m=4$.

`simple_moment_12.json` gives a second exact 12-crossing sequence, starting
from the particularly simple columns $(-1)^i(1,i,i^2)$, for
$i=1,\ldots,7$. `verify_moment.py` independently replays its 219 positive
rational moves and checks all maximal minors at every sampled vertex.
`simple_moment_28.json` is an optimal 28-move continuation from the same columns;
`simple_moment_40.json` is an example produced by the general constructor.

The proof report also shows that every generic start with $m=4$ has a
12-crossing continuation, by a separate positive-cell transfer argument.
The `positive_cell/` directory contains two independent symbolic checks of
the 12-parameter chart and a Jacobian check. These optional checks require
SymPy, while `run.sh` and the data-free upper certificate do not:

```sh
python3 positive_cell/verify_positive_cell_factorization.py
python3 positive_cell/audit_positive_cell_elimination.py
python3 positive_cell/verify_positive_cell_jacobian.py
python3 positive_cell/verify_elementary_positive_chart.py
```

For any rational generic matrix in the maximum-variation chamber,
`positive_cell/exact_ift_construct_rational.py` constructs a 40-move
12-crossing trajectory using exact arithmetic. It takes a JSON file of seven
three-component columns, the included 28-move witness, and an output path:

```sh
python3 positive_cell/exact_ift_construct_rational.py \
  positive_cell/matrix_simple.json optimal_witness.json result.json
python3 positive_cell/verify_constructed_12.py result.json
```

The script checks every output sign exactly. Its termination argument and
input format are described in `positive_cell/exact_rational_constructor.md`.
The constructor and its verifier use Python's standard library only.

The [dual-bound note](Dual_Bound_and_Conjecture.md) describes an additional
exhaustive table, its positive-move dual convention, and how to reproduce it.
The dual checker adapts this package’s graph construction; it is not a
third-party reimplementation.
