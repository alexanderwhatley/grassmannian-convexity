#!/usr/bin/env python3
"""Construct exact generic reduced-word witnesses attaining all s3 flips.

The path construction is from Gaetz--Gao--Jiradilok--Nenashev--Postnikov,
arXiv:2204.03033v2, sec. 3.5.  Here the new ingredients are an explicit
rank-three sign realization and a uniform-time totally-nonnegative flag
minor limit, explained in reduced_sharp_all_n.md.  This program independently
replays finite cases with integer arithmetic.
"""
from itertools import combinations
import sys


BASE = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, -1, -1),
        (-1, 0, 1), (-2, -1, -1)]
REPEAT = ('123', '124', '125', '145', '245', '345', '456',
          '457', '567', '578', '678', '789')
TAIL = {
    4: ('123', '124', '134', '234'),
    5: ('123', '124', '125', '145', '245', '345'),
    6: ('123', '124', '125', '145', '245', '345', '456'),
    7: ('123', '124', '125', '145', '245', '345', '456', '457', '567'),
    8: ('123', '124', '125', '145', '245', '345', '456', '457',
        '567', '578', '678'),
    9: ('123', '124', '125', '145', '245', '345', '456', '457',
        '567', '578', '579', '589', '789'),
}


def triples(strings, offset):
    return [tuple(int(c) + offset for c in s) for s in strings]


def path(n):
    m = (n - 4) // 6
    r = n - 6 * m
    out = [(1, 2, 3)]
    for h in range(m):
        out += triples(REPEAT, 6*h)[1:]
    out += triples(TAIL[r], 6*m)[1:]
    assert len(out) - 1 == (11*n + 5)//6 - 5
    return out


def det(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - b[0]*(a[1]*c[2]-a[2]*c[1])
            + c[0]*(a[1]*b[2]-a[2]*b[1]))


def raw_source(n):
    columns = []
    for i in range(n):
        h, r = divmod(i, 6)
        x, y, z = BASE[r]
        columns.append((x, y * (-1)**h, z))
    if n % 6 == 3:
        h = (n - 9)//6
        assert h >= 0
        columns[-1] = (0, (-1)**h, -1)
    return columns


def det_at(columns, I):
    return det(*(columns[i-1] for i in I))


def source(n, P):
    base = raw_source(n)
    assert all(det_at(base, I) * (-1)**j > 0 for j, I in enumerate(P))
    # A single perturbation by a moment curve avoids every cubic determinant
    # polynomial's roots.  Multiply through by the common denominator R.
    R = 100*n*n
    all_I = tuple(combinations(range(1, n + 1), 3))
    while True:
        V = [(R*x+1, R*y+i, R*z+i*i)
             for i, (x, y, z) in enumerate(base, 1)]
        if (all(det_at(V, I) for I in all_I)
                and all(det_at(V, I) * (-1)**j > 0
                        for j, I in enumerate(P))):
            return V, R, all_I
        R += 1


def lift_word(n, P):
    # The constructive path lift of Proposition 2.3 in the cited paper.
    perm = list(range(1, n + 1))
    word = []
    def add(j):
        assert perm[j-1] < perm[j], (j, perm, P)
        perm[j-1], perm[j] = perm[j], perm[j-1]
        word.append(j)
    for first, second in zip(P, P[1:]):
        gone = (set(first) - set(second)).pop()
        new = (set(second) - set(first)).pop()
        a, b = perm.index(gone) + 1, perm.index(new) + 1
        assert a <= 3 < b and gone < new
        for j in range(a, 3):
            add(j)
        for j in range(b - 1, 3, -1):
            add(j)
        add(3)
        assert set(perm[:3]) == set(second)
    assert word.count(3) == len(P)-1
    # The final 3-set contains the largest labels, so the only remaining
    # ascending adjacent pairs lie within the prefix or within the suffix.
    # Sorting these blocks descending completes the word to w0 without s3.
    while True:
        ascents = [j for j in range(1, n)
                   if j != 3 and perm[j-1] < perm[j]]
        if not ascents:
            break
        add(ascents[0])
    assert perm == list(range(n, 0, -1))
    assert len(word) == n*(n-1)//2 and word.count(3) == len(P)-1
    return word


def replay(V, word, t, all_I):
    columns = V.copy()
    previous = det(*columns[:3])
    flips = 0
    for j in word:
        columns[j-1] = tuple(a + t*b
                             for a, b in zip(columns[j-1], columns[j]))
        current = det(*columns[:3])
        if j == 3 and previous * current >= 0:
            return None
        if j != 3 and previous != current:
            raise AssertionError('target changed at a non-U3 generator')
        flips += previous * current < 0
        previous = current
        if not all(det_at(columns, I) for I in all_I):
            return None
    assert flips == word.count(3)
    return flips


def check(n):
    P = path(n)
    V, R, all_I = source(n, P)
    word = lift_word(n, P)
    t = 2
    while True:
        flips = replay(V, word, t, all_I)
        if flips is not None:
            break
        t *= 2
    print(f'n={n:2d} q={n-3:2d} flips={flips:3d} '
          f'shears={len(word):4d} R={R} common_t={t}', flush=True)


if __name__ == '__main__':
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    for n in range(4, limit + 1):
        check(n)
