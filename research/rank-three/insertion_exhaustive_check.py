#!/usr/bin/env python3
"""Exact signed-determinant census for the planar one-ray lemma.

Enumerates every oriented projective order through seven rays, after fixing
one ray by a rotation. Rank is computed from determinant-sign changes along
monotone turn runs; no sorting of ray angles is used.
"""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import importlib.util
import random


def det(a,b): return a[0]*b[1]-a[1]*b[0]
def sign(a): return (a>0)-(a<0)


def rank(W):
    assert len(W)>=2 and all(det(a,b) for a,b in combinations(W,2))
    turns=[sign(det(a,b)) for a,b in zip(W,W[1:])]
    changes=sum(a!=b for a,b in zip(turns,turns[1:]))
    starts=[0]+[i for i in range(1,len(turns)) if turns[i]!=turns[i-1]]
    stops=starts[1:]+[len(turns)]
    carry=0
    for a,b in zip(starts,stops):
        signs=[sign(det(W[a],W[j])) for j in range(a+1,b+1)]
        carry+=sum(x!=y for x,y in zip(signs,signs[1:]))
    return changes+2*carry


if __name__=='__main__':
    spec=importlib.util.spec_from_file_location('earlier',Path(__file__).resolve().parent/'planar_rank_crosscheck.py')
    earlier=importlib.util.module_from_spec(spec);spec.loader.exec_module(earlier)
    rng=random.Random(41192);overall=0;deletions=0;crosschecks=0
    for n in range(3,8):
        counts=Counter();maxima=Counter();cases=0
        for p in permutations(range(1,n)):
            for signs in product((-1,1),repeat=n-1):
                W=[(0,1)]+[(s,s*j) for s,j in zip(signs,p)]
                r=rank(W)
                if n<=5 or rng.randrange(100)==0:
                    assert r==earlier.rho(W)
                    crosschecks+=1
                for i in range(n):
                    drop=r-rank(W[:i]+W[i+1:]);counts[drop]+=1
                    tag='endpoint' if i in (0,n-1) else 'near' if i in (1,n-2) else 'interior'
                    maxima[tag]=max(maxima[tag],drop)
                    assert 0<=drop<=4
                    if tag=='endpoint':assert drop<=2
                    if tag=='near':assert drop<=3
                reflected=[(-W[0][0],-W[0][1])]+W[1:]
                assert abs(rank(reflected)-r)<=1
                cases+=1
        overall+=cases;deletions+=n*cases
        print(f'n={n} orders={cases} deletion_drops={dict(sorted(counts.items()))} '
              f'category_max={dict(maxima)} PASS',flush=True)
    print(f'PASS total orders={overall} deletions={deletions} independent_rank_crosschecks={crosschecks}')
