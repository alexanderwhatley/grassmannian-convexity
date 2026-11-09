#!/usr/bin/env python3
"""Independent finite checks of the subset-order, Pluecker, and sign-ideal lemma.

Random integer Vandermonde-reoriented starts, exact minors, and exact rational
samples from every open interval cut by a shear's determinant roots. These
checks illustrate the general proof; they do not replace it.
"""
from itertools import combinations, product
from fractions import Fraction
from math import prod
import random

random.seed(491793)

def det(a):
    n=len(a)
    if n==0:return 1
    if n==1:return a[0][0]
    return sum((-1)**j*a[0][j]*det([row[:j]+row[j+1:] for row in a[1:]]) for j in range(n))

def leq(a,b):return all(x<=y for x,y in zip(a,b))
def advance(a,j):
    if j not in a or j+1 in a:return None
    return tuple(sorted((set(a)-{j})|{j+1}))
def covers(a,subsets):
    return [b for j in a if (b:=advance(a,j)) is not None and b in subsets]
def check_ideal(signs,edges):
    return all(not (signs[a]>0>signs[b]) for a,b in edges)

structural=identities=ideal_starts=intervals=simultaneous=0
for n in range(3,10):
 for k in range(1,min(n,5)):
    S=list(combinations(range(1,n+1),k)); ss=set(S)
    edges=[(a,b) for a in S for b in covers(a,ss)]
    X=[[random.randint(-8,8) for _ in range(n)] for _ in range(k)]
    p={a:det([[X[i][j-1] for j in a] for i in range(k)]) for a in S}
    for a,b in edges:
     for j in range(1,n):
        aa,bb=advance(a,j),advance(b,j)
        structural+=1
        if aa and not bb:assert aa==b
        if aa and bb:
            assert leq(aa,bb)
            x=next(iter(set(a)-set(b)))
            assert abs(x-j)>=2
            common=set(a)&set(b)-{j}
            low=tuple(sorted(common|{min(x,j),min(x,j)+1}))
            high=tuple(sorted(common|{max(x,j),max(x,j)+1}))
            assert len(low)==len(high)==k and leq(low,a) and leq(aa,high)
            assert p[a]*p[bb]-p[b]*p[aa]==p[low]*p[high]
            identities+=1
    if n>7 or k>4:continue
    # Alternate many exactly computable sign ideals. For random moment curves,
    # every maximal minor is Vandermonde-positive prior to column sign flips.
    for rep in range(5):
     x=sorted(random.sample(range(-30,31),n))
     for colors in product((-1,1),repeat=n):
        # Multiplying an arbitrary row by -1 normalizes p_last>0.
        z={a:prod(colors[v-1] for v in a)*prod(x[v-1]-x[u-1] for u,v in combinations(a,2)) for a in S}
        if z[S[-1]]<0:z={a:-v for a,v in z.items()}
        if not check_ideal(z,edges):continue
        ideal_starts+=1
        for j in range(1,n):
            active={a:advance(a,j) for a in S}
            roots=sorted({Fraction(-z[a],z[b]) for a,b in active.items() if b is not None and -z[a]*z[b]>0})
            bounds=[Fraction(0)]+roots
            samples=[(x+y)/2 for x,y in zip(bounds,bounds[1:])]+[bounds[-1]+1]
            for t in samples:
                zz={a:z[a]+t*z[b] if (b:=active[a]) is not None else z[a] for a in S}
                assert all(zz.values()) and check_ideal(zz,edges)
                assert all(zz[a]>0 for a in S if z[a]>0)
                intervals+=1
            simultaneous+=sum(1 for a,b in edges if active[a] is not None and active[b] is not None and z[a]<0 and z[b]<0 and z[active[a]]>0 and z[active[b]]>0)
print('structural cover/shear cases:',structural)
print('exact 2-by-2 Pluecker identity cases:',identities)
print('generic sign-ideal starting points:',ideal_starts)
print('strict intervals checked:',intervals)
print('hard-case simultaneous-active covers:',simultaneous)
print('PASS')
