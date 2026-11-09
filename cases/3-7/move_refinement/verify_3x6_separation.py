#!/usr/bin/env python3
"""Exact nonzero event-root resultants for 20 minors and 3 F features.

All 23 polynomials are multilinear in the six columns, hence affine in each
elementary shear parameter. An exact nonzero value of every pairwise root
resultant proves that it is a nonzero polynomial and generic events separate.
"""
from itertools import combinations
if not __debug__:raise RuntimeError('assertions are required')
M=[[38,7,-22],[27,-26,4],[23,-28,13],[33,-15,-4],[19,12,-24],[-3,34,31]]
triples=list(combinations(range(1,7),3))
def det(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def vals(cols):
    p={I:det(*(cols[i-1] for i in I)) for I in triples}
    return [p[I] for I in triples]+[
        p[1,2,5]*p[3,4,6]-p[1,2,6]*p[3,4,5],
        p[1,3,5]*p[2,4,6]-p[1,3,6]*p[2,4,5],
        p[1,4,5]*p[2,3,6]-p[1,4,6]*p[2,3,5]]
value=vals(M);assert all(value);counts=[];checked=0
for j in range(1,6):
    one=[x[:] for x in M];one[j-1]=[a+b for a,b in zip(one[j-1],one[j])]
    two=[x[:] for x in M];two[j-1]=[a+2*b for a,b in zip(two[j-1],two[j])]
    a=vals(one);b=vals(two);delta=[v-u for u,v in zip(value,a)]
    assert all(z-2*y+x==0 for x,y,z in zip(value,a,b))
    active=[i for i,d in enumerate(delta) if d];counts.append(len(active))
    for i,k in combinations(active,2):
        assert value[i]*delta[k]-value[k]*delta[i]!=0
        checked+=1
assert counts==[8,8,8,9,6] and checked==135
print('PASS all 23 feature values nonzero; active events',counts)
print('PASS',checked,'exact nonzero same-generator event-root resultants')
