#!/usr/bin/env python3
"""Circular-order planar rank cross-check extracted from the compression verifier."""
from fractions import Fraction as Q
from functools import cmp_to_key
from itertools import permutations, combinations


def det(u,v): return u[0]*v[1]-u[1]*v[0]
def sgn(a): return (a>0)-(a<0)

def rho(V):
    assert all(v!=(0,0) for v in V)
    assert all(det(u,v) for u,v in combinations(V,2))
    n=len(V)
    rays=[(x,y,i,1) for i,(x,y) in enumerate(V)] + [(-x,-y,i,-1) for i,(x,y) in enumerate(V)]
    def half(z): return 0 if z[1]>0 or z[1]==0 and z[0]>0 else 1
    def cmp(a,b):
        h,k=half(a),half(b)
        return h-k if h!=k else -sgn(det(a,b))
    rays.sort(key=cmp_to_key(cmp))
    pos=[0]*n
    for i,(_,_,j,e) in enumerate(rays):
        if e==1: pos[j]=i
    runs=[]
    for i in range(n-1):
        d=(pos[i+1]-pos[i])%(2*n)
        assert d not in (0,n)
        s,l=(1,d) if d<n else (-1,2*n-d)
        if runs and runs[-1][0]==s: runs[-1][1]+=l
        else: runs.append([s,l])
    return len(runs)-1+2*sum(l//n for _,l in runs)

