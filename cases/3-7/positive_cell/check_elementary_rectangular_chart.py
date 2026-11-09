#!/usr/bin/env python3
"""Exact finite checks for the elementary rectangular-chart argument.

Paths are enumerated by dynamic programming over k-subsets, not by the
symbolic product/flag-factorization programs used elsewhere. A count of one
also records the unique set of shear positions selected by that path.
"""
from itertools import combinations

cases=selected=0
for n in range(2,11):
    for k in range(1,n):
        r=n-k
        word=[j for i in range(k) for j in range(r+i,i,-1)]
        start=tuple(range(r+1,n+1))
        counts={start:(1,0)}
        for step,j in enumerate(word):
            out={}
            def add(I,c,mask):
                if I not in out:out[I]=(c,mask if c==1 else None)
                else:out[I]=(out[I][0]+c,None)
            for I,(c,mask) in counts.items():
                add(I,c,mask)
                if j+1 in I and j not in I:
                    target=tuple(sorted((set(I)-{j+1})|{j}))
                    add(target,c,(mask|(1<<step)) if mask is not None else None)
            counts=out
        assert len(counts)==len(list(combinations(range(1,n+1),k)))
        assert all(c>0 for c,mask in counts.values())
        for h in range(k+1):
            for ell in range(r+1):
                I=tuple(range(1,h+1))+tuple(range(h+ell+1,k+ell+1))
                expected=0
                for i in range(k):
                    for j in range(r):
                        if i<h or j<r-ell:expected|=1<<(i*r+j)
                assert counts[I]==(1,expected),(n,k,h,ell,I,counts[I],expected)
                selected+=1
        cases+=1
print('PASS',cases,'rectangular words: every maximal minor has a positive path')
print('PASS',selected,'initial solid minors: exactly one path with the predicted monomial')
