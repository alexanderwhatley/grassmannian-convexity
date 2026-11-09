#!/usr/bin/env python3
"""Check the general positive-chart inverse on exact Vandermonde inputs."""
from fractions import Fraction as Q
from itertools import combinations

def det(A):
    A=[list(map(Q,row)) for row in A];n=len(A);sign=1;value=Q(1)
    for j in range(n):
        p=next((i for i in range(j,n) if A[i][j]),None)
        if p is None:return Q(0)
        if p!=j:A[p],A[j]=A[j],A[p];sign=-sign
        a=A[j][j];value*=a
        for i in range(j+1,n):
            c=A[i][j]/a
            for s in range(j+1,n):A[i][s]-=c*A[j][s]
    return sign*value

def minors(A,n,k):
    return {I:det([[A[a][j-1] for j in I] for a in range(k)])
            for I in combinations(range(1,n+1),k)}

cases=checks=0
for n in range(2,10):
    for k in range(1,n):
        r=n-k
        X=[[Q(j)**i for j in range(1,n+1)] for i in range(k)]
        p=minors(X,n,k);normal=p[tuple(range(r+1,n+1))]
        p={I:q/normal for I,q in p.items()}
        def C(h,l):
            return p[tuple(range(1,h+1))+tuple(range(h+l+1,k+l+1))]
        parameters=[]
        for i in range(1,k+1):
            for j in range(1,r+1):
                l=r-j+1;h=i-1
                t=C(h,l-1)*C(h+1,l)/(C(h,l)*C(h+1,l-1))
                assert t>0;parameters.append(t)
        Y=[[Q(int(j==r+i)) for j in range(1,n+1)] for i in range(1,k+1)]
        word=[j for b in range(k) for j in range(r+b,b,-1)]
        for j,t in zip(word,parameters):
            for i in range(k):Y[i][j-1]+=t*Y[i][j]
        q=minors(Y,n,k)
        assert q==p,(n,k)
        checks+=len(q);cases+=1
print('PASS exact rational inverse and all',checks,'Pluecker coordinates in',cases,'Vandermonde charts')
