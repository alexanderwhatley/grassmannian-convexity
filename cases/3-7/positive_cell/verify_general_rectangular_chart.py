#!/usr/bin/env python3
"""Exact tests of conventions in general_positive_chart_and_global_rank.md.

Finite examples only; the proof for arbitrary dimensions is in the note.
"""
import itertools
import random
import sympy as s


def rect(n,k):
    r=n-k
    return [j for h in range(k) for j in range(r+h,h,-1)]


def parabolic(n,k):
    r=n-k
    return ([j for h in range(r-1) for j in range(1,r-h)] +
            [r+j for h in range(k-1) for j in range(1,k-h)])


def permutation(n,word):
    p=list(range(n))
    for j in word:
        p[j-1],p[j]=p[j],p[j-1]
    return p


def shear(n,j,t):
    u=s.eye(n)
    u[j,j-1]=t
    return u


def invcount(a):
    return sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))


nperm=0
for n in range(2,31):
    for k in range(1,n):
        r=n-k
        p=parabolic(n,k)
        w=rect(n,k)
        assert len(w)==k*r
        assert invcount(permutation(n,w))==len(w)
        assert permutation(n,w)==list(range(r,n))+list(range(r))
        assert len(p)+len(w)==n*(n-1)//2
        assert permutation(n,p+w)==list(reversed(range(n)))
        nperm+=1
print('PASS permutation, reduced-word and block conventions:',nperm,'pairs')

a,b,c=s.symbols('a b c',positive=True)
assert (shear(3,1,a)*shear(3,2,b)*shear(3,1,c)-shear(3,2,b*c/(a+c))*shear(3,1,a+c)*shear(3,2,a*b/(a+c))).applyfunc(s.simplify)==s.zeros(3)
print('PASS symbolic lower-shear braid transition')

ncharts=nminors=0
for n in range(2,8):
    for k in range(1,n):
        r=n-k
        w=rect(n,k)
        ts=[s.Integer((3*j+2*n+k)%7+1) for j in range(len(w))]
        pref=[s.eye(n)]
        for j,t in zip(w,ts):
            pref.append(pref[-1]*shear(n,j,t))
        suff=[None]*(len(w)+1)
        suff[-1]=s.eye(n)
        for q in range(len(w)-1,-1,-1):
            suff[q]=shear(n,w[q],ts[q])*suff[q+1]
        assert pref[-1]==suff[0]
        Y=pref[-1][r:n,:]
        Z=Y[:,r:n]
        A=Z.inv()*Y[:,:r]
        for I in itertools.combinations(range(n),k):
            assert Y[:,I].det()>0
            nminors+=1
        cols=[]
        for q,j in enumerate(w):
            e=s.zeros(n);e[j,j-1]=1
            dy=(pref[q]*e*suff[q+1])[r:n,:]
            dA=Z.inv()*(dy[:,:r]-dy[:,r:n]*A)
            cols.append(s.Matrix(list(dA)))
        jac=s.Matrix.hstack(*cols)
        assert jac.det()!=0
        if k==1:
            x=[s.prod(ts[:n-i-1]) for i in range(n-1)]+[s.Integer(1)]
            assert list(Y)==x
        ncharts+=1
print('PASS positive maximal minors:',nminors,'at',ncharts,'chart points')
print('PASS exact nonzero full chart Jacobian at all',ncharts,'points')
print('PASS exact k=1 chart formula at all',6,'sizes')

# Positivity of the coefficient used in the general backward construction.
nflag=0
for n in range(2,9):
    for k in range(1,n):
        d=k*(n-k)
        block=s.eye(n)
        for j in reversed(range(1,n)):
            block=block*shear(n,j,1)
        rev=block**d
        for I in itertools.combinations(range(n),k):
            assert rev[n-k:n,I].det()>0
            nflag+=1
print('PASS positive backward dominant coefficients:',nflag,'row minors')

# Numerical instances of Cauchy--Binet domination, with rational arithmetic.
rng=random.Random(37723)
nexamples=0
for n in range(2,8):
    for k in range(1,n):
        r=n-k
        index=list(itertools.combinations(range(n),k))
        while True:
            M=s.Matrix(k,n,lambda i,j:rng.randrange(-8,9))
            minors=[M[:,I].det() for I in index]
            if all(minors):
                break
        if minors[-1]<0:
            M[0,:]=-M[0,:]
        d=k*r
        L=s.eye(n)
        for sweep in range(d):
            for j in range(1,n):
                L=L*shear(n,j,1)
        eps=s.Rational(1,10**30)
        H=s.diag(*([1]*r+[eps]*k))
        M0=M*H.inv()*L.inv()*H
        sum_last=sum(range(r+1,n+1))
        for I in index:
            value=M0[:,I].det()
            sign=(-1)**(sum(i+1 for i in I)+sum_last)
            assert sign*value>0
        nexamples+=1
print('PASS backward epsilon dominance from generic matrices:',nexamples,'pairs')
print('All finite checks passed.')
