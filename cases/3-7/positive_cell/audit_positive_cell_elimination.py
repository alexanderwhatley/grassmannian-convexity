#!/usr/bin/env python3
"""Independent direct elimination audit of the positive-cell 12-word chart.

Given symbolic X=[A|I], compute the proposed positive Plucker-ratio
parameters, right-multiply by the inverse elementary word, and verify that
the first four columns vanish. This differs from reconstructing X by a
symbolic forward product. SymPy is used for exact rational cancellation.
"""
import sympy as s
x=s.symbols('z0:12');X=s.Matrix(3,4,x).row_join(s.eye(3))
def p(a,b,c):return s.expand(X[:,[a-1,b-1,c-1]].det())
t=[p(1,6,7)*p(4,5,6)/p(1,5,6),
 p(1,5,6)*p(3,4,5)/(p(1,4,5)*p(4,5,6)),
 p(1,4,5)*p(2,3,4)/(p(1,3,4)*p(3,4,5)),
 p(1,3,4)/p(2,3,4),
 p(1,2,7)*p(1,5,6)/(p(1,2,6)*p(1,6,7)),
 p(1,2,6)*p(1,4,5)/(p(1,2,5)*p(1,5,6)),
 p(1,2,5)*p(1,3,4)/(p(1,2,4)*p(1,4,5)),
 p(1,2,4)/p(1,3,4),p(1,2,6)/p(1,2,7),
 p(1,2,5)/p(1,2,6),p(1,2,4)/p(1,2,5),p(1,2,3)/p(1,2,4)]
word=[4,3,2,1,5,4,3,2,6,5,4,3]
Y=X.copy()
for j,a in reversed(list(zip(word,t))):
 for i in range(3):Y[i,j-1]=s.factor(s.cancel(Y[i,j-1]-a*Y[i,j]))
for i in range(3):
 for j in range(4):assert Y[i,j]==0,(i,j,Y[i,j])
assert s.factor(Y[:,4:7].det())==1
print('PASS symbolic inverse word eliminates all first-four columns; remaining determinant 1')
# Independently check inverse ratios on a symbolic forward product.
a=s.symbols('a0:12');T=s.eye(7)
for j,z in zip(word,a):
 for i in range(7):T[i,j-1]+=z*T[i,j]
P=T[4:7,:]
def q(i,j,k):return s.factor(P[:,[i-1,j-1,k-1]].det())
r=[q(1,6,7)*q(4,5,6)/q(1,5,6),
 q(1,5,6)*q(3,4,5)/(q(1,4,5)*q(4,5,6)),
 q(1,4,5)*q(2,3,4)/(q(1,3,4)*q(3,4,5)),
 q(1,3,4)/q(2,3,4),
 q(1,2,7)*q(1,5,6)/(q(1,2,6)*q(1,6,7)),
 q(1,2,6)*q(1,4,5)/(q(1,2,5)*q(1,5,6)),
 q(1,2,5)*q(1,3,4)/(q(1,2,4)*q(1,4,5)),
 q(1,2,4)/q(1,3,4),q(1,2,6)/q(1,2,7),
 q(1,2,5)/q(1,2,6),q(1,2,4)/q(1,2,5),q(1,2,3)/q(1,2,4)]
for i in range(12):assert s.cancel(r[i]-a[i])==0,i
print('PASS all 12 forward Plucker-ratio identities; chart map has explicit positive rational inverse')

from itertools import combinations
for I in combinations(range(1,8),3):
 polynomial=s.Poly(s.expand(q(*I)),*a)
 assert polynomial and all(c>0 for c in polynomial.coeffs()),I
print('PASS every one of the 35 maximal minors is a positive-coefficient polynomial')
