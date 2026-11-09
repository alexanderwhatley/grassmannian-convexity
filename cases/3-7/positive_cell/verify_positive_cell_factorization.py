#!/usr/bin/env python3
"""Symbolically verify a 12-parameter factorization of Gr(3,7)>0.

For X=[A|I3], every ordered 3x3 minor positive, the displayed monomial
minor ratios are positive and reconstruct X as E times the indicated
positive lower elementary product, up to the usual row normalization.
"""
import sympy as s
from itertools import combinations
x=s.symbols('x0:12');A=s.Matrix(3,4,x);X=A.row_join(s.eye(3))
p=lambda a,b,c:s.expand(X[:,[a-1,b-1,c-1]].det())
T=[p(1,6,7)*p(4,5,6)/p(1,5,6),
   p(1,5,6)*p(3,4,5)/(p(1,4,5)*p(4,5,6)),
   p(1,4,5)*p(2,3,4)/(p(1,3,4)*p(3,4,5)),
   p(1,3,4)/p(2,3,4),
   p(1,2,7)*p(1,5,6)/(p(1,2,6)*p(1,6,7)),
   p(1,2,6)*p(1,4,5)/(p(1,2,5)*p(1,5,6)),
   p(1,2,5)*p(1,3,4)/(p(1,2,4)*p(1,4,5)),
   p(1,2,4)/p(1,3,4), p(1,2,6)/p(1,2,7),
   p(1,2,5)/p(1,2,6),p(1,2,4)/p(1,2,5),p(1,2,3)/p(1,2,4)]
word=[4,3,2,1,5,4,3,2,6,5,4,3]
# Use a symbolic product with independent variables first, then substitute.
a=s.symbols('a0:12');U=s.eye(7)
for j,t in zip(word,a):
 for i in range(7):U[i,j-1]+=t*U[i,j]
Y=U[4:7,:];Q=s.simplify(Y[:,4:7].inv()*Y[:,:4]);sub=dict(zip(a,T))
for i in range(3):
 for j in range(4):
  diff=s.factor(s.cancel(Q[i,j].subs(sub,simultaneous=True)-A[i,j]))
  assert diff==0,(i,j,diff)
  print('PASS coordinate',i,j,flush=True)
print('PASS all 12 rational inverse identities for the positive-cell elementary factorization')
