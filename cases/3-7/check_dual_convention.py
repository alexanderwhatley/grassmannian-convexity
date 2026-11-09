#!/usr/bin/env python3
"""Exact sample check of the signed-reversed Gale-dual convention; requires SymPy.
Diagnostic only: the exhaustive state check is verify_dual_bound.cpp.
"""
import random,itertools
import sympy as s
from math import comb

def var(seq):
 return sum(a!=b for a,b in zip(seq,seq[1:]))
def sign(x):
 return 1 if x>0 else -1 if x<0 else 0

def karp(N):
 k,n=N.shape
 vals=[]
 for js in itertools.combinations(range(n),k-1):
  seq=[]
  for i in range(n):
   if i in js: continue
   p=N[:,sorted((*js,i))].det()
   seq.append(sign(p)*(-1 if (i+1)%2 else 1))
  if 0 in seq: return None
  vals.append(var(seq))
 return n-k-max(vals)

def circuit_m(N):
 k,n=N.shape
 signs={p:sign(N[:,p].det()) for p in itertools.combinations(range(n),k)}
 if not all(signs.values()):return None
 for v in range(n):
  for mask in range(1<<(n-1)):
   x=[-1 if mask>>i&1 else 1 for i in range(n-1)]+[1]
   if var(x)!=v:continue
   if all(len({x[i]*(-1 if j%2 else 1)*signs[c[:j]+c[j+1:]] for j,i in enumerate(c)})>1 for c in itertools.combinations(range(n),k+1)):
    return v
 raise RuntimeError('No covectors')
# The change of coordinates converts every elementary move exactly.
J=s.zeros(7)
for i in range(7): J[i,6-i]=1
D=s.diag(*[(-1)**i for i in range(1,8)])
C=J*D
for j in range(1,7):
 U=s.eye(7);U[j,j-1]=s.Rational(5,3)
 Udual=s.eye(7);Udual[7-j,6-j]=s.Rational(5,3)
 assert C.inv()*U.inv().T*C==Udual
print('PASS signed-reversed dual conjugates each of the six positive moves')
random.seed(764)
count=0
while count<30:
 M=s.Matrix([[random.randrange(-7,8) for i in range(7)] for j in range(3)])
 if M.rank()!=3 or any(M[:,I].det()==0 for I in itertools.combinations(range(7),3)):continue
 B=s.Matrix.hstack(*M.nullspace()).T
 N=B[:,::-1]*s.diag(*[(-1)**i for i in range(1,8)])
 ratios=set()
 for I in itertools.combinations(range(1,8),4):
  Cset={8-i for i in I}
  complement=[i-1 for i in range(1,8) if i not in Cset]
  ratios.add(sign(N[:,[i-1 for i in I]].det())*sign(M[:,complement].det()))
 assert len(ratios)==1
 m1=karp(N);m2=circuit_m(N)
 # The direct complementary-minor formula used by the full-state check.
 vals=[]
 for j in itertools.combinations(range(1,8),3):
  v=[]
  for i in range(1,8):
   if i in j:continue
   subset={8-x for x in (*j,i)}
   I=[x-1 for x in range(1,8) if x not in subset]
   v.append(sign(M[:,I].det())*(-1 if i%2 else 1))
  vals.append(var(v))
 assert m1==m2==3-max(vals),(m1,m2,vals)
 count+=1
print('PASS 30 exact generic matrices: dual minor signs, nullspace circuits, and Karp/complementary-minor formulas agree')
M=s.Matrix([[(-1)**i*i**j for i in range(1,8)] for j in range(3)])
b=s.Matrix([comb(6,i) for i in range(7)])
assert M*b==s.zeros(3,1) and all(x>0 for x in b)
print('PASS ordinary kernel has a strictly positive vector for the 28-move witness; raw dual variation is 0')
