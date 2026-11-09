#!/usr/bin/env python3
"""Exact integer-column audit of one trajectory sharp in all five classes."""
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path
import json

data=json.loads(Path(__file__).with_name('sharp_trajectory.json').read_text())
v=[list(map(int,a)) for a in data['initial_columns']]
scale=[Fraction(1)]*7
triples=list(combinations(range(7),3))
quads=list(combinations(range(7),4))
features=[(1,2,3,4),(1,3,2,4),(1,4,2,3),(1,2,3,5),(1,3,2,5),
 (1,5,2,3),(1,2,4,5),(1,4,2,5),(1,5,2,4),(1,3,4,5),
 (1,4,3,5),(1,5,3,4),(2,3,4,5),(2,4,3,5),(2,5,3,4)]
def det(a,b,c):
 return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def evaluate():
 minors={I:det(*(v[j] for j in I)) for I in triples}
 assert all(minors.values())
 def p(a,b,c):return minors[a-1,b-1,c-1]
 for a,b,c,d in features:
  assert p(a,b,6)*p(c,d,7)-p(a,b,7)*p(c,d,6)
 assert p(1,2,5)*p(3,4,6)-p(1,2,6)*p(3,4,5)
 assert p(1,2,5)*p(3,4,7)-p(1,2,7)*p(3,4,5)
 sign={I:1 if a>0 else -1 for I,a in minors.items()}
 feasible=[]
 for bits in range(64):
  for a,b,c,d in quads:
   co=[sign[b,c,d],-sign[a,c,d],sign[a,b,d],-sign[a,b,c]]
   s=[co[k]*(-1 if (bits>>i)&1 else 1) for k,i in enumerate((a,b,c,d))]
   if s[0]==s[1]==s[2]==s[3]:break
  else:feasible.append((bits^(bits>>1)).bit_count())
 assert feasible
 return sign[0,1,2],min(feasible)

oldsign,oldm=evaluate();assert oldm==4
flips=[];drops=[]
for step,(j,t) in enumerate(data['moves'],1):
 j=int(j)-1;t=Fraction(t);assert 0<=j<6 and t>0
 ratio=t*scale[j+1]/scale[j];a,b=ratio.numerator,ratio.denominator
 w=[b*v[j][i]+a*v[j+1][i] for i in range(3)]
 common=gcd(gcd(abs(w[0]),abs(w[1])),abs(w[2]));assert common>0
 v[j]=[x//common for x in w];scale[j]*=Fraction(common,b)
 newsign,newm=evaluate();assert newm<=oldm
 if newsign!=oldsign:flips.append(step)
 if newm!=oldm:drops.append((step,oldm,newm,len(flips)))
 oldsign,oldm=newsign,newm
assert len(data['moves'])==207
assert flips==[15,36,57,78,92,111,132,153,174,195,204,207]
assert drops==[(15,4,3,1),(36,3,2,2),(57,2,1,3),(89,1,0,4)]
print('PASS 207 positive rational moves; all 52 features strict at all 208 vertices')
print('PASS 12 target crossings; one trajectory is sharp in every m=0,1,2,3,4 class')
