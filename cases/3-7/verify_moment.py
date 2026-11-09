#!/usr/bin/env python3
"""Independent exact exterior-power replay of the moment-curve 12 example.

Evolves all 35 maximal minors, not the underlying 3x7 columns or fitting
code. Each of 228 moves is strictly positive, every recorded maximal minor
is nonzero, and p123 changes sign exactly twelve times.
"""
from itertools import combinations,permutations
from fractions import Fraction as Q
from pathlib import Path
import json,sys
sys.set_int_max_str_digits(0)
root=Path(__file__).resolve().parent
key=sys.argv[1] if len(sys.argv)>1 else '18204'
path=Path(sys.argv[2]) if len(sys.argv)>2 else root/('moment_curve_12_'+key+'.json')
z=json.loads(path.read_text())
v=[tuple(map(Q,c)) for c in z['initial_columns']]
assert len(v)==7
nodes=[]
for i,c in enumerate(v,1):
 assert c[0]==(-1)**i
 t=c[1]/c[0];assert c[2]/c[0]==t*t;nodes.append(t)
assert all(a<b for a,b in zip(nodes,nodes[1:]))
def det(a,b,c):
 v=(a,b,c);total=Q(0)
 for p in permutations(range(3)):
  inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
  total+=(-1 if inv%2 else 1)*v[0][p[0]]*v[1][p[1]]*v[2][p[2]]
 return total
I=list(combinations(range(7),3));P={a:det(*(v[i] for i in a)) for a in I};assert all(P.values())
for (a,b,c),p in P.items():assert (p>0)==((a+b+c+3)%2==0)
assert P[(4,5,6)]>0
p=lambda a,b,c:P[tuple(sorted((a-1,b-1,c-1)))]
features=[]
for a,b,c,d in combinations(range(1,6),4):
 for i,j,k,l in ((a,b,c,d),(a,c,b,d),(a,d,b,c)):
  features.append(p(i,j,6)*p(k,l,7)-p(i,j,7)*p(k,l,6))
assert len(features)==15 and all(features)
fword=sum((x>0)<<i for i,x in enumerate(features));assert fword==int(key)
g56=p(1,2,5)*p(3,4,6)-p(1,2,6)*p(3,4,5)
g57=p(1,2,5)*p(3,4,7)-p(1,2,7)*p(3,4,5)
assert g56<0<g57 and features[0]<0
ratios=[p(1,2,i)/p(3,4,i) for i in (5,6,7)];assert ratios[0]>ratios[1]>ratios[2]>0
old=P[(0,1,2)];flips=[]
for n,move in enumerate(z['moves'],1):
 j=move['move']-1;t=Q(move['time']);assert 0<=j<6 and t>0
 changes={}
 for a in I:
  if j in a and j+1 not in a:
   b=tuple(sorted((set(a)-{j})|{j+1}));changes[a]=P[a]+t*P[b]
 P.update(changes);assert all(P.values()),n
 cur=P[(0,1,2)]
 if (old>0)!=(cur>0):flips.append(n)
 old=cur
 if n%50==0:print('steps',n,'flips',len(flips),flush=True)
assert len(z['moves'])>=12 and len(flips)==12 and flips==z['target_flip_steps']
assert all(i>z['prefix_moves'] for i in flips)
print('PASS Fword',fword,'G56<0<G57, r5>r6>r7>0; exact',len(z['moves']),'positive moves from moment-curve nodes',nodes,flush=True)
print('PASS all 35 minors strict at all vertices; exactly 12 target changes at',flips,flush=True)
