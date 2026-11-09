#!/usr/bin/env python3
"""Independently replay any output of the exact rational IFT constructor."""
from fractions import Fraction as Q
from itertools import combinations,permutations
from pathlib import Path
import json,sys
if not __debug__:raise RuntimeError('assertions are required; do not run Python with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
z=json.loads(Path(sys.argv[1]).read_text());v=[tuple(map(Q,c)) for c in z['initial_columns']]
assert len(v)==7 and all(len(c)==3 for c in v)
def det(a,b,c):
 cols=(a,b,c);total=Q(0)
 for p in permutations(range(3)):
  inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
  total+=(-1 if inv%2 else 1)*cols[0][p[0]]*cols[1][p[1]]*cols[2][p[2]]
 return total
I=list(combinations(range(7),3));P={a:det(*(v[i] for i in a)) for a in I};assert all(P.values())
orient=1 if P[(4,5,6)]>0 else -1
for (a,b,c),p in P.items():assert (p>0)==(orient*((-1)**(a+b+c+3))>0)
old=P[(0,1,2)];flips=[]
for step,move in enumerate(z['moves'],1):
 j=move['move']-1;t=Q(move['time']);assert 0<=j<6 and t>0
 changed={}
 for a in I:
  if j in a and j+1 not in a:
   b=tuple(sorted((set(a)-{j})|{j+1}));changed[a]=P[a]+t*P[b]
 P.update(changed);assert all(P.values()),step
 now=P[(0,1,2)]
 if (old>0)!=(now>0):flips.append(step)
 old=now
assert len(z['moves'])>=24 and z['prefix_moves']==12 and len(flips)==12
assert all(i>12 for i in flips)
if 'target_flip_steps' in z:assert z['target_flip_steps']==flips
print('PASS generic alternating input;',len(z['moves']),'positive rational moves; 35 minors strict at every vertex')
print('PASS exactly 12 target changes at',flips,'with no change in the first 12 moves')
