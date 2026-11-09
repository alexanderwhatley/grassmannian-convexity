#!/usr/bin/env python3
"""A second, rank-two cyclic-word proof of the 3x5 prerank jump bound."""
from itertools import permutations,product,combinations
N=5
P=list(combinations(range(1,6),3));idx={I:i for i,I in enumerate(P)}
def orient(w,i,j):
 d=(w.index(j)-w.index(i))%(2*N);assert d not in (0,N)
 return 1 if d<N else -1
def rank(w):
 d=[]
 for i in range(1,N):
  t=(w.index(i+1)-w.index(i))%(2*N);assert t not in (0,N)
  if t>N:t-=2*N
  d.append(t)
 runs=[]
 for t in d:
  if not runs or t*runs[-1][0]<0:runs.append([t,abs(t)])
  else:runs[-1][1]+=abs(t)
 return 2*sum(z//N for s,z in runs)+len(runs)-1
G={};seen=set()
for perm in permutations(range(1,6)):
 if perm[0]!=1:continue
 for flips in product((0,1),repeat=4):
  flip={1:0,**{i:flips[i-2] for i in range(2,6)}}
  a=[(-x if flip[x] else x) for x in perm];w=tuple(a+[-x for x in a])
  if orient(w,4,5)<0:continue
  # Under 3x5 Gale duality, each pair sign equals the complementary 3-minor
  # sign after reversal, up to one common positive scalar because p345>0.
  signs={}
  for i,j in combinations(range(1,6),2):
   I=tuple(x for x in range(1,6) if x not in (6-j,6-i));signs[I]=orient(w,i,j)
  assert len(signs)==10 and signs[3,4,5]>0
  key=tuple(signs[I] for I in P);assert key not in seen;seen.add(key)
  pattern=''.join('+' if signs[I]>0 else '-' for I in combinations(range(1,5),3))
  G.setdefault(pattern,[]).append(rank(w))
assert len(seen)==192 and len(G)==16
for pattern,values in sorted(G.items()):
 assert len(values)==12 and max(values)-min(values)==2
 print(pattern,'cyclic words=',len(values),'ranks=',sorted(set(values)))
print('PASS all 192 rank-two cyclic words; each fixed-sign group has rank range exactly 2')
