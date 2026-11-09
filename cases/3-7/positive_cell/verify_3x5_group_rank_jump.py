#!/usr/bin/env python3
"""Direct 3x5 sign-prerank and five-sign jump-difference certificate.

No prior graph or rank files are read. Every edge and every jump pair is
checked. All 384 states are necessary Pluecker sign patterns, not claimed
realizable.
"""
from itertools import combinations
P=list(combinations(range(1,6),3));idx={I:i for i,I in enumerate(P)}
def sgn(s,a,b,c):
 q=(a,b,c);par=sum(q[i]>q[j] for i in range(3) for j in range(i+1,3));z=1 if s>>idx[tuple(sorted(q))]&1 else -1
 return -z if par%2 else z
def admissible(s):
 for z in range(1,6):
  a,b,c,d=(q for q in range(1,6) if q!=z)
  q=(sgn(s,z,a,b)*sgn(s,z,c,d),-sgn(s,z,a,c)*sgn(s,z,b,d),sgn(s,z,a,d)*sgn(s,z,b,c))
  if len(set(q))==1:return False
 return True
states=[s for s in range(1024) if admissible(s)];S=set(states);assert len(states)==384
arcs={s:[] for s in states}
for s in states:
 for j in range(1,5):
  for i,I in enumerate(P):
   if j not in I or j+1 in I:continue
   J=tuple(sorted(set(I)-{j}|{j+1}));d=idx[J]
   if ((s>>i)^(s>>d))&1:
    u=s^(1<<i)
    if u in S:arcs[s].append((u,int(I==(1,2,3)),j))
edgecount=sum(map(len,arcs.values()));assert edgecount==1152
# Longest future target weight, by relaxation; termination certifies there is
# no reachable positive cycle. Every local inequality is checked separately.
R={s:0 for s in states}
for iteration in range(len(states)+1):
 changed=False
 for s in states:
  q=max([0]+[R[u]+w for u,w,j in arcs[s]])
  if q>R[s]:R[s]=q;changed=True
 if not changed:break
else:raise AssertionError('positive target-weight cycle')
assert all(R[s]>=R[u]+w for s in states for u,w,j in arcs[s])
alt=sum((sum(I)%2==0)<<i for i,I in enumerate(P));assert R[alt]==6 and max(R.values())==6
fixed=[idx[I] for I in combinations(range(1,5),3)]+[idx[(3,4,5)]]
groups={}
for s in states:
 key=tuple((s>>i)&1 for i in fixed);groups.setdefault(key,[]).append(s)
assert len(groups)==32 and all(len(v)==12 for v in groups.values())
pairs=0
for values in groups.values():
 for a in values:
  for b in values:
   assert R[b]-R[a]<=2;pairs+=1
print('PASS',len(states),'necessary sign states,',edgecount,'move-labeled edges, prerank maximum',max(R.values()),'and alternating rank',R[alt])
print('PASS',len(groups),'five-sign groups of size 12;',pairs,'ordered jump pairs all satisfy R_after−R_before<=2')
print('PASS ≤2 U5 moves preserving p345 allow at most 6+2+2=10 target changes')
