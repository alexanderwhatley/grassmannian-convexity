#!/usr/bin/env python3
"""Independent finite verifier for the <=2 U5, one-U1-run obstruction.

No prior graph code or rank data is imported. States are ten Pluecker signs,
with all five necessary 3-term Grassmann-Pluecker sign relations imposed.
"""
from collections import deque,Counter
from itertools import combinations
from random import Random

P=list(combinations(range(1,6),3));idx={I:i for i,I in enumerate(P)}
def bit(s,I):return (s>>idx[tuple(sorted(I))])&1
def oriented(s,a,b,c):
 q=(a,b,c);par=sum(q[i]>q[j] for i in range(3) for j in range(i+1,3))
 return (-1 if par%2 else 1)*(1 if bit(s,q) else -1)
def necessary(s):
 for z in range(1,6):
  a,b,c,d=(x for x in range(1,6) if x!=z)
  q=(oriented(s,z,a,b)*oriented(s,z,c,d),
    -oriented(s,z,a,c)*oriented(s,z,b,d),
     oriented(s,z,a,d)*oriented(s,z,b,c))
  if len(set(q))==1:return False
 return True
states=[s for s in range(1<<10) if necessary(s)]
assert len(states)==384
S=set(states);fixed=[idx[x] for x in combinations(range(1,5),3)]
groups={}
for s in states:
 q=sum(((s>>j)&1)<<i for i,j in enumerate(fixed))
 groups.setdefault(q,[]).append(s)
assert len(groups)==16
arcs={s:[] for s in states};baseedges=0
for s in states:
 for j in range(1,5):
  for I in P:
   if j not in I or j+1 in I:continue
   successor=tuple(sorted((set(I)-{j})|{j+1}))
   if bit(s,I)==bit(s,successor):continue
   b=s^(1<<idx[I])
   if b not in S:continue
   arcs[s].append((b,int(I==(1,2,3)),j));baseedges+=1
assert baseedges==1152

def group(s):return sum(((s>>j)&1)<<i for i,j in enumerate(fixed))
def rankstate(s,used,mode):return (s,used,mode)
start=sum(int((-1)**sum(I)>0)<<idx[I] for I in P)
assert start in S
source=rankstate(start,0,0);dist={source:0};queue=deque([source]);inq={source}
while queue:
 s,used,mode=queue.popleft();inq.remove((s,used,mode));score=dist[s,used,mode]
 possible=[]
 # mode 0=before the sole U1 run, 1=inside, 2=after it
 if mode<2:possible.append((s,used,2,0))
 for b,w,j in arcs[s]:
  if j==1 and mode!=2:possible.append((b,used,1,w))
  elif j!=1 and mode!=1:possible.append((b,used,mode,w))
 if used<2 and mode!=1:
  for b in groups[group(s)]:possible.append((b,used+1,mode,0))
 for b,nb,nm,w in possible:
  dest=(b,nb,nm);value=score+w
  assert value<=18,'positive cycle or unexpected rank-two bound violation'
  if value>dist.get(dest,-1):
   dist[dest]=value
   if dest not in inq:queue.append(dest);inq.add(dest)
maxima=[max((v for (s,u,m),v in dist.items() if u==i),default=-1) for i in range(3)]
assert maxima==[6,8,10],maxima
print('PASS independent 384-state,',baseedges,'labeled edges; <=2 U5 and one U1 run maxima',maxima,'reached',len(dist),'automaton states')

# A fixed witness that all same-generator active-root resultants are nonzero.
v=[[19,16,2],[0,-3,0],[-2,19,9],[4,-6,5],[9,-6,-3]]
def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
p=[det(*(v[i-1] for i in I)) for I in P];assert all(p)
ct=0
for j in range(1,5):
 w=[x[:] for x in v];w[j-1]=[a+b for a,b in zip(w[j-1],w[j])]
 q=[det(*(w[i-1] for i in I)) for I in P];d=[q[i]-p[i] for i in range(10)]
 active=[idx[I] for I in P if j in I and j+1 not in I]
 assert len(active)==3 and all(d[i] for i in active)
 for a,b in combinations(active,2):assert p[a]*d[b]-p[b]*d[a]!=0;ct+=1
assert ct==12
print('PASS exact witness separates all 12 same-generator minor-root pairs')
