#!/usr/bin/env python3
"""Exact over-approximation for (3,5) using just ten minor signs."""
from itertools import combinations,permutations
from random import Random

n=5;k=3
P=list(combinations(range(1,n+1),k));idx={I:i for i,I in enumerate(P)}
def ori(a,b,c):
    x=(a,b,c);return (-1 if sum(x[i]>x[j] for i in range(3) for j in range(i+1,3))%2 else 1),idx[tuple(sorted(x))]
def sg(key,I):return 1 if key>>idx[I]&1 else -1
def gos(key,*I):
    s,i=ori(*I);return s*(1 if key>>i&1 else -1)
def admissible(key):
    if sg(key,(3,4,5))<0:return False
    for s in range(1,6):
        a,b,c,d=[x for x in range(1,6) if x!=s]
        z=[gos(key,s,a,b)*gos(key,s,c,d),-gos(key,s,a,c)*gos(key,s,b,d),gos(key,s,a,d)*gos(key,s,b,c)]
        if abs(sum(z))==3:return False
    return True
def mvariation(key):
    values=[]
    for mask in range(1<<(n-1)):
        signs=[-1 if mask>>i&1 else 1 for i in range(n-1)]+[1]
        for a,b,c,d in combinations(range(1,6),4):
            circ=[sg(key,(b,c,d)),-sg(key,(a,c,d)),sg(key,(a,b,d)),-sg(key,(a,b,c))]
            z=[circ[j]*signs[i-1] for j,i in enumerate((a,b,c,d))]
            if len(set(z))==1:break
        else:values.append(sum(a!=b for a,b in zip(signs,signs[1:])))
    return min(values) if values else 100
states=[a for a in range(1<<10) if admissible(a)];S=set(states)
arcs={a:[] for a in states}
for key in states:
    for j in range(1,5):
        for i,I in enumerate(P):
            if j in I and j+1 not in I:
                J=tuple(sorted((set(I)-{j})|{j+1}));d=sg(key,J);s=sg(key,I)
                if s!=d and key^(1<<i) in S:arcs[key].append((key^(1<<i),int(I==(1,2,3)),j,I))
rank={a:0 for a in states}
for turn in range(len(states)+1):
    changed=False
    for a in states:
        x=max([rank[b]+w for b,w,j,I in arcs[a]]+[0])
        if x>rank[a]:rank[a]=x;changed=True
    if not changed:break
else:raise AssertionError('positive cycle in over-approximation')
maxima={}
for a in states:
    m=mvariation(a);assert 0<=m<=2
    maxima[m]=max(maxima.get(m,0),rank[a])
assert maxima=={0:4,1:5,2:6},maxima
assert all(rank[a]<=4+mvariation(a) for a in states)
print('PASS',len(states),'necessary sign states;',sum(map(len,arcs.values())),'edges; max rank by m',maxima)

def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def vals(V):return [det(*(V[i-1] for i in I)) for I in P]
rng=Random(351)
for trial in range(1000):
    v=[[rng.randrange(-20,21) for _ in range(3)] for _ in range(5)];a=vals(v)
    if not all(a):continue
    ok=True
    for j in range(1,5):
        w=[x[:] for x in v];w[j-1]=[x+y for x,y in zip(w[j-1],w[j])];b=vals(w)
        active=[i for i,I in enumerate(P) if j in I and j+1 not in I]
        d=[b[i]-a[i] for i in range(10)]
        if any(not d[i] for i in active) or any(a[i]*d[h]-a[h]*d[i]==0 for i,h in combinations(active,2)):
            ok=False;break
    if ok:
        print('PASS exact generic event separation witness',v);break
else:raise AssertionError('no separation witness')
