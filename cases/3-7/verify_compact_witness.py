#!/usr/bin/env python3
"""Replay a compact 12-flip witness by integer columns and exact scales.

Unlike the exterior-power replay, this program evolves integer representatives
of the columns and takes 3-by-3 determinants directly at each vertex.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import json,sys

if not __debug__:raise RuntimeError('assertions are required; do not run Python with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
if len(sys.argv)!=2:raise SystemExit('usage: verify_compact_witness.py WITNESS.json')
data=json.loads(Path(sys.argv[1]).read_text(),parse_float=Q)
v=[];scale=[]
for column in data['initial_columns']:
    q=list(map(Q,column));assert len(q)==3
    common=lcm(*(x.denominator for x in q))
    z=[int(x*common) for x in q];g=gcd(*(abs(x) for x in z));assert g
    v.append([x//g for x in z]);scale.append(Q(g,common))
assert len(v)==7
triples=list(combinations(range(7),3));quads=list(combinations(range(7),4))
features=[]
for a,b,c,d in combinations(range(1,6),4):
    features.extend([(a,b,c,d),(a,c,b,d),(a,d,b,c)])
def det(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def signs_and_m():
    minors={I:det(*(v[i] for i in I)) for I in triples}
    assert all(minors.values())
    def p(a,b,c):return minors[a-1,b-1,c-1]
    vals=[p(a,b,6)*p(c,d,7)-p(a,b,7)*p(c,d,6) for a,b,c,d in features]
    vals += [p(1,2,5)*p(3,4,j)-p(1,2,j)*p(3,4,5) for j in (6,7)]
    assert all(vals)
    signs={I:(1 if q>0 else -1) for I,q in minors.items()}
    feasible=[]
    for mask in range(64):
        for a,b,c,d in quads:
            circuit=[signs[b,c,d],-signs[a,c,d],signs[a,b,d],-signs[a,b,c]]
            x=[circuit[i]*(-1 if mask>>j&1 else 1) for i,j in enumerate((a,b,c,d))]
            if x[0]==x[1]==x[2]==x[3]:break
        else:feasible.append((mask^(mask>>1)).bit_count())
    assert feasible
    return signs[0,1,2],min(feasible)

old,m=signs_and_m();assert m==4
flips=[];drops=[]
for step,move in enumerate(data['moves'],1):
    if isinstance(move,dict):j=int(move['move'])-1;t=Q(move['time'])
    else:j=int(move[0])-1;t=Q(move[1])
    assert 0<=j<6 and t>0
    ratio=t*scale[j+1]/scale[j];a,b=ratio.numerator,ratio.denominator
    z=[b*v[j][i]+a*v[j+1][i] for i in range(3)]
    g=gcd(*(abs(x) for x in z));assert g
    v[j]=[x//g for x in z];scale[j]*=Q(g,b)
    now,mnew=signs_and_m();assert mnew<=m
    if old!=now:flips.append(step)
    if mnew!=m:drops.append((step,m,mnew,len(flips)))
    old,m=now,mnew
assert len(flips)==12
if isinstance(data.get('target_flip_steps'),list):assert data['target_flip_steps']==flips
print('PASS',len(data['moves']),'positive rational moves; 52 strict features at every vertex')
print('PASS exactly 12 target changes at',flips)
print('PASS minimum-variation drops:',drops)
