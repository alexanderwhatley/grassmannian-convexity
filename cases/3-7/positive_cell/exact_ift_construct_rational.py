#!/usr/bin/env python3
"""Construct a positive rational 12-crossing path from rational generic m=4 M.

Usage:
  python3 exact_ift_construct_rational.py MATRIX.json WITNESS.json OUTPUT.json

MATRIX.json is a list of seven 3-component columns, or an object with an
"initial_columns" or "columns" field. WITNESS.json is a generic
alternating-chamber 12-crossing witness. All arithmetic and all 35-minor
checks are exact.
The accompanying proof explains why the geometric epsilon search terminates.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json,sys,time
if not __debug__:raise RuntimeError('assertions are required; do not run Python with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=Path(__file__).resolve().parent
if len(sys.argv)!=4:
 raise SystemExit('usage: exact_ift_construct_rational.py MATRIX.json WITNESS.json OUTPUT.json')
raw=json.loads(Path(sys.argv[1]).read_text(),parse_float=Q)
if isinstance(raw,dict):raw=raw.get('initial_columns',raw.get('columns'))
if not isinstance(raw,list) or len(raw)!=7 or any(not isinstance(c,list) or len(c)!=3 for c in raw):
 raise ValueError('input must contain seven 3-component columns')
original=[list(map(Q,c)) for c in raw]
def det_initial(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
d567=det_initial(*original[4:7])
if not d567:raise ValueError('p567 must be nonzero')
seed=[([ -c[0],c[1],c[2] ] if d567<0 else c[:]) for c in original]
witness=json.loads(Path(sys.argv[2]).read_text(),parse_float=Q)
W=[list(map(Q,c)) for c in witness['initial_columns']]
if len(W)!=7 or any(len(c)!=3 for c in W):raise ValueError('witness must contain seven 3-component columns')
w567=det_initial(*W[4:7])
if not w567:raise ValueError('witness p567 must be nonzero')
if w567<0:W=[[-c[0],c[1],c[2]] for c in W]
witness_moves=[]
for move in witness['moves']:
 if isinstance(move,dict):j,t=move['move'],Q(move['time'])
 else:j,t=move[0],Q(move[1])
 if not isinstance(j,int) or not 1<=j<=6 or t<=0:raise ValueError('witness moves must be positive and have generator in 1..6')
 witness_moves.append((j,t))
I=list(combinations(range(7),3));D=[(-1)**i for i in range(1,8)];word=[4,3,2,1,5,4,3,2,6,5,4,3]
def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def pvals(cols):return {a:det(*(cols[i] for i in a)) for a in I}
def inv(A):
 n=len(A);assert all(len(row)==n for row in A)
 B=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
 for c in range(n):
  r=next(i for i in range(c,n) if B[i][c]);B[c],B[r]=B[r],B[c];q=B[c][c];B[c]=[x/q for x in B[c]]
  for i in range(n):
   if i!=c and B[i][c]:q=B[i][c];B[i]=[x-q*y for x,y in zip(B[i],B[c])]
 return [row[n:] for row in B]
def chart(cols):
 B=[[cols[4+j][i] for j in range(3)] for i in range(3)];C=inv(B)
 return [[sum(C[i][r]*cols[c][r] for r in range(3)) for c in range(4)] for i in range(3)]
Mbar=[[D[i]*x for x in col] for i,col in enumerate(seed)];Wbar=[[D[i]*x for x in col] for i,col in enumerate(W)]
p0=pvals(Mbar)
if not all(q>0 for q in p0.values()):raise ValueError('input is not generic with m=4 / alternating Pluecker signs')
if not all(q>0 for q in pvals(Wbar).values()):raise ValueError('witness must start in the generic alternating chamber')
if not witness_moves:raise ValueError('witness must contain positive moves')
norm=p0[(4,5,6)]
p=lambda a,b,c:p0[(a-1,b-1,c-1)]/norm
T0=[p(1,6,7)*p(4,5,6)/p(1,5,6),p(1,5,6)*p(3,4,5)/(p(1,4,5)*p(4,5,6)),
 p(1,4,5)*p(2,3,4)/(p(1,3,4)*p(3,4,5)),p(1,3,4)/p(2,3,4),
 p(1,2,7)*p(1,5,6)/(p(1,2,6)*p(1,6,7)),p(1,2,6)*p(1,4,5)/(p(1,2,5)*p(1,5,6)),
 p(1,2,5)*p(1,3,4)/(p(1,2,4)*p(1,4,5)),p(1,2,4)/p(1,3,4),
 p(1,2,6)/p(1,2,7),p(1,2,5)/p(1,2,6),p(1,2,4)/p(1,2,5),p(1,2,3)/p(1,2,4)]
assert len(T0)==12 and all(t>0 for t in T0)
X=chart(Mbar);Wc=chart(Wbar)
# G = Z_first - Z_last*X where Z=[εWc|I] T(t), in row-major order.
def gvalue(t):
 Z=[[Q(0) for c in range(7)] for r in range(3)]
 for r in range(3):Z[r][4+r]=1
 for j,a in zip(word,t):
  for r in range(3):Z[r][j-1]+=a*Z[r][j]
 return [Z[r][c]-sum(Z[r][4+k]*X[k][c] for k in range(3)) for r in range(3) for c in range(4)]
g0=gvalue(T0);assert not any(g0)
# G is affine in each individual t_i, so unit finite differences are exact.
J=[]
for j in range(12):
 tt=T0[:];tt[j]+=1;J.append(gvalue(tt))
J=[[J[c][r] for c in range(12)] for r in range(12)];Jinv=inv(J)
# Truncated Taylor arithmetic to degree three in ε.
L=4
def add(a,b):return [a[i]+b[i] for i in range(L)]
def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(L)]
def muls(a,b):return [q*b for q in a]
def gseries(T):
 Z=[[[Q(0) for i in range(L)] for c in range(7)] for r in range(3)]
 for r in range(3):
  for c in range(4):Z[r][c][1]=Wc[r][c]
  Z[r][4+r][0]=1
 for j,a in zip(word,T):
  for r in range(3):Z[r][j-1]=add(Z[r][j-1],mul(a,Z[r][j]))
 out=[]
 for r in range(3):
  for c in range(4):
   q=Z[r][c][:]
   for a in range(3):q=add(q,muls(Z[r][4+a],-X[a][c]))
   out.append(q)
 return out
start=time.perf_counter();T=[[t,Q(0),Q(0),Q(0)] for t in T0]
for degree in (1,2,3):
 G=gseries(T);known=[-g[degree] for g in G]
 solution=[sum(Jinv[i][j]*known[j] for j in range(12)) for i in range(12)]
 for i in range(12):T[i][degree]=solution[i]
 assert all(g[degree]==0 for g in gseries(T))
 print('PASS exact Taylor coefficient',degree,'elapsed',round(time.perf_counter()-start,3),flush=True)
print('max Taylor coefficient bits',max(max(abs(q.numerator).bit_length(),q.denominator.bit_length()) for a in T for q in a),flush=True)
# Independent-style minor-only update for witness sign templates and candidates.
def apply(P,j,t):
 j-=1;updates={}
 for a in I:
  if j in a and j+1 not in a:
   b=tuple(sorted((set(a)-{j})|{j+1}));updates[a]=P[a]+t*P[b]
 P.update(updates)
def signmask(P):
 if not all(P.values()):return None
 return sum((P[a]>0)<<i for i,a in enumerate(I))
P=pvals(W);wmasks=[signmask(P)];assert wmasks[0] is not None
for j,t in witness_moves:
 apply(P,j,t);wmasks.append(signmask(P));assert wmasks[-1] is not None
assert sum(((a^b)&1)!=0 for a,b in zip(wmasks,wmasks[1:]))==12
initial=pvals(seed);imask=signmask(initial);assert imask is not None
assert (imask&1)==(wmasks[0]&1)
power=1
while True:
 eps=Q(1,2**power);tt=[sum(a[r]*eps**r for r in range(4)) for a in T]
 if not all(t>0 for t in tt):
  print('epsilon exponent',power,'nonpositive parameter',flush=True);power*=2;continue
 P=initial.copy();prefix=[];reason=None
 for j,t in reversed(list(zip(word,tt))):
  apply(P,j,t);prefix.append({'move':j,'time':str(t)})
  if signmask(P)!=imask:reason='prefix sign mismatch';break
 if reason:
  print('epsilon exponent',power,reason,flush=True);power*=2;continue
 moves=prefix[:];flips=[];old=P[(0,1,2)]
 for index,(j,s) in enumerate(witness_moves,1):
  t=s*(eps if j==4 else 1);apply(P,j,t);mask=signmask(P)
  if mask!=wmasks[index]:reason='witness sign mismatch at '+str(index);break
  cur=P[(0,1,2)]
  if (old>0)!=(cur>0):flips.append(12+index)
  old=cur;moves.append({'move':j,'time':str(t)})
 if reason:
  print('epsilon exponent',power,reason,flush=True);power*=2;continue
 assert len(moves)==12+len(witness_moves) and len(flips)==12
 out={'initial_columns':[[str(q) for q in col] for col in original], 'moves':moves,'target_changes':12,
      'target_flip_steps':flips,'prefix_moves':12,'epsilon':str(eps),
      'description':'Exact rational cubic-IFT prefix plus scaled 12-crossing witness'}
 out.pop('seed_id',None)
 path=Path(sys.argv[3]);path.write_text(json.dumps(out,indent=2)+'\n')
 bits=max(max(abs(Q(m['time']).numerator).bit_length(),Q(m['time']).denominator.bit_length()) for m in moves)
 print('PASS exact rational construction','epsilon exponent',power,'moves',len(moves),'flips',flips,'max parameter bits',bits,'elapsed',round(time.perf_counter()-start,3),flush=True)
 break

