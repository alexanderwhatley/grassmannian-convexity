#!/usr/bin/env python3
"""Check the elementary all-(k,n) rectangular-chart proof on finite examples."""
import itertools
import random
import sympy as s


def check_paths(n,k):
    r=n-k
    word=[r+i-j for i in range(1,k+1) for j in range(1,r+1)]
    # For each occupied set: number of monomials and a witness mask.
    states={tuple(range(r+1,n+1)):(1,0)}
    for z,j in enumerate(word):
        nxt=dict(states)
        for occupied,(count,mask) in states.items():
            if j+1 in occupied and j not in occupied:
                dest=tuple(sorted((set(occupied)-{j+1})|{j}))
                old,oldmask=nxt.get(dest,(0,0))
                nxt[dest]=(old+count,oldmask if old else mask|(1<<z))
        states=nxt
    assert len(states)==s.binomial(n,k)
    special=0
    for h in range(k+1):
        for l in range(r+1):
            dest=tuple(list(range(1,h+1))+list(range(h+l+1,k+l+1)))
            count,mask=states[dest]
            expected=0
            for i in range(1,k+1):
                for j in range(1,r+1):
                    if i<=h or j<=r-l:
                        expected|=1<<((i-1)*r+j-1)
            assert count==1 and mask==expected,(n,k,h,l,count,mask,expected)
            special+=1
    return special,len(states)

specials=states=0
for n in range(2,15):
    for k in range(1,n):
        a,b=check_paths(n,k);specials+=a;states+=b
print('PASS exactly one claimed monomial:',specials,'solid minors, n<=14')
print('PASS a nonzero path to every one of',states,'positive chart minors')

rng=random.Random(332701)
examples=minorchecks=0
for n in range(2,9):
  for k in range(1,n):
    r=n-k
    subsets=list(itertools.combinations(range(n),k))
    while True:
        X=s.Matrix(k,n,lambda i,j:rng.randrange(-12,13))
        if all(X[:,I].det()!=0 for I in subsets):break
    X=X[:,:k].inv()*X
    B=X[:,k:]
    def C(h,l):
        I=list(range(h))+list(range(h+l,k+l))
        assert len(I)==k
        return X[:,I].det()
    recovered=s.zeros(k,r)
    for h in range(k-1,-1,-1):
        for l in range(1,r+1):
            trial=s.eye(k).row_join(recovered)
            I=list(range(h))+list(range(h+l,k+l))
            det0=trial[:,I].det()
            trial[h,k+l-1]=1
            coeff=trial[:,I].det()-det0
            assert coeff!=0
            recovered[h,l-1]=s.cancel((C(h,l)-det0)/coeff)
            assert recovered[h,l-1]==B[h,l-1]
            q=min(k-h,l)
            Q=B[h:h+q,l-q:l]
            assert Q.det()!=0
            minorchecks+=1
    assert recovered==B
    examples+=1
print('PASS cofactor reconstruction:',minorchecks,'entries across',examples,'generic matrices')

# Chart inverse and all border monomials, checked without symbolic determinants.
for n in range(2,10):
  for k in range(1,n):
    r=n-k
    t=[[s.Rational(rng.randrange(1,10),rng.randrange(1,10)) for j in range(r)] for i in range(k)]
    L=s.eye(n)
    for i in range(1,k+1):
      for j in range(1,r+1):
        gen=r+i-j
        L[:,gen-1]=L[:,gen-1]+t[i-1][j-1]*L[:,gen]
    X=L[r:n,:]
    def C(h,l):
      I=list(range(h))+list(range(h+l,k+l))
      return X[:,I].det()
    for h in range(k+1):
      for l in range(r+1):
        expected=s.prod(t[i][j] for i in range(k) for j in range(r) if i<h or j<r-l)
        assert C(h,l)==expected
    for h in range(k):
      for l in range(1,r+1):
        inverse=C(h,l-1)*C(h+1,l)/(C(h,l)*C(h+1,l-1))
        assert inverse==t[h][r-l]
print('PASS exact border monomials and rational inverse at all n<=9 samples')

# Source-free maximal sign-variation characterization: check parity and
# alternating-subsequence circuit-sign arguments on every small sign sequence.
paritychecks=sequences=0
for n in range(2,11):
  for k in range(1,n):
    r=n-k
    for J in itertools.combinations(range(1,n+1),k-1):
      remain=[i for i in range(1,n+1) if i not in J]
      values={(i+sum(j>i for j in J)-(s+1)-(k-1))%2 for s,i in enumerate(remain)}
      assert values=={0}
      paritychecks+=1
    for bits in range(1<<n):
      sig=[1 if bits>>i&1 else -1 for i in range(n)]
      variation=sum(sig[i]!=sig[i-1] for i in range(1,n))
      if variation>r-1:continue
      tau=[((-1)**(i+1))*sig[i] for i in range(n)]
      assert sum(tau[i]!=tau[i-1] for i in range(1,n))>=k
      selected=[1]
      for i in range(1,n):
        if tau[i]!=tau[i-1] and len(selected)<k+1:selected.append(i+1)
      assert len(selected)==k+1
      circuit=[((-1)**(a+1))*((-1)**(sum(selected)-idx)) for a,idx in enumerate(selected)]
      assert len({circuit[a]*sig[idx-1] for a,idx in enumerate(selected)})==1
      sequences+=1
print('PASS maximal-variation parity:',paritychecks,'J choices and',sequences,'forbidden low-variation sign sequences')
print('All finite checks passed.')
