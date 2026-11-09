#!/usr/bin/env python3
"""Check an explicit two-valued approximation to the 3x5 graph rank."""
from pathlib import Path
import runpy
G = runpy.run_path(str(Path(__file__).with_name('verify_3x5_group_rank_jump.py')))
chain = [(1,2,3),(1,2,4),(1,3,4),(2,3,4),(3,4,5)]
counts = [0,0]
for state in G['states']:
    bits = [(state >> G['idx'][I]) & 1 for I in chain]
    V = sum(a != b for a,b in zip(bits,bits[1:]))
    d = G['R'][state] - V
    assert d in (0,2), (state,V,G['R'][state])
    counts[d//2] += 1
assert sum(counts)==384
print('PASS 3x5 graph rank is V or V+2 on all 384 states; counts',counts)
