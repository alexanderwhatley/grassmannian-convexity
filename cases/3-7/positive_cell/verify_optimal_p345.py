#!/usr/bin/env python3
"""Check a necessary sign condition for every shortest 12-crossing path.

The neighboring independent checker reconstructs the necessary five-column
sign graph and verifies its event-separation witness. This script performs a
separate reachability search in that graph. It deliberately allows arbitrary
U5 jumps preserving the first-four minors and p345, and arbitrary minor-event
closures in a physical U1, U2, U3, or U4 move.
"""
from collections import deque
from pathlib import Path
import runpy

base = runpy.run_path(str(Path(__file__).with_name('verify_minimum_28_independent.py')))
states, arcs, idx = base['states'], base['arcs'], base['idx']
groups, group, start = base['groups'], base['group'], base['start']
b123, b124, b345 = (1 << idx[x] for x in ((1,2,3),(1,2,4),(3,4,5)))

def closure(source, generator):
    seen = {source}
    todo = [source]
    while todo:
        v = todo.pop()
        for w, target, j in arcs[v]:
            if j == generator and w not in seen:
                seen.add(w)
                todo.append(w)
    return tuple(seen)

closures = {(s, j): closure(s, j) for s in states for j in range(1, 5)}
# phase: 0 = require a target-changing U3 next; 1 = require a U2/U4
# separator before the next U3. U1/U5 may be inserted in either phase.
source = (start, 0, 0, 0, 0)
seen = {source}
queue = deque([source])
best = 0
while queue:
    s, nflip, phase, n1, n5 = queue.popleft()
    best = max(best, nflip)
    assert nflip < 12, 'an abstract 12-flip path violates the asserted bound'
    candidates = []
    if phase == 0:
        candidates.extend((t, nflip+1, 1, n1, n5) for t in closures[s, 3]
                          if (s ^ t) & b123)
    if phase == 1 and nflip < 12:
        for j in (2, 4):
            candidates.extend((t, nflip, 0, n1, n5) for t in closures[s, j]
                              if (s ^ t) & b124)
    if n1 < 2:
        candidates.extend((t, nflip, phase, n1+1, n5) for t in closures[s, 1])
    if n5 < 2:
        candidates.extend((t, nflip, phase, n1, n5+1) for t in groups[group(s)]
                          if not ((s ^ t) & b345))
    for q in candidates:
        if q not in seen:
            seen.add(q)
            queue.append(q)
assert best < 12
print(f'PASS no 12-flip path with at most two U1 and two U5 moves if p345 is unchanged; reachable={len(seen)}, max_flips={best}')
