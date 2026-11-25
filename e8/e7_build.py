"""
Build the E7 rung of the ladder: the 56 minuscule weights, their Coxeter-plane
projection, and the caches the E7 drawers read.

Minuscule descent: the E7 weights are the E8 roots at inner product +1 (60 deg)
to the fixed root r0 = e1+e2.  The E7 root system itself is the E8 roots
orthogonal to r0 (126 of them); its Coxeter element has order h = 18.

Writes /tmp/e7_W.npy (the 56 weights), /tmp/e7_P.npy (their projection),
and /tmp/e7_r0.npy (the descent vector).  Run after e8_mandala.py (or it
rebuilds the E8 roots itself).
"""
import numpy as np
import ladder as L

R = L.e8_roots()
r0 = L.R0

# the 56 minuscule weights = roots at ip +1 to r0
W = R[np.round(R @ r0).astype(int) == 1]
assert len(W) == 56, len(W)

# the E7 root subsystem = roots orthogonal to r0, and its Coxeter plane
E7 = R[np.abs(R @ r0) < 1e-9]
assert len(E7) == 126, len(E7)
order, u1, u2 = L.coxeter_plane(L.simple_roots(E7), 18)
P = L.project(W, u1, u2)

rings = L.ring_histogram(P)
print("E7: 56 weights, Coxeter order", order)
print("ring histogram (radius: count):", rings)
# expected: 2 central points (radius 0) + three rings of 18
counts = sorted(rings.values())
assert counts == [2, 18, 18, 18], f"unexpected rings {counts} (want 2 + 3x18)"
print("ring structure 2 + 3x18 verified.")

np.save("/tmp/e7_W.npy", W)
np.save("/tmp/e7_P.npy", P)
np.save("/tmp/e7_r0.npy", r0)
print("wrote /tmp/e7_W.npy, /tmp/e7_P.npy, /tmp/e7_r0.npy")
