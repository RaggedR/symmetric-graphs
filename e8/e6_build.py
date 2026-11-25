"""
Build the E6 rung of the ladder: the 27 lines (minuscule weights), their
Coxeter-plane projection, and the caches the E6 drawers read.

Minuscule descent: the E6 weights are the E7 weights at inner product +1 to
the fixed weight w0 = e1+e3.  The E6 root system is the E8 roots orthogonal to
BOTH r0 = e1+e2 and w0 (72 of them -- the ad-hoc "reflections preserving the
27-set" test over-counts to 74, do not use it); its Coxeter element has order
h = 12.

Writes /tmp/e6_W.npy (the 27 weights) and /tmp/e6_P.npy (their projection).
Run after e7_build.py (it reads /tmp/e7_W.npy), or it rebuilds E7 itself.
"""
import numpy as np
import ladder as L

r0, w0 = L.R0, L.W0

# the 56 E7 weights (reuse the cache if present, else rebuild)
try:
    W7 = np.load("/tmp/e7_W.npy")
    assert len(W7) == 56
except (FileNotFoundError, AssertionError):
    R = L.e8_roots()
    W7 = R[np.round(R @ r0).astype(int) == 1]

# the 27 E6 weights = E7 weights at ip +1 to w0
W = W7[np.round(W7 @ w0).astype(int) == 1]
assert len(W) == 27, len(W)

# the E6 root subsystem = E8 roots orthogonal to both r0 and w0
R = L.e8_roots()
E6 = R[(np.abs(R @ r0) < 1e-9) & (np.abs(R @ w0) < 1e-9)]
assert len(E6) == 72, len(E6)
order, v1, v2 = L.coxeter_plane(L.simple_roots(E6), 12)
P = L.project(W, v1, v2)

rings = L.ring_histogram(P)
print("E6: 27 weights, Coxeter order", order)
print("ring histogram (radius: count):", rings)
# expected: 3 central points (radius 0) + two rings of 12
counts = sorted(rings.values())
assert counts == [3, 12, 12], f"unexpected rings {counts} (want 3 + 2x12)"
print("ring structure 3 + 2x12 verified.")

np.save("/tmp/e6_W.npy", W)
np.save("/tmp/e6_P.npy", P)
print("wrote /tmp/e6_W.npy, /tmp/e6_P.npy")
