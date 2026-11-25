import numpy as np, itertools

# ---------- 1. the 240 roots of E8 in R^8 (standard Euclidean inner product) ----------
roots = []
# type A: (+-1,+-1,0^6), all positions & signs  -> C(8,2)*4 = 112
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = [0]*8; v[i] = si; v[j] = sj; roots.append(v)
# type B: (+-1/2)^8 with an EVEN number of minus signs -> 2^7 = 128
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append([s*0.5 for s in signs])
R = np.array(roots, dtype=float)
print("number of roots:", len(R), "(expect 240)")
print("squared lengths all == 2:", np.allclose((R*R).sum(1), 2.0))

# inner-product signature from a fixed root
g = R @ R[0]
from collections import Counter
sig = Counter(np.round(g, 6))
print("inner-product histogram from one root:", dict(sorted(sig.items())))
print("  expect {-2:1, -1:56, 0:126, 1:56, 2:1}")

# ---------- 2. Coxeter element & Coxeter plane ----------
# Bourbaki simple roots of E8
e = np.eye(8)
a1 = 0.5*(e[0] - e[1] - e[2] - e[3] - e[4] - e[5] - e[6] + e[7])
a2 = e[0] + e[1]
a3 = e[1] - e[0]
a4 = e[2] - e[1]
a5 = e[3] - e[2]
a6 = e[4] - e[3]
a7 = e[5] - e[4]
a8 = e[6] - e[5]
simple = [a1, a2, a3, a4, a5, a6, a7, a8]
print("simple roots all length^2==2:", all(abs(a@a-2)<1e-9 for a in simple))

def refl(alpha):
    # s_alpha(x) = x - (x.alpha) alpha   (since alpha.alpha = 2)
    return np.eye(8) - np.outer(alpha, alpha)

C = np.eye(8)
for a in simple:            # Coxeter element = product of the 8 simple reflections
    C = C @ refl(a)
# order of C
Cn = np.eye(8); order = 0
for k in range(1, 100):
    Cn = Cn @ C
    if np.allclose(Cn, np.eye(8), atol=1e-8): order = k; break
print("order of Coxeter element:", order, "(expect 30 = Coxeter number of E8)")

# eigenvector for eigenvalue e^{2 pi i /30}  -> the Coxeter plane
w, V = np.linalg.eig(C)
theta = 2*np.pi/30
idx = np.argmin(np.abs(np.angle(w) - theta))
v = V[:, idx]
u1 = np.real(v); u2 = np.imag(v)
# orthonormalise (Gram-Schmidt) for a clean projection
u1 = u1/np.linalg.norm(u1)
u2 = u2 - (u2@u1)*u1; u2 = u2/np.linalg.norm(u2)
print("Coxeter-plane basis orthonormal:", abs(u1@u2)<1e-9, np.allclose([u1@u1,u2@u2],1))

# ---------- 3. project & verify 8 rings of 30 ----------
P = np.column_stack([R @ u1, R @ u2])
radii = np.round(np.hypot(P[:,0], P[:,1]), 4)
ring_hist = Counter(radii)
print("\nprojected-radius histogram (radius : count):")
for r in sorted(ring_hist): print(f"   {r:.4f} : {ring_hist[r]}")
print("number of distinct rings:", len(ring_hist), "(expect 8)")
print("every ring has 30 points:", all(c == 30 for c in ring_hist.values()))

# verify 30-fold: within each ring the angles are spaced by 12 deg
ok = True
for r in sorted(ring_hist):
    ang = np.sort(np.degrees(np.arctan2(P[radii==r,1], P[radii==r,0])) % 360)
    gaps = np.diff(np.concatenate([ang, [ang[0]+360]]))
    if not np.allclose(gaps, 12.0, atol=1e-6): ok = False
print("every ring is 30-fold symmetric (12 deg spacing):", ok)

np.save('/tmp/e8_P.npy', P)
np.save('/tmp/e8_R.npy', R)
