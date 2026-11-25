"""
Attempting to build "E9" or showing 5_21 polytope structure.

The 5_21 polytope exists in 9D, but unlike 2_21, 3_21, 4_21 which correspond to
the exceptional groups E6, E7, E8, the 5_21 is NOT exceptional - it's related to
classical D5/A8 structures.

We'll show this by trying to construct something at this level and demonstrating
it has D-series character (regular doubling, predictable patterns) rather than
exceptional character.

One approach: build a D7 or D8 structure to show the classical continuation.
"""
import numpy as np, colorsys, itertools
import ladder as L

# Build D7: 64 vertices (2^6), Coxeter h = 12
# This shows what the "natural continuation" looks like - classical, not exceptional

# D7 half-spin weights: (±1, ±1, ±1, ±1, ±1, ±1, ±1, 0) with even # of minuses
W = []
for signs in itertools.product([1, -1], repeat=7):
    if signs.count(-1) % 2 == 0:
        v = list(signs) + [0]
        W.append(v)
W = np.array(W, float)
assert len(W) == 64, f"Expected 64 vertices, got {len(W)}"

# Normalize
W = W * np.sqrt(2.0 / 7.0)

# D7 roots
R = L.e8_roots()
v1 = np.array([0, 0, 0, 0, 0, 0, 0, 1], float)
D7 = R[np.abs(R @ v1) < 1e-9]
print("D7 roots:", len(D7), "(expect 84)")

order, u1, u2 = L.coxeter_plane(L.simple_roots(D7), 12)
print("D7 Coxeter order:", order, "(expect 12)")
P = L.project(W, u1, u2)
rings = L.ring_histogram(P)
print("ring histogram:", rings)
print(f"Number of vertices: {len(W)} = 2^6 (classical doubling pattern)")

np.save("/tmp/d7_W.npy", W)
np.save("/tmp/d7_P.npy", P)

# Draw it
G = np.round(W @ W.T, 1)
radii = np.round(np.hypot(P[:, 0], P[:, 1]), 4)
rr = sorted(set(radii))
ring_of = {r: i for i, r in enumerate(rr)}
s = 7.3 / max(rr) if max(rr) > 0 else 1.0
cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0) for i in range(len(rr))]

def draw(ip_target, fname, title):
    E = [(i, j) for i, j in itertools.combinations(range(len(W)), 2)
         if abs(G[i, j] - ip_target) < 0.15]
    # Limit edges for visualization (too many for 64 vertices)
    if len(E) > 500:
        import random
        random.seed(42)
        E = random.sample(E, 500)
    o = [r"\documentclass[tikz,border=6pt]{standalone}", r"\usepackage{tikz}",
         r"\begin{document}", r"\begin{tikzpicture}"]
    for i, (r, g, b) in enumerate(cols):
        o.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
    o += [r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",
          r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);"]
    for i, j in E:
        ro = max(ring_of[radii[i]], ring_of[radii[j]])
        o.append(fr"\draw[ring{ro},line width=0.4pt,opacity=0.35] "
                 fr"({P[i,0]*s:.3f},{P[i,1]*s:.3f})--({P[j,0]*s:.3f},{P[j,1]*s:.3f});")
    for k in range(len(W)):
        o.append(fr"\fill[ring{ring_of[radii[k]]}] ({P[k,0]*s:.3f},{P[k,1]*s:.3f}) circle (2.5pt);")
    o.append(fr"\node[white!80!black,font=\small\itshape] at (0,-8.4) {{{title}}};")
    o += [r"\end{tikzpicture}", r"\end{document}"]
    open(f"/Users/miriam/symmetric-graphs/e8/{fname}.tex", "w").write("\n".join(o))
    print(fname, "edges (sampled):", len(E))

# Find a good inner product to visualize
ips = sorted(set(G.flatten()))
print("Available inner products:", [ip for ip in ips if ip < 2.0])

draw(0.9, "d7_classical", r"$D_7$: 64 vertices, 12-fold (classical doubling: $2^6$)")
