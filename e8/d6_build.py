"""
Build D6 to show the classical D-series pattern continuing beyond D5.
D6 has 32 half-spin minuscule weights, Coxeter number h = 10.

This is NOT part of the exceptional ladder (there's no E5.5 containing D6),
but serves to show what the classical pattern looks like when extended
parallel to E6/E7/E8.

The 32 vertices are constructed from E7 weights that are orthogonal to a
suitable vector, or directly as vectors in R^8.
"""
import numpy as np, colorsys, itertools
import ladder as L

# Build D6 weights directly: take D5 extended with its "sibling" spinor rep
# D6 has 32 vertices (half-spin weights)
# Construct from E8 roots orthogonal to two vectors (forms D6 subsystem)

R = L.e8_roots()

# Find a suitable pair of orthogonal roots to fix D6
# We want D6, so we need to fix 2 directions in E8 (leaving 6D)
v1 = np.array([1, 0, 0, 0, 0, 0, 0, 1], float)  # a suitable root
v2 = np.array([0, 1, 0, 0, 0, 0, 1, 0], float)  # orthogonal to v1

# D6 roots = E8 roots orthogonal to v1 and v2
D6_roots = R[(np.abs(R @ v1) < 1e-9) & (np.abs(R @ v2) < 1e-9)]
print("D6 roots:", len(D6_roots), "(expect 60)")

# Build the half-spin weights
# For D6, we can construct the 32 half-spin weights from the root system
# They form the vertices of the 6-demihypercube

# Method: take all vectors (±1, ±1, ±1, ±1, ±1, ±1, 0, 0) with even # of minuses
W = []
for signs in itertools.product([1, -1], repeat=6):
    if signs.count(-1) % 2 == 0:
        v = list(signs) + [0, 0]
        W.append(v)
W = np.array(W, float)
assert len(W) == 32, f"Expected 32 vertices, got {len(W)}"

# Normalize to match E8 conventions (length^2 = 2)
# Current vectors have length^2 = 6, so scale by sqrt(2/6)
W = W * np.sqrt(2.0 / 6.0)

# Get Coxeter plane
order, u1, u2 = L.coxeter_plane(L.simple_roots(D6_roots), 10)
print("D6 Coxeter order:", order, "(expect 10)")
P = L.project(W, u1, u2)
rings = L.ring_histogram(P)
print("ring histogram:", rings)

np.save("/tmp/d6_W.npy", W)
np.save("/tmp/d6_P.npy", P)

# --- draw the orbital graphs ---
G = np.round(W @ W.T, 1)  # inner products
radii = np.round(np.hypot(P[:, 0], P[:, 1]), 4)
rr = sorted(set(radii))
ring_of = {r: i for i, r in enumerate(rr)}
s = 7.3 / max(rr) if max(rr) > 0 else 1.0
cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0) for i in range(len(rr))]


def draw(ip_target, fname, title):
    E = [(i, j) for i, j in itertools.combinations(range(len(W)), 2)
         if abs(G[i, j] - ip_target) < 0.15]
    o = [r"\documentclass[tikz,border=6pt]{standalone}", r"\usepackage{tikz}",
         r"\begin{document}", r"\begin{tikzpicture}"]
    for i, (r, g, b) in enumerate(cols):
        o.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
    o += [r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",
          r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);"]
    for i, j in E:
        ro = max(ring_of[radii[i]], ring_of[radii[j]])
        o.append(fr"\draw[ring{ro},line width=0.5pt,opacity=0.55] "
                 fr"({P[i,0]*s:.3f},{P[i,1]*s:.3f})--({P[j,0]*s:.3f},{P[j,1]*s:.3f});")
    for k in range(len(W)):
        o.append(fr"\fill[ring{ring_of[radii[k]]}] ({P[k,0]*s:.3f},{P[k,1]*s:.3f}) circle (3.0pt);")
    o.append(fr"\node[white!80!black,font=\small\itshape] at (0,-8.4) {{{title}}};")
    o += [r"\end{tikzpicture}", r"\end{document}"]
    open(f"/Users/miriam/symmetric-graphs/e8/{fname}.tex", "w").write("\n".join(o))
    print(fname, "edges:", len(E))


# Draw main adjacency graph at ip ≈ 0.67 (the 240-edge graph)
draw(0.67, "d6_graph", r"$D_6$: 32 vertices, 10-fold symmetry (classical series)")
