"""
Build the D4 rung: W(D4)/W(A3), the 8 half-spin minuscule weights (vertices of
the 16-cell / demihypercube in 4D).  Descent: fix a weight w2 in the 16 (D5
rung); the 8 are the D5-weights at inner product +1 to w2.  The scheme is
rank 3 (strongly regular).  D4 Coxeter number h = 6, giving rings on the
Coxeter plane.

Writes /tmp/d4_W.npy, /tmp/d4_P.npy and drawings.
"""
import numpy as np, colorsys, itertools
import ladder as L

# Build D4 weights directly: 8 vertices of the 16-cell (demihypercube in 4D)
# These are vectors (±1, ±1, ±1, ±1, 0, 0, 0, 0) with an even number of minuses
W = []
for signs in itertools.product([1, -1], repeat=4):
    if signs.count(-1) % 2 == 0:
        v = list(signs) + [0, 0, 0, 0]
        W.append(v)
W = np.array(W, float)
assert len(W) == 8, f"Expected 8 vertices, got {len(W)}"

# Normalize to match E8 conventions (length^2 = 2)
# Current vectors have length^2 = 4, so scale by sqrt(2/4) = 1/sqrt(2)
W = W / np.sqrt(2.0)

# D4 root subsystem
R = L.e8_roots()
# Fix 4 directions to get D4 (remaining 4D)
v1 = np.array([0, 0, 0, 0, 1, 1, 0, 0], float)
v2 = np.array([0, 0, 0, 0, 1, -1, 0, 0], float)
v3 = np.array([0, 0, 0, 0, 0, 0, 1, 1], float)
v4 = np.array([0, 0, 0, 0, 0, 0, 1, -1], float)
D4 = R[(np.abs(R @ v1) < 1e-9) & (np.abs(R @ v2) < 1e-9) &
       (np.abs(R @ v3) < 1e-9) & (np.abs(R @ v4) < 1e-9)]
print("D4 roots:", len(D4), "(expect 24)")
order, u1, u2 = L.coxeter_plane(L.simple_roots(D4), 6)
print("D4 Coxeter order:", order, "(expect 6)")
P = L.project(W, u1, u2)
rings = L.ring_histogram(P)
print("ring histogram:", rings)

np.save("/tmp/d4_W.npy", W)
np.save("/tmp/d4_P.npy", P)

# --- draw the orbital graphs ---
G = np.round(W @ W.T).astype(int)
radii = np.round(np.hypot(P[:, 0], P[:, 1]), 4)
rr = sorted(set(radii))
ring_of = {r: i for i, r in enumerate(rr)}
s = 7.3 / max(rr)
cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0) for i in range(len(rr))]


def draw(ip, fname, title):
    E = [(i, j) for i, j in itertools.combinations(range(len(W)), 2) if G[i, j] == ip]
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


draw(0, "d4_graph", r"$D_4$: 8 vertices, 6-fold symmetry")
draw(1, "d4_complement", r"$D_4$: complement graph")
