"""
Build the D5 (= E5) rung: W(D5)/W(A4), the 16 half-spin minuscule weights.
Descent: fix a weight w1 in the 27 (E6 rung); the 16 are the E6-weights at
inner product +1 to w1.  The scheme is rank 3 (strongly regular): the Clebsch
graph SRG(16,5,0,2) (ip 0) and its complement, the halved 5-cube SRG(16,10,6,6)
(ip +1).  D5 Coxeter number h = 8, giving 2 rings of 8 on the Coxeter plane.

Writes /tmp/e5_W.npy, /tmp/e5_P.npy and the two drawings e5_clebsch.tex,
e5_halved5cube.tex.
"""
import numpy as np, colorsys, itertools
import ladder as L

W6 = np.load("/tmp/e6_W.npy")
w1 = W6[0]
W = W6[np.round(W6 @ w1).astype(int) == 1]          # the 16
assert len(W) == 16, len(W)

# D5 root subsystem = E8 roots orthogonal to r0, w0, w1
R = L.e8_roots()
D5 = R[(np.abs(R @ L.R0) < 1e-9) & (np.abs(R @ L.W0) < 1e-9) & (np.abs(R @ w1) < 1e-9)]
print("D5 roots:", len(D5), "(expect 40)")
order, u1, u2 = L.coxeter_plane(L.simple_roots(D5), 8)
print("D5 Coxeter order:", order, "(expect 8)")
P = L.project(W, u1, u2)
rings = L.ring_histogram(P)
print("ring histogram:", rings, "(expect 2 rings of 8)")

np.save("/tmp/e5_W.npy", W)
np.save("/tmp/e5_P.npy", P)

# --- draw the two orbital graphs, same style as the other rungs ---
G = np.round(W @ W.T).astype(int)
radii = np.round(np.hypot(P[:, 0], P[:, 1]), 4)
rr = sorted(set(radii))
ring_of = {r: i for i, r in enumerate(rr)}
s = 7.3 / max(rr)
cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0) for i in range(len(rr))]


def draw(ip, fname, title):
    E = [(i, j) for i, j in itertools.combinations(range(16), 2) if G[i, j] == ip]
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
    for k in range(16):
        o.append(fr"\fill[ring{ring_of[radii[k]]}] ({P[k,0]*s:.3f},{P[k,1]*s:.3f}) circle (3.0pt);")
    o.append(fr"\node[white!80!black,font=\small\itshape] at (0,-8.4) {{{title}}};")
    o += [r"\end{tikzpicture}", r"\end{document}"]
    open(f"/Users/miriam/symmetric-graphs/e8/{fname}.tex", "w").write("\n".join(o))
    print(fname, "edges:", len(E))


draw(0, "e5_clebsch", r"$D_5$: the Clebsch graph SRG(16,5,0,2)")
draw(1, "e5_halved5cube", r"$D_5$: complement = halved 5-cube SRG(16,10,6,6)")
