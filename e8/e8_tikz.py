import numpy as np, colorsys, itertools
P = np.load('/tmp/e8_P.npy'); R = np.load('/tmp/e8_R.npy')
n = len(R)
radii = np.round(np.hypot(P[:,0], P[:,1]), 4)
rings = sorted(set(radii))
ring_of = {r:i for i,r in enumerate(rings)}
scale = 7.3/max(rings)                     # outer ring -> 7.3 cm

# 8 rainbow colours (inner=red .. outer=violet)
cols = []
for i in range(8):
    h = 0.02 + 0.80*i/7
    r,g,b = colorsys.hsv_to_rgb(h, 0.75, 1.0)
    cols.append((r,g,b))

# edges: inner-product == 1
G = R @ R.T
edges = [(i,j) for i,j in itertools.combinations(range(n),2) if abs(G[i,j]-1.0)<1e-6]

def xy(k): return P[k,0]*scale, P[k,1]*scale

out = []
out.append(r"\documentclass[tikz,border=6pt]{standalone}")
out.append(r"\usepackage{tikz}")
out.append(r"\begin{document}")
out.append(r"\begin{tikzpicture}")
for i,(r,g,b) in enumerate(cols):
    out.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
out.append(r"\definecolor{bg}{rgb}{0.04,0.05,0.09}")
# background
out.append(r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);")
# edges first (coloured by the OUTER endpoint's ring, low opacity)
for i,j in edges:
    ro = max(ring_of[radii[i]], ring_of[radii[j]])
    x1,y1 = xy(i); x2,y2 = xy(j)
    out.append(fr"\draw[ring{ro},line width=0.09pt,opacity=0.22] ({x1:.4f},{y1:.4f})--({x2:.4f},{y2:.4f});")
# vertices on top
for k in range(n):
    x,y = xy(k); c = ring_of[radii[k]]
    out.append(fr"\fill[ring{c}] ({x:.4f},{y:.4f}) circle (1.5pt);")
out.append(r"\node[white!80!black,font=\small\itshape] at (0,-8.4)"
           r" {$E_8$: 240 roots on the Coxeter plane, 30-fold symmetry};")
out.append(r"\end{tikzpicture}")
out.append(r"\end{document}")
open('/Users/robin/git/graphviz/e8/e8_mandala.tex','w').write("\n".join(out))
print("edges:", len(edges), "(expect 6720)")
print("wrote e8_mandala.tex")
