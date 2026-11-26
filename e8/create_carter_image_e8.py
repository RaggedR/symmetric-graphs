"""
Create E8 Carter class image: orbital graph (mandala) + Dynkin diagram in top left.
"""
import numpy as np, colorsys, itertools

# Load the precomputed E8 root projections
P = np.load('/tmp/e8_P.npy')  # 240 x 2 (Coxeter plane projection)
R = np.load('/tmp/e8_R.npy')  # 240 x 8 (roots in R^8)
n = len(R)

# Compute adjacency by inner product
G = np.round(R @ R.T).astype(int)

# Get radii and assign to rings
radii = np.round(np.hypot(P[:,0], P[:,1]), 4)
rings = sorted(set(radii))
ring_of = {r: i for i, r in enumerate(rings)}

# Color scheme: 8 rings with hue gradient
cols = [colorsys.hsv_to_rgb(0.02 + 0.80*i/7, 0.75, 1.0) for i in range(8)]

# Which edges to draw: ip = +1 (the mandala, 56-regular)
edges = [(i, j) for i, j in itertools.combinations(range(n), 2) if G[i, j] == 1]

# Scaling to fit nicely
s = 3.5  # scale factor for good visual size

# --- Generate TikZ ---
out = [
    r"\documentclass[tikz,border=6pt]{standalone}",
    r"\usepackage{tikz}",
    r"\usetikzlibrary{positioning}",
    r"\begin{document}",
    r"\begin{tikzpicture}"
]

# Define ring colors
for i, (r, g, b) in enumerate(cols):
    out.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")

# Background
out.append(r"\definecolor{bg}{rgb}{0.04,0.05,0.09}")
out.append(r"\fill[bg] (-5.5,-5.5) rectangle (5.5,5.5);")

# Draw edges (the mandala: ip=+1)
for i, j in edges:
    ro = max(ring_of[radii[i]], ring_of[radii[j]])
    x1, y1 = P[i, 0]*s, P[i, 1]*s
    x2, y2 = P[j, 0]*s, P[j, 1]*s
    out.append(fr"\draw[ring{ro},line width=0.22pt,opacity=0.08] ({x1:.3f},{y1:.3f})--({x2:.3f},{y2:.3f});")

# Draw vertices
for k in range(n):
    x, y = P[k, 0]*s, P[k, 1]*s
    out.append(fr"\fill[ring{ring_of[radii[k]]}] ({x:.3f},{y:.3f}) circle (0.7pt);")

# --- Add E8 Dynkin diagram in top left corner ---
# E8 Dynkin diagram:
#       o---o---o---o---o---o---o
#                   |
#                   o
# Nodes 1-7 in a line, node 8 attached to node 3

out.append(r"% E8 Dynkin diagram")
out.append(r"\begin{scope}[shift={(-4.5,4.5)},scale=0.25]")

# Node positions (compact layout)
node_positions = [
    (0, 0),   # 1
    (1, 0),   # 2
    (2, 0),   # 3
    (3, 0),   # 4
    (4, 0),   # 5
    (5, 0),   # 6
    (6, 0),   # 7
    (2, -1),  # 8 (hanging below node 3)
]

# Draw edges
edges_dynkin = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (2,7)]
for i, j in edges_dynkin:
    x1, y1 = node_positions[i]
    x2, y2 = node_positions[j]
    out.append(fr"\draw[white!60!bg,line width=0.6pt] ({x1},{y1})--({x2},{y2});")

# Draw nodes
for x, y in node_positions:
    out.append(fr"\fill[white!90!bg] ({x},{y}) circle (0.15);")

out.append(r"\end{scope}")

out += [r"\end{tikzpicture}", r"\end{document}"]

# Write to file
tex_path = '/Users/miriam/symmetric-graphs/e8/carter_image_e8.tex'
with open(tex_path, 'w') as f:
    f.write("\n".join(out))

print(f"✓ Created {tex_path}")
print(f"  Edges drawn: {len(edges)} (ip=+1, the mandala)")
