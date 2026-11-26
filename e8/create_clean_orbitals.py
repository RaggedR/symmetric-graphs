"""
Create clean versions of orbital graphs without captions for use in Carter class images.
"""
import numpy as np, colorsys, itertools
from pathlib import Path

# Load E8 root data
P = np.load('/tmp/e8_P.npy')  # Coxeter plane projections
R = np.load('/tmp/e8_R.npy')  # 240 roots in R^8
n = len(R)

# Compute adjacency by inner product
G = np.round(R @ R.T).astype(int)

# Get radii and assign to rings
radii = np.round(np.hypot(P[:,0], P[:,1]), 4)
rings = sorted(set(radii))
ring_of = {r: i for i, r in enumerate(rings)}

# Color scheme
cols = [colorsys.hsv_to_rgb(0.02 + 0.80*i/7, 0.75, 1.0) for i in range(8)]
s = 3.5  # scale

def create_clean_orbital(ip, filename, opacity=0.08, linewidth=0.22):
    """Create a clean orbital graph without captions."""

    edges = [(i, j) for i, j in itertools.combinations(range(n), 2) if G[i, j] == ip]

    out = [
        r"\documentclass[tikz,border=6pt]{standalone}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}"
    ]

    # Define colors
    for i, (r, g, b) in enumerate(cols):
        out.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")

    # Background
    out.append(r"\definecolor{bg}{rgb}{0.04,0.05,0.09}")
    out.append(r"\fill[bg] (-5.5,-5.5) rectangle (5.5,5.5);")

    # Draw edges
    for i, j in edges:
        ro = max(ring_of[radii[i]], ring_of[radii[j]])
        x1, y1 = P[i, 0]*s, P[i, 1]*s
        x2, y2 = P[j, 0]*s, P[j, 1]*s
        out.append(fr"\draw[ring{ro},line width={linewidth}pt,opacity={opacity}] ({x1:.3f},{y1:.3f})--({x2:.3f},{y2:.3f});")

    # Draw vertices
    for k in range(n):
        x, y = P[k, 0]*s, P[k, 1]*s
        out.append(fr"\fill[ring{ring_of[radii[k]]}] ({x:.3f},{y:.3f}) circle (0.7pt);")

    out += [r"\end{tikzpicture}", r"\end{document}"]

    # Write
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / f"{filename}.tex"
    with open(tex_path, 'w') as f:
        f.write("\n".join(out))

    return tex_path

# Create clean versions of the main orbitals
print("Creating clean orbital graphs (no captions)...\n")

orbitals = [
    (+2, "clean_e8_id", 0.08, 0.22),         # identity (no edges actually)
    (+1, "clean_e8_mandala", 0.08, 0.22),    # mandala
    (0, "clean_e8_bloom", 0.09, 0.07),       # bloom
    (-1, "clean_e8_weave", 0.08, 0.22),      # weave
    (-2, "clean_e8_antipodal", 0.60, 0.22),  # antipodal
]

created = []
for ip, name, opacity, linewidth in orbitals:
    print(f"Creating {name} (ip={ip:+d})...", end=" ")
    tex_path = create_clean_orbital(ip, name, opacity, linewidth)
    created.append((name, tex_path))
    print("✓")

print(f"\n✓ Created {len(created)} clean orbital TEX files")
print("Ready to compile with pdflatex")

# Save list
with open("/Users/miriam/symmetric-graphs/e8/clean_orbitals_list.txt", "w") as f:
    for name, tex in created:
        f.write(f"{name}\t{tex}\n")
