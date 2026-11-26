"""
Compute all Hecke algebra orbitals needed for the 74 Carter classes.

For each root system (G2, F4, A2, A3, D4, D6, etc.), we:
1. Load the root system
2. Compute all orbital graphs by inner products
3. Generate clean PDFs for each orbital

This gives us enough unique graphs for a 1-1 mapping.
"""

import numpy as np
import colorsys
import itertools
import subprocess
from pathlib import Path

def create_orbital_graph(roots_P, roots_R, ip_value, output_name, opacity=0.08, linewidth=0.22):
    """Create a single orbital graph from root projections."""

    n = len(roots_R)
    G = np.round(roots_R @ roots_R.T).astype(int)

    # Get radii for coloring
    radii = np.round(np.hypot(roots_P[:,0], roots_P[:,1]), 4)
    rings = sorted(set(radii))
    ring_of = {r: i for i, r in enumerate(rings)}

    # Scale
    s = 3.5

    # Colors
    cols = [colorsys.hsv_to_rgb(0.02 + 0.80*i/max(1, len(rings)-1), 0.75, 1.0)
            for i in range(len(rings))]

    # Edges for this orbital
    edges = [(i, j) for i, j in itertools.combinations(range(n), 2)
             if G[i, j] == ip_value]

    # Generate TikZ
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
        x1, y1 = roots_P[i, 0]*s, roots_P[i, 1]*s
        x2, y2 = roots_P[j, 0]*s, roots_P[j, 1]*s
        out.append(fr"\draw[ring{ro},line width={linewidth}pt,opacity={opacity}] ({x1:.3f},{y1:.3f})--({x2:.3f},{y2:.3f});")

    # Draw vertices
    for k in range(n):
        x, y = roots_P[k, 0]*s, roots_P[k, 1]*s
        out.append(fr"\fill[ring{ring_of[radii[k]]}] ({x:.3f},{y:.3f}) circle (0.7pt);")

    out += [r"\end{tikzpicture}", r"\end{document}"]

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / f"{output_name}.tex"
    with open(tex_path, 'w') as f:
        f.write("\n".join(out))

    return tex_path, len(edges)


def compile_pdf(tex_path):
    """Compile a tex file to PDF."""
    result = subprocess.run(
        ['pdflatex', '-interaction=nonstopmode', str(tex_path)],
        cwd='/Users/miriam/symmetric-graphs/e8',
        capture_output=True,
        timeout=30
    )
    return result.returncode == 0


# Check which root systems we need orbitals for
print("Computing Hecke algebra orbitals...\n")

# We already have E8/E7 orbitals (5 of them)
# We need orbitals for the root systems that have multiple classes

root_systems_needed = {
    "G2": 10,  # 10 G2 classes
    "F4": 7,   # 7 F4 classes
    "D4": 7,   # 7 D4 classes
    "D6": 3,   # 3 D6 classes
    "D7": 3,   # 3 D7/D8 classes
    "A3": 10,  # 10 A3-dominant classes
    "A2": 9,   # 9 A2-dominant classes
    "E5/D5": 4, # 4 D5 classes
    "E6": 4,   # 4 E6-related classes
    "E7": 2,   # 2 E7-related classes
}

print("Root systems needing multiple orbitals:")
for rs, count in root_systems_needed.items():
    print(f"  {rs}: need {count} orbitals")

print("\nStrategy:")
print("1. For each root system, compute inner product matrix")
print("2. Extract all unique ip values (these are the orbitals)")
print("3. Generate a clean graph for each orbital")
print("4. This gives us enough unique graphs for 1-1 mapping")

print("\nNOTE: We need to actually run the build scripts (d4_build.py, etc.)")
print("      to get the root projections, then extract multiple orbitals.")
print("\nFor now, let's check if we have the data files:")

import os
data_files = [
    "/tmp/e8_P.npy",
    "/tmp/e8_R.npy",
    "/tmp/e7_W.npy",
    "/tmp/e7_P.npy",
]

for f in data_files:
    exists = "✓" if os.path.exists(f) else "✗"
    print(f"  {exists} {f}")

print("\n" + "="*60)
print("NEXT STEPS:")
print("="*60)
print("1. Run each build script to generate root projections:")
print("   python3 e8/d4_build.py  → /tmp/d4_*.npy")
print("   python3 e8/d6_build.py  → /tmp/d6_*.npy")
print("   python3 e8/exceptional_build.py → F4, G2, A2, A3")
print("")
print("2. For each root system, extract ALL orbitals (all ip values)")
print("   instead of just one")
print("")
print("3. Generate ~50 new orbital PDFs to get 74 unique total")
print("")
print("4. Update mapping to use unique orbital for each class")
