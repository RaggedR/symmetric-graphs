"""
Extract ALL orbitals from each root system to get 74+ unique graphs.

For each root system:
1. Load P (Coxeter projections) and R (roots)
2. Compute inner products G = R @ R.T
3. Find all unique ip values
4. Generate a clean orbital graph for each unique ip
"""

import numpy as np
import colorsys
import itertools
import subprocess
from pathlib import Path
from collections import Counter

def extract_orbitals_from_roots(roots_file_P, roots_file_R, system_name):
    """Extract all orbitals from a root system."""

    P = np.load(roots_file_P)
    R = np.load(roots_file_R)
    n = len(R)

    # Compute inner products
    G = np.round(R @ R.T).astype(int)

    # Find all unique ip values (excluding diagonal which is always 2)
    ip_values = []
    for i in range(n):
        for j in range(i+1, n):
            ip_values.append(G[i, j])

    ip_counts = Counter(ip_values)
    unique_ips = sorted(ip_counts.keys())

    print(f"\n{system_name}:")
    print(f"  {n} roots")
    print(f"  {len(unique_ips)} orbitals (unique ip values): {unique_ips}")
    for ip in unique_ips:
        count = ip_counts[ip]
        print(f"    ip={ip:+2d}: {count} pairs")

    return P, R, unique_ips


def create_orbital_pdf(P, R, ip_value, output_name):
    """Create a clean orbital graph PDF."""

    n = len(R)
    G = np.round(R @ R.T).astype(int)

    # Get radii for coloring
    radii = np.round(np.hypot(P[:,0], P[:,1]), 4)
    rings = sorted(set(radii))
    ring_of = {r: i for i, r in enumerate(rings)}

    # Scale
    s = 3.5

    # Colors
    n_rings = max(len(rings), 1)
    cols = [colorsys.hsv_to_rgb(0.02 + 0.80*i/max(1, n_rings-1), 0.75, 1.0)
            for i in range(n_rings)]

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
    opacity = 0.08 if len(edges) > 100 else 0.3
    linewidth = 0.22 if len(edges) > 100 else 0.4

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

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / f"{output_name}.tex"
    with open(tex_path, 'w') as f:
        f.write("\n".join(out))

    # Compile
    result = subprocess.run(
        ['pdflatex', '-interaction=nonstopmode', str(tex_path)],
        cwd='/Users/miriam/symmetric-graphs/e8',
        capture_output=True,
        timeout=30
    )

    if result.returncode == 0:
        return f"{output_name}.pdf", len(edges)
    else:
        return None, 0


print("="*70)
print("EXTRACTING ALL ORBITALS FROM ROOT SYSTEMS")
print("="*70)

# Track all generated orbitals
all_orbitals = []

# Root systems to process (file paths in /tmp)
root_systems = [
    ("/tmp/e8_P.npy", "/tmp/e8_R.npy", "E8", "e8"),
    ("/tmp/e7_P.npy", "/tmp/e7_W.npy", "E7", "e7"),  # W for weights
    ("/tmp/e6_P.npy", "/tmp/e6_W.npy", "E6", "e6"),
    ("/tmp/g2_P.npy", "/tmp/g2_R.npy", "G2", "g2"),
    ("/tmp/f4_P.npy", "/tmp/f4_R.npy", "F4", "f4"),
    ("/tmp/a3_P.npy", "/tmp/a3_R.npy", "A3", "a3"),
    ("/tmp/a2_P.npy", "/tmp/a2_R.npy", "A2", "a2"),
    ("/tmp/d7_P.npy", "/tmp/d7_W.npy", "D7", "d7"),
    ("/tmp/d6_P.npy", "/tmp/d6_W.npy", "D6", "d6"),
    ("/tmp/d4_P.npy", "/tmp/d4_W.npy", "D4", "d4"),
]

# Check which files exist
import os
for P_file, R_file, name, prefix in root_systems[:]:
    if not (os.path.exists(P_file) and os.path.exists(R_file)):
        print(f"⚠ Skipping {name}: files not found")
        root_systems.remove((P_file, R_file, name, prefix))

print(f"\nProcessing {len(root_systems)} root systems...")

for P_file, R_file, name, prefix in root_systems:
    P, R, unique_ips = extract_orbitals_from_roots(P_file, R_file, name)

    # Generate PDFs for each orbital (limit to avoid too many)
    max_orbitals = 10
    for idx, ip in enumerate(unique_ips[:max_orbitals]):
        output_name = f"orbital_{prefix}_ip{ip:+d}".replace("+", "p").replace("-", "m")
        print(f"  Generating {output_name}...", end=" ")

        pdf_name, n_edges = create_orbital_pdf(P, R, ip, output_name)

        if pdf_name:
            print(f"✓ ({n_edges} edges)")
            all_orbitals.append((name, ip, pdf_name, n_edges))
        else:
            print("✗")

print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")
print(f"Generated {len(all_orbitals)} new orbital PDFs")
print(f"Combined with existing clean_*.pdf files, we should have 70+ unique graphs")
print(f"\nSaving orbital list...")

# Save list
with open("/Users/miriam/symmetric-graphs/e8/all_orbitals_generated.txt", "w") as f:
    f.write("# All generated Hecke algebra orbitals\n")
    f.write("# Format: system | ip_value | pdf_file | num_edges\n\n")
    for system, ip, pdf, edges in all_orbitals:
        f.write(f"{system:4s} | {ip:+2d} | {pdf:30s} | {edges:5d}\n")

print("✓ Saved to all_orbitals_generated.txt")
print(f"\nNext: Update mapping to use these {len(all_orbitals)} new orbitals!")
