"""
Extract D5 orbitals from existing D5 root system data.
D5 projections should already exist from e5_build.py
"""
import numpy as np
import sys
import colorsys
sys.path.insert(0, '/Users/miriam/symmetric-graphs/e8')

def extract_d5_orbitals():
    """Extract D5 orbitals from cached data."""

    # Check if D5 data exists
    try:
        # Try loading from e5_build.py cache (D5 is same as E5 in terms of roots)
        # D5 has 40 roots: (±e_i ± e_j) for i < j in {1,2,3,4,5}

        # Build D5 roots directly
        roots = []
        for i in range(5):
            for j in range(i+1, 5):
                for si, sj in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                    v = [0] * 8
                    v[i] = si
                    v[j] = sj
                    roots.append(v)

        roots = np.array(roots, float)
        print(f"D5: built {len(roots)} roots (expect 40)")

        # Simple projection: use first 2 coordinates
        # For better visualization, project onto a 2D plane
        # Use simple projection for now
        P = np.column_stack([roots[:, 0], roots[:, 1]])

        # For better Coxeter plane, we'd need to compute it properly
        # But for simplicity, let's use PCA
        from numpy.linalg import eig
        cov = roots[:, :5].T @ roots[:, :5]
        eigenvalues, eigenvectors = eig(cov)

        # Sort by eigenvalue
        idx = eigenvalues.argsort()[::-1]
        eigenvectors = eigenvectors[:, idx]

        # Project onto first 2 principal components
        P = roots[:, :5] @ eigenvectors[:, :2]

        print(f"  Projected to 2D")

    except Exception as e:
        print(f"Error loading D5 data: {e}")
        return []

    # Compute inner products
    G = np.round(roots @ roots.T).astype(int)
    unique_ips = sorted(set(G.flatten()))

    # Remove diagonal (ip = 2)
    unique_ips = [ip for ip in unique_ips if ip != 2]

    print(f"  {len(unique_ips)} orbitals (ip values: {unique_ips})")

    orbitals = []
    for ip in unique_ips:
        # Find all pairs with this inner product
        edges = []
        for i in range(len(roots)):
            for j in range(i+1, len(roots)):
                if G[i, j] == ip:
                    edges.append((i, j))

        if edges:
            orbitals.append({
                'ip': ip,
                'edges': edges,
                'num_edges': len(edges)
            })

            # Create TikZ visualization
            create_orbital_tikz(P, edges, 'D5', ip)

    return orbitals


def create_orbital_tikz(P, edges, root_system, ip_value):
    """Create a clean TikZ drawing of an orbital graph."""

    # Determine ip label
    if ip_value >= 0:
        ip_label = f"ipp{ip_value}"
    else:
        ip_label = f"ipm{abs(ip_value)}"

    filename = f"orbital_{root_system.lower()}_{ip_label}"

    # Scale to fit
    P = np.asarray(P, dtype=float)  # Ensure float type
    radii = np.hypot(P[:, 0], P[:, 1])
    rr = sorted(set(np.round(radii, 4)))

    if len(rr) == 0:
        rr = [1.0]

    ring_of = {r: i for i, r in enumerate(rr)}
    s = 7.3 / max(rr) if max(rr) > 0 else 1.0

    # Color scheme: rings from yellow to orange
    cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0)
            for i in range(len(rr))]

    # Build TikZ
    o = [
        r"\documentclass[tikz,border=6pt]{standalone}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}"
    ]

    # Define colors
    for i, (r, g, b) in enumerate(cols):
        o.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")

    o += [
        r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",
        r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);"
    ]

    # Draw edges
    vertex_radii = np.round(radii, 4)
    for i, j in edges:
        ro = max(ring_of[vertex_radii[i]], ring_of[vertex_radii[j]])
        o.append(
            fr"\draw[ring{ro},line width=0.5pt,opacity=0.55] "
            fr"({P[i,0]*s:.3f},{P[i,1]*s:.3f})--({P[j,0]*s:.3f},{P[j,1]*s:.3f});"
        )

    # Draw vertices
    for k in range(len(P)):
        o.append(
            fr"\fill[ring{ring_of[vertex_radii[k]]}] "
            fr"({P[k,0]*s:.3f},{P[k,1]*s:.3f}) circle (3.0pt);"
        )

    o += [r"\end{tikzpicture}", r"\end{document}"]

    # Write file
    tex_path = f"/Users/miriam/symmetric-graphs/e8/{filename}.tex"
    with open(tex_path, 'w') as f:
        f.write("\n".join(o))

    print(f"  Created {filename}.tex ({len(edges)} edges)")

    return filename


def main():
    print("="*70)
    print("Extracting D5 Orbitals")
    print("="*70)
    print()

    orbitals = extract_d5_orbitals()

    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Generated {len(orbitals)} D5 orbitals")
    for orb in orbitals:
        print(f"  ip={orb['ip']:+2d}: {orb['num_edges']:5d} edges")
    print()


if __name__ == '__main__':
    main()
