"""
Build A1, A4, and A5 root systems and extract all orbitals.
Following the pattern from exceptional_build.py
"""
import numpy as np
import itertools
import sys
sys.path.insert(0, '/Users/miriam/symmetric-graphs/e8')
import ladder as L

def build_a_n(n):
    """
    Build A_n root system in R^(n+1).

    A_n has (n+1)*n roots: e_i - e_j for all i ≠ j
    where e_i are standard basis vectors in R^(n+1)

    For embedding in R^8, we use the first (n+1) coordinates.
    """
    roots = []

    # e_i - e_j for all i ≠ j (i, j = 0, ..., n)
    for i, j in itertools.permutations(range(n+1), 2):
        v = [0] * 8
        v[i] = 1
        v[j] = -1
        roots.append(v)

    roots = np.array(roots, float)
    print(f"A{n}: built {len(roots)} roots (expect {(n+1)*n})")

    # Verify orthogonality properties
    G = np.round(roots @ roots.T).astype(int)
    unique_ips = sorted(set(G.flatten()))
    print(f"  Inner product values: {unique_ips}")

    return roots


def extract_orbitals(roots, name, h):
    """Extract all orbitals from a root system."""
    # Get Coxeter plane projection
    simple = L.simple_roots(roots)
    order, u1, u2 = L.coxeter_plane(simple, h)
    print(f"{name}: Coxeter order {order} (expect {h})")

    P = L.project(roots, u1, u2)

    # Save projections
    np.save(f"/tmp/{name}_R.npy", roots)
    np.save(f"/tmp/{name}_P.npy", P)

    # Compute inner products
    G = np.round(roots @ roots.T).astype(int)
    unique_ips = sorted(set(G.flatten()))

    # Remove diagonal (ip = 2) and duplicates
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
            create_orbital_tikz(P, edges, name, ip)

    return orbitals


def create_orbital_tikz(P, edges, root_system, ip_value):
    """Create a clean TikZ drawing of an orbital graph."""
    import colorsys

    # Determine ip label
    if ip_value >= 0:
        ip_label = f"ipp{ip_value}"
    else:
        ip_label = f"ipm{abs(ip_value)}"

    filename = f"orbital_{root_system.lower()}_{ip_label}"

    # Handle NaN values in projection (happens for low-dimensional systems)
    # Replace NaN with small random offsets
    if np.any(np.isnan(P)):
        print(f"  WARNING: NaN in projection for {root_system}, using fallback layout")
        n = len(P)
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        P = np.column_stack([np.cos(angles), np.sin(angles)])

    # Scale to fit
    radii = np.hypot(P[:, 0], P[:, 1])
    rr = sorted(set(np.round(radii, 4)))

    # Handle case where all radii are the same (or very close)
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

    # No caption for clean orbitals
    o += [r"\end{tikzpicture}", r"\end{document}"]

    # Write file
    tex_path = f"/Users/miriam/symmetric-graphs/e8/{filename}.tex"
    with open(tex_path, 'w') as f:
        f.write("\n".join(o))

    print(f"  Created {filename}.tex ({len(edges)} edges)")

    return filename


def main():
    """Build A1, A4, A5 and extract all orbitals."""
    print("="*70)
    print("Building A-series root systems (A1, A4, A5)")
    print("="*70)
    print()

    systems = [
        (1, 2),   # A1: 2 roots, Coxeter h=2
        (4, 5),   # A4: 20 roots, Coxeter h=5
        (5, 6),   # A5: 30 roots, Coxeter h=6
    ]

    total_orbitals = 0

    for n, h in systems:
        print(f"\n{'='*70}")
        print(f"A{n} Root System")
        print(f"{'='*70}\n")

        # Build roots
        roots = build_a_n(n)

        # Extract orbitals
        orbitals = extract_orbitals(roots, f"A{n}", h)

        print(f"\nA{n} Summary:")
        print(f"  Roots: {len(roots)}")
        print(f"  Orbitals: {len(orbitals)}")
        for orb in orbitals:
            print(f"    ip={orb['ip']:+2d}: {orb['num_edges']:5d} edges")

        total_orbitals += len(orbitals)

    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Generated {total_orbitals} new orbitals from A-series")
    print("Files created:")
    print("  - orbital_a1_*.tex/pdf")
    print("  - orbital_a4_*.tex/pdf")
    print("  - orbital_a5_*.tex/pdf")
    print()
    print("Next: Compile with pdflatex to generate PDFs")


if __name__ == '__main__':
    main()
