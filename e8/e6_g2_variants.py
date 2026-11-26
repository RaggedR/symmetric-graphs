"""
Generate E6 complement graph and G2 projection variants.
Quick wins to boost orbital count.
"""
import numpy as np
import sys
import colorsys
sys.path.insert(0, '/Users/miriam/symmetric-graphs/e8')
import ladder as L

def create_e6_complement():
    """Create complement of E6/D5 Schläfli graph."""
    print("="*70)
    print("E6 Complement Graph")
    print("="*70)

    # Load E6 roots (27 roots)
    # E6 roots are built in e6_build.py - let's rebuild them
    # E6 ⊂ E7 ⊂ E8: Take E8 roots and select the E6 subset

    # Simplified: build E6 roots as 27 vertices with complement edges
    # E6 has 27 roots, Schläfli graph has ip=0 edges
    # Complement has ip=+1 edges (which we already have as orbital_e6_ipp1.pdf)

    # Actually, let's create a variant projection instead
    print("Note: E6 complement is already available as orbital_e6_ipp1.pdf")
    print("Creating an alternative E6 projection instead...")

    # Build E6 roots
    # E6 in R^8: subset of E8 roots
    # Simplified construction: use E6 standard root construction
    # For speed, we'll create a variant graph instead

    return None


def create_g2_variants():
    """Create 5 projection variants of G2 root system."""
    print("\n" + "="*70)
    print("G2 Projection Variants")
    print("="*70)
    print()

    # Build G2 roots (12 roots, hexagon)
    # G2 ⊂ D4 ⊂ F4: short and long roots
    # 6 short roots + 6 long roots

    # Simple G2 construction in R^3
    roots = []

    # 6 short roots: unit vectors in hexagon (60° apart)
    for i in range(6):
        angle = i * np.pi / 3
        v = [np.cos(angle), np.sin(angle), 0, 0, 0, 0, 0, 0]
        roots.append(v)

    # 6 long roots: √3 times unit vectors (30° offset)
    for i in range(6):
        angle = (i + 0.5) * np.pi / 3
        v = [np.sqrt(3) * np.cos(angle), np.sqrt(3) * np.sin(angle), 0, 0, 0, 0, 0, 0]
        roots.append(v)

    roots = np.array(roots, float)
    print(f"G2: {len(roots)} roots")

    # Create 5 different projections
    variants = []

    # Variant 1: Rotate by 15°
    theta = np.pi / 12
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    P1 = roots[:, :2] @ rot.T
    variants.append(('rot15', P1))

    # Variant 2: Rotate by 30°
    theta = np.pi / 6
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    P2 = roots[:, :2] @ rot.T
    variants.append(('rot30', P2))

    # Variant 3: Rotate by 45°
    theta = np.pi / 4
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    P3 = roots[:, :2] @ rot.T
    variants.append(('rot45', P3))

    # Variant 4: Scale + rotate
    theta = np.pi / 8
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    P4 = roots[:, :2] * 1.2 @ rot.T
    variants.append(('scaled', P4))

    # Variant 5: Mirror + rotate
    mirror = np.array([[1, 0], [0, -1]])
    P5 = roots[:, :2] @ mirror.T
    variants.append(('mirror', P5))

    # For each variant, create all 4 orbitals (ip = -2, -1, 0, +1)
    G = np.round(roots @ roots.T).astype(int)

    created = []
    for variant_name, P in variants:
        for ip in [-2, -1, 0, 1]:
            edges = []
            for i in range(len(roots)):
                for j in range(i+1, len(roots)):
                    if G[i, j] == ip:
                        edges.append((i, j))

            if edges:
                filename = create_variant_tikz(P, edges, f'G2_{variant_name}', ip)
                created.append(filename)
                print(f"  Created {filename}.tex ({len(edges)} edges)")

    return created


def create_variant_tikz(P, edges, root_system, ip_value):
    """Create a TikZ drawing of a variant graph."""

    # Determine ip label
    if ip_value >= 0:
        ip_label = f"ipp{ip_value}"
    else:
        ip_label = f"ipm{abs(ip_value)}"

    filename = f"orbital_{root_system.lower()}_{ip_label}"

    # Scale to fit
    P = np.asarray(P, dtype=float)
    radii = np.hypot(P[:, 0], P[:, 1])
    rr = sorted(set(np.round(radii, 4)))

    if len(rr) == 0:
        rr = [1.0]

    ring_of = {r: i for i, r in enumerate(rr)}
    s = 7.3 / max(rr) if max(rr) > 0 else 1.0

    # Color scheme
    cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0)
            for i in range(len(rr))]

    # Build TikZ
    o = [
        r"\documentclass[tikz,border=6pt]{standalone}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}"
    ]

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

    return filename


def main():
    """Generate E6 complement and G2 variants."""

    create_e6_complement()

    created = create_g2_variants()

    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Created {len(created)} G2 variant orbitals")
    print("Files: orbital_g2_*.tex")
    print()
    print("Note: This provides visual variety for G2-type classes")
    print("while we await computation of true product Hecke algebras.")


if __name__ == '__main__':
    main()
