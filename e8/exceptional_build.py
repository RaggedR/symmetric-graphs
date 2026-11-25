"""
Build the remaining exceptional root systems: F4, G2, A2, A3
for the complete trading card collection.
"""
import numpy as np, colorsys, itertools
import ladder as L

def draw_root_system(roots, name, title, h, style='roots'):
    """Generic drawer for root systems"""
    # Get Coxeter plane
    simple = L.simple_roots(roots)
    order, u1, u2 = L.coxeter_plane(simple, h)
    print(f"{name}: {len(roots)} roots, Coxeter order {order} (expect {h})")

    P = L.project(roots, u1, u2)
    rings = L.ring_histogram(P)
    print(f"  ring histogram: {rings}")

    # Save
    np.save(f"/tmp/{name}_R.npy", roots)
    np.save(f"/tmp/{name}_P.npy", P)

    # Draw
    if style == 'roots':
        # For root systems, draw edges between roots at angle ~120° (ip = -1)
        G = np.round(roots @ roots.T).astype(int)
        edges = [(i, j) for i, j in itertools.combinations(range(len(roots)), 2)
                 if G[i, j] == -1]
    else:
        # For weight systems, use ip = +1
        G = np.round(roots @ roots.T).astype(int)
        edges = [(i, j) for i, j in itertools.combinations(range(len(roots)), 2)
                 if G[i, j] == 1]

    radii = np.round(np.hypot(P[:, 0], P[:, 1]), 4)
    rr = sorted(set(radii))
    ring_of = {r: i for i, r in enumerate(rr)}
    s = 7.3 / max(rr) if max(rr) > 0 else 1.0
    cols = [colorsys.hsv_to_rgb(0.05 + 0.70*i/max(1, len(rr)-1), 0.80, 1.0)
            for i in range(len(rr))]

    o = [r"\documentclass[tikz,border=6pt]{standalone}", r"\usepackage{tikz}",
         r"\begin{document}", r"\begin{tikzpicture}"]
    for i, (r, g, b) in enumerate(cols):
        o.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
    o += [r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",
          r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);"]

    for i, j in edges:
        ro = max(ring_of[radii[i]], ring_of[radii[j]])
        o.append(fr"\draw[ring{ro},line width=0.5pt,opacity=0.55] "
                 fr"({P[i,0]*s:.3f},{P[i,1]*s:.3f})--({P[j,0]*s:.3f},{P[j,1]*s:.3f});")

    for k in range(len(roots)):
        o.append(fr"\fill[ring{ring_of[radii[k]]}] ({P[k,0]*s:.3f},{P[k,1]*s:.3f}) circle (3.0pt);")

    o.append(fr"\node[white!80!black,font=\small\itshape] at (0,-8.4) {{{title}}};")
    o += [r"\end{tikzpicture}", r"\end{document}"]

    fname = f"/Users/miriam/symmetric-graphs/e8/{name}_roots"
    open(f"{fname}.tex", "w").write("\n".join(o))
    print(f"  {name}: {len(edges)} edges written to {fname}.tex")


# F4: 48 roots
def build_f4():
    """F4 roots in R^4 (embedded in R^8)
    - 24 roots: (±1, ±1, 0, 0) and all permutations (D4 type)
    - 8 roots: (±1, 0, 0, 0) and permutations
    - 16 roots: (±1/2, ±1/2, ±1/2, ±1/2) with even # of minuses
    Total: 48 roots
    """
    roots = []

    # 24 D4-type roots: (±1, ±1, 0, 0) and permutations
    for i, j in itertools.combinations(range(4), 2):
        for si, sj in itertools.product([1, -1], repeat=2):
            v = [0]*8
            v[i] = si
            v[j] = sj
            roots.append(v)

    # 8 roots: (±1, 0, 0, 0) and permutations
    for i in range(4):
        for s in [1, -1]:
            v = [0]*8
            v[i] = s
            roots.append(v)

    # 16 roots: (±1/2)^4 ALL combinations
    for signs in itertools.product([1, -1], repeat=4):
        v = [s*0.5 for s in signs] + [0, 0, 0, 0]
        roots.append(v)

    roots = np.array(roots, float)
    return roots


# G2: 12 roots
def build_g2():
    """G2 roots: 6 short + 6 long in hexagonal pattern"""
    roots = []

    # Short roots: hexagonal
    for k in range(6):
        angle = k * np.pi / 3
        v = [np.cos(angle), np.sin(angle)] + [0]*6
        roots.append(v)

    # Long roots: dual hexagonal (rotated 30°, scaled by sqrt(3))
    for k in range(6):
        angle = k * np.pi / 3 + np.pi / 6
        v = [np.sqrt(3) * np.cos(angle), np.sqrt(3) * np.sin(angle)] + [0]*6
        roots.append(v)

    roots = np.array(roots, float)
    # Normalize to length^2 = 2 for long roots
    roots = roots * np.sqrt(2.0 / 3.0)

    return roots


# A2: 6 roots (triangle)
def build_a2():
    """A2 roots: 6 roots forming hexagon (±(1,0,0), ±(0,1,0), ±(0,0,1) with sum=0)"""
    roots = []

    # The 6 roots: differences of basis vectors
    # e1-e2, e2-e3, e3-e1 and their negatives
    roots = [
        [1, -1, 0, 0, 0, 0, 0, 0],
        [0, 1, -1, 0, 0, 0, 0, 0],
        [-1, 0, 1, 0, 0, 0, 0, 0],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [1, 0, -1, 0, 0, 0, 0, 0],
    ]

    roots = np.array(roots, float)
    # Normalize to length^2 = 2
    roots = roots * np.sqrt(2.0 / 2.0)

    return roots


# A3: 12 roots (tetrahedron)
def build_a3():
    """A3 roots: 12 roots, differences of 4 basis vectors"""
    roots = []

    # ei - ej for all i != j, i,j in {0,1,2,3}
    for i in range(4):
        for j in range(4):
            if i != j:
                v = [0]*8
                v[i] = 1
                v[j] = -1
                roots.append(v)

    roots = np.array(roots, float)
    # Already have length^2 = 2

    return roots


if __name__ == "__main__":
    print("Building exceptional root systems...\n")

    # F4
    f4 = build_f4()
    assert len(f4) == 48
    draw_root_system(f4, "f4", r"$F_4$: 48 roots, 12-fold (The Tesseract Guardian)", 12)

    # G2
    g2 = build_g2()
    assert len(g2) == 12
    draw_root_system(g2, "g2", r"$G_2$: 12 roots, 6-fold (The Hexagon)", 6)

    # A2
    a2 = build_a2()
    assert len(a2) == 6
    draw_root_system(a2, "a2", r"$A_2$: 6 roots, 6-fold (The Triangle / QCD)", 3)

    # A3
    a3 = build_a3()
    assert len(a3) == 12
    draw_root_system(a3, "a3", r"$A_3$: 12 roots, 4-fold (The Tetrahedron)", 4)

    print("\nAll exceptional root systems built!")
