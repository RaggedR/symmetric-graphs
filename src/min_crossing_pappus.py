#!/usr/bin/env python3
"""
Find a minimum-crossing drawing of the Pappus graph (cr = 5).

Uses simulated annealing to optimize vertex positions, starting from
our symmetry-aware layout (7 crossings) and relaxing to find 5.

The Pappus graph: 18 vertices, 27 edges, cubic, bipartite.
Vertices 1-9 are "points", 10-18 are "lines" of the Pappus configuration.
"""

import math
import random
import subprocess
import os
import sys

# Pappus graph adjacency (1-indexed)
# Standard construction as Levi graph of Pappus configuration
PAPPUS_EDGES = [
    (1, 2), (1, 6), (1, 18),
    (2, 3), (2, 9),
    (3, 4), (3, 14),
    (4, 5), (4, 11),
    (5, 6), (5, 16),
    (6, 7),
    (7, 8), (7, 12),
    (8, 9), (8, 15),
    (9, 10),
    (10, 11), (10, 17),
    (11, 12),
    (12, 13),
    (13, 14), (13, 18),
    (14, 15),
    (15, 16),
    (16, 17),
    (17, 18),
]

N = 18

def build_adj():
    adj = [[] for _ in range(N)]
    for u, v in PAPPUS_EDGES:
        adj[u-1].append(v)
        adj[v-1].append(u)
    return adj

ADJ = build_adj()


def segments_cross(p1, p2, p3, p4):
    """Test if segments p1-p2 and p3-p4 properly cross."""
    def cross2d(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    d1 = cross2d(p3, p4, p1)
    d2 = cross2d(p3, p4, p2)
    d3 = cross2d(p1, p2, p3)
    d4 = cross2d(p1, p2, p4)

    denom1 = d1 - d2
    denom2 = d3 - d4

    if abs(denom1) < 1e-12 or abs(denom2) < 1e-12:
        return False

    t = d1 / denom1
    s = d3 / denom2

    EPS = 1e-9
    return EPS < t < 1 - EPS and EPS < s < 1 - EPS


def count_crossings(coords):
    """Count edge crossings."""
    crossings = 0
    for i in range(len(PAPPUS_EDGES)):
        u1, v1 = PAPPUS_EDGES[i]
        p1, p2 = coords[u1], coords[v1]
        for j in range(i + 1, len(PAPPUS_EDGES)):
            u2, v2 = PAPPUS_EDGES[j]
            if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                continue
            p3, p4 = coords[u2], coords[v2]
            if segments_cross(p1, p2, p3, p4):
                crossings += 1
    return crossings


def find_crossing_edges(coords):
    """Return list of crossing edge pairs."""
    crossing_pairs = []
    for i in range(len(PAPPUS_EDGES)):
        u1, v1 = PAPPUS_EDGES[i]
        p1, p2 = coords[u1], coords[v1]
        for j in range(i + 1, len(PAPPUS_EDGES)):
            u2, v2 = PAPPUS_EDGES[j]
            if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                continue
            p3, p4 = coords[u2], coords[v2]
            if segments_cross(p1, p2, p3, p4):
                crossing_pairs.append(((u1, v1), (u2, v2)))
    return crossing_pairs


def edge_length_variance(coords):
    """Measure how uniform edge lengths are (lower = more uniform)."""
    lengths = []
    for u, v in PAPPUS_EDGES:
        dx = coords[u][0] - coords[v][0]
        dy = coords[u][1] - coords[v][1]
        lengths.append(math.sqrt(dx*dx + dy*dy))
    mean = sum(lengths) / len(lengths)
    return sum((l - mean)**2 for l in lengths) / len(lengths)


def min_vertex_distance(coords):
    """Minimum distance between any two vertices."""
    min_d = float('inf')
    for i in range(1, N + 1):
        for j in range(i + 1, N + 1):
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            d = math.sqrt(dx*dx + dy*dy)
            if d < min_d:
                min_d = d
    return min_d


def initial_coords_circle():
    """Start with vertices on a circle."""
    coords = {}
    for v in range(1, N + 1):
        theta = 2 * math.pi * (v - 1) / N
        coords[v] = (4.0 * math.cos(theta), 4.0 * math.sin(theta))
    return coords


def initial_coords_from_gap():
    """Use our symmetry-aware layout as starting point."""
    # From auto-pappus.tex — the 7-crossing layout
    coords = {
        1:  ( 0.000,  4.500),
        2:  (-3.897,  2.250),
        3:  (-3.897, -2.250),
        4:  (-0.000, -4.500),
        5:  ( 3.897, -2.250),
        6:  ( 3.897,  2.250),
        7:  ( 0.000,  3.038),
        8:  ( 0.000,  1.575),
        9:  (-2.631,  1.519),
        10: (-1.364,  0.787),
        11: (-2.631, -1.519),
        12: (-1.364, -0.788),
        13: (-0.000, -1.575),
        14: (-0.000, -3.038),
        15: ( 1.364, -0.788),
        16: ( 2.631, -1.519),
        17: ( 1.364,  0.788),
        18: ( 2.631,  1.519),
    }
    return coords


def simulated_annealing(coords, max_iter=2_000_000, target=5):
    """Simulated annealing to minimize crossings."""
    best_coords = dict(coords)
    best_crossings = count_crossings(coords)
    current_crossings = best_crossings
    print(f"  Start: {current_crossings} crossings")

    if best_crossings <= target:
        return best_coords, best_crossings

    T = 2.0
    T_min = 0.001
    alpha = 0.999997

    for iteration in range(max_iter):
        # Pick a random vertex
        v = random.randint(1, N)

        # Perturbation size decreases with temperature
        scale = T * 0.5
        dx = random.gauss(0, scale)
        dy = random.gauss(0, scale)

        old_pos = coords[v]
        new_pos = (old_pos[0] + dx, old_pos[1] + dy)

        # Keep in reasonable bounds
        if abs(new_pos[0]) > 8 or abs(new_pos[1]) > 8:
            continue

        coords[v] = new_pos
        new_crossings = count_crossings(coords)
        delta = new_crossings - current_crossings

        # Accept or reject
        if delta <= 0 or random.random() < math.exp(-delta / T):
            current_crossings = new_crossings
            if current_crossings < best_crossings:
                best_crossings = current_crossings
                best_coords = dict(coords)
                print(f"  iter {iteration}: {best_crossings} crossings (T={T:.4f})")
                if best_crossings <= target:
                    print(f"  Target reached!")
                    return best_coords, best_crossings
        else:
            coords[v] = old_pos

        T = max(T_min, T * alpha)

        if iteration % 500000 == 0 and iteration > 0:
            print(f"  iter {iteration}: best={best_crossings}, current={current_crossings}, T={T:.6f}")

    return best_coords, best_crossings


def refine_layout(coords, target=5):
    """Post-SA refinement: gradient-free local optimization."""
    crossings = count_crossings(coords)
    if crossings <= target:
        return coords, crossings

    # Try small adjustments to each vertex
    improved = True
    step = 0.1
    while improved and step > 0.001:
        improved = False
        for v in range(1, N + 1):
            for dx, dy in [(step, 0), (-step, 0), (0, step), (0, -step),
                           (step, step), (step, -step), (-step, step), (-step, -step)]:
                old = coords[v]
                coords[v] = (old[0] + dx, old[1] + dy)
                new_cross = count_crossings(coords)
                if new_cross < crossings:
                    crossings = new_cross
                    improved = True
                    if crossings <= target:
                        return coords, crossings
                else:
                    coords[v] = old
        if not improved:
            step *= 0.5
            improved = True  # try again with smaller step

    return coords, crossings


def beautify(coords, crossings_target=5, iterations=5000):
    """
    Improve aesthetics while preserving crossing count.
    Uses spring-electrical model with crossing constraint.
    """
    base_crossings = count_crossings(coords)

    for iteration in range(iterations):
        # Compute forces
        forces = {v: [0.0, 0.0] for v in range(1, N + 1)}

        # Spring forces (edges want ideal length)
        ideal_len = 2.0
        spring_k = 0.01
        for u, v in PAPPUS_EDGES:
            dx = coords[v][0] - coords[u][0]
            dy = coords[v][1] - coords[u][1]
            dist = math.sqrt(dx*dx + dy*dy) + 1e-10
            force = spring_k * (dist - ideal_len) / dist
            forces[u][0] += force * dx
            forces[u][1] += force * dy
            forces[v][0] -= force * dx
            forces[v][1] -= force * dy

        # Repulsion between all pairs
        repulsion_k = 0.5
        for i in range(1, N + 1):
            for j in range(i + 1, N + 1):
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]
                dist = math.sqrt(dx*dx + dy*dy) + 1e-10
                force = -repulsion_k / (dist * dist)
                forces[i][0] += force * dx / dist
                forces[i][1] += force * dy / dist
                forces[j][0] -= force * dx / dist
                forces[j][1] -= force * dy / dist

        # Apply forces, but reject if crossings increase
        step = max(0.01, 0.1 * (1 - iteration / iterations))
        new_coords = {}
        for v in range(1, N + 1):
            new_coords[v] = (
                coords[v][0] + step * forces[v][0],
                coords[v][1] + step * forces[v][1],
            )

        if count_crossings(new_coords) <= base_crossings:
            coords = new_coords

    return coords


def generate_tikz(coords, crossing_edges=None):
    """Generate TikZ for the minimum-crossing drawing."""
    lines = []
    lines.append(r"\documentclass[border=15pt]{standalone}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\usetikzlibrary{arrows.meta}")
    lines.append("")
    lines.append(r"\definecolor{crossedge}{HTML}{E53935}")
    lines.append(r"\definecolor{normaledge}{HTML}{4A4A4A}")
    lines.append("")
    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}[")
    lines.append(r"  vertex/.style={circle, draw, thick, minimum size=6.5mm, inner sep=0pt,")
    lines.append(r"                 font=\scriptsize\bfseries, fill=blue!12},")
    lines.append(r"  edge/.style={thick, normaledge!60},")
    lines.append(r"  crossedge/.style={thick, crossedge!80},")
    lines.append(r"  title/.style={font=\Large\bfseries},")
    lines.append(r"  subtitle/.style={font=\small, text=gray},")
    lines.append(r"]")
    lines.append("")

    n_crossings = count_crossings(coords)
    lines.append(f"\\node[title] at (0, 6.5) {{Pappus graph}};")
    lines.append(f"\\node[subtitle] at (0, 5.8) {{Minimum crossing number: {n_crossings}}};")
    lines.append("")

    # Determine which edges are involved in crossings
    crossing_edge_set = set()
    if crossing_edges:
        for (e1, e2) in crossing_edges:
            crossing_edge_set.add(e1)
            crossing_edge_set.add(e2)

    # Vertices (define first so edges can reference them)
    for v in range(1, N + 1):
        x, y = coords[v]
        lines.append(f"\\node[vertex] (v{v}) at ({x:.3f}, {y:.3f}) {{{v}}};")

    lines.append("")

    # Edges
    for u, v in PAPPUS_EDGES:
        edge = (min(u, v), max(u, v))
        style = "crossedge" if edge in crossing_edge_set else "edge"
        lines.append(f"\\draw[{style}] (v{u}) -- (v{v});")

    lines.append("")
    lines.append(r"\end{tikzpicture}")
    lines.append(r"\end{document}")
    return "\n".join(lines)


def main():
    random.seed(42)

    print("=== Pappus graph: finding cr(G) = 5 embedding ===\n")

    # Try multiple starting points
    best_overall = None
    best_overall_crossings = float('inf')

    starts = [
        ("symmetry-aware (7-crossing)", initial_coords_from_gap),
        ("circle", initial_coords_circle),
    ]

    for name, init_fn in starts:
        print(f"\n--- Starting from {name} ---")
        for trial in range(5):
            random.seed(42 + trial * 137)
            coords = init_fn()
            coords, crossings = simulated_annealing(coords, max_iter=2_000_000, target=5)

            if crossings > 5:
                coords, crossings = refine_layout(coords, target=5)

            if crossings < best_overall_crossings:
                best_overall_crossings = crossings
                best_overall = dict(coords)
                if crossings <= 5:
                    break

        if best_overall_crossings <= 5:
            break

    print(f"\n=== Best: {best_overall_crossings} crossings ===")

    if best_overall_crossings <= 5:
        # Beautify while preserving crossing count
        print("Beautifying layout...")
        best_overall = beautify(best_overall, crossings_target=best_overall_crossings)
        final_crossings = count_crossings(best_overall)
        print(f"After beautification: {final_crossings} crossings")

    # Find which edges cross
    crossing_pairs = find_crossing_edges(best_overall)
    print(f"\nCrossing edge pairs:")
    for (e1, e2) in crossing_pairs:
        print(f"  {e1} x {e2}")

    # Generate output
    tikz = generate_tikz(best_overall, crossing_pairs)

    out_dir = os.path.join(os.path.dirname(__file__), "individual_graphs")
    os.makedirs(out_dir, exist_ok=True)
    tex_file = os.path.join(out_dir, "pappus-min-crossing.tex")
    pdf_file = os.path.join(out_dir, "pappus-min-crossing.pdf")

    with open(tex_file, 'w') as f:
        f.write(tikz)

    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "pappus-min-crossing.tex"],
        capture_output=True, text=True, cwd=out_dir
    )

    if result.returncode != 0:
        print("LaTeX compilation failed!")
        print(result.stdout[-500:])
        sys.exit(1)

    for ext in ['aux', 'log']:
        try:
            os.remove(os.path.join(out_dir, f"pappus-min-crossing.{ext}"))
        except FileNotFoundError:
            pass

    print(f"\nOutput: {pdf_file}")
    if "--no-open" not in sys.argv:
        subprocess.run(["open", pdf_file])


if __name__ == "__main__":
    main()
