#!/usr/bin/env python3
"""
Draw the collinearity graph of the twisted triality hexagon GH(2,8).

819 vertices, valency 18, distance-regular, diameter 3.
Automorphism group: ³D₄(2), order 211,341,312.

Layout: 117 concentric rings of 7 vertices each (C₇ quotient),
with angular offsets optimized to reduce visual clutter.
Edges colored by inter-ring distance to reveal structure.
"""

import math
import subprocess
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_adjacency(path):
    adj = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or ':' not in line:
                continue
            v_str, nbrs_str = line.split(':', 1)
            v = int(v_str)
            nbrs = [int(x) for x in nbrs_str.split(',') if x.strip()]
            adj[v] = nbrs
    return adj


def load_orbits(path):
    orbits = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or ':' not in line:
                continue
            _, verts_str = line.split(':', 1)
            orbit = [int(x) for x in verts_str.split(',') if x.strip()]
            orbits.append(orbit)
    return orbits


def reorder_orbits_by_connectivity(orbits, adj):
    """BFS reorder: heavily-connected orbits go on adjacent rings."""
    n = len(orbits)
    v_to_orbit = {}
    for oi, orbit in enumerate(orbits):
        for v in orbit:
            v_to_orbit[v] = oi

    edge_count = [[0] * n for _ in range(n)]
    for oi, orbit in enumerate(orbits):
        for v in orbit:
            for nb in adj.get(v, []):
                oj = v_to_orbit.get(nb)
                if oj is not None and oj != oi:
                    edge_count[oi][oj] += 1

    total_edges = [sum(row) for row in edge_count]
    start = max(range(n), key=lambda i: total_edges[i])

    order = [start]
    visited = {start}
    while len(order) < n:
        best = -1
        best_score = -1
        for oi in range(n):
            if oi in visited:
                continue
            score = sum(edge_count[oi][oj] for oj in order)
            if score > best_score:
                best_score = score
                best = oi
        if best >= 0:
            order.append(best)
            visited.add(best)
        else:
            for oi in range(n):
                if oi not in visited:
                    order.append(oi)
                    visited.add(oi)

    return [orbits[i] for i in order]


def layout_concentric_optimized(orbits, adj, n_verts):
    """Place vertices on concentric rings with greedy angular optimization."""
    n_rings = len(orbits)
    R_max = 8.0
    R_min = 0.5

    coords = {}

    for ring_idx, orbit in enumerate(orbits):
        k = len(orbit)
        r = R_max - ring_idx * (R_max - R_min) / max(1, n_rings - 1)

        if ring_idx == 0:
            for vi, v in enumerate(orbit):
                theta = math.radians(90 + vi * (360.0 / k))
                coords[v] = (r * math.cos(theta), r * math.sin(theta))
        else:
            best_offset = 0
            best_cost = float('inf')
            for offset in range(k):
                cost = 0
                for vi, v in enumerate(orbit):
                    theta = math.radians(90 + ((vi + offset) % k) * (360.0 / k))
                    x, y = r * math.cos(theta), r * math.sin(theta)
                    for nb in adj.get(v, []):
                        if nb in coords:
                            nx, ny = coords[nb]
                            cost += math.hypot(x - nx, y - ny)
                if cost < best_cost:
                    best_cost = cost
                    best_offset = offset

            for vi, v in enumerate(orbit):
                theta = math.radians(90 + ((vi + best_offset) % k) * (360.0 / k))
                coords[v] = (r * math.cos(theta), r * math.sin(theta))

    return coords


def generate_tikz(coords, adj, orbits):
    """Generate TikZ with edges colored by inter-ring distance."""
    v_to_ring = {}
    for ri, orbit in enumerate(orbits):
        for v in orbit:
            v_to_ring[v] = ri

    n_rings = len(orbits)

    # Classify edges by ring gap
    edges_by_gap = {}  # gap -> list of (v, w)
    drawn = set()
    for v, nbrs in adj.items():
        for w in nbrs:
            edge = (min(v, w), max(v, w))
            if edge not in drawn:
                drawn.add(edge)
                gap = abs(v_to_ring[v] - v_to_ring[w])
                edges_by_gap.setdefault(gap, []).append(edge)

    gaps = sorted(edges_by_gap.keys())
    print(f"  Edge ring-gaps: {[(g, len(edges_by_gap[g])) for g in gaps[:10]]}...")
    total_edges = sum(len(v) for v in edges_by_gap.values())
    print(f"  Total edges: {total_edges}")

    lines = []
    lines.append(r"\documentclass[border=10pt]{standalone}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append("")

    # Color palette for rings (blue outer to red inner)
    ring_palette = [
        ("0033CC", "ring0"),
        ("0066FF", "ring1"),
        ("0099FF", "ring2"),
        ("00CCCC", "ring3"),
        ("00CC66", "ring4"),
        ("33CC00", "ring5"),
        ("99CC00", "ring6"),
        ("CCCC00", "ring7"),
        ("FF9900", "ring8"),
        ("FF6600", "ring9"),
        ("FF3300", "ring10"),
        ("CC0000", "ring11"),
    ]
    for hex_c, name in ring_palette:
        lines.append(f"\\definecolor{{{name}}}{{HTML}}{{{hex_c}}}")

    # Edge colors by gap
    lines.append(r"\definecolor{gap0}{HTML}{E53935}")   # same ring: red
    lines.append(r"\definecolor{gap1}{HTML}{1E88E5}")   # adjacent rings: blue
    lines.append(r"\definecolor{gapfar}{HTML}{CCCCCC}") # far: light grey
    lines.append("")

    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}[")
    lines.append(r"  title/.style={font=\Large\bfseries},")
    lines.append(r"  subtitle/.style={font=\small, text=gray},")
    lines.append(r"]")
    lines.append("")

    lines.append(r"\node[title] at (0, 10) {Collinearity graph of $\mathrm{GH}(2,8)$};")
    lines.append(r"\node[subtitle] at (0, 9.3) {819 vertices, valency 18, "
                 r"$\mathrm{Aut} = {}^3\!D_4(2)$, $|\mathrm{Aut}| = 211{,}341{,}312$};")
    lines.append(r"\node[subtitle] at (0, 8.7) {Distance-regular, diameter 3 "
                 r"--- the twisted triality hexagon};")
    lines.append("")

    # Draw edges: far edges first (background), close edges on top
    for gap in sorted(edges_by_gap.keys(), reverse=True):
        edge_list = edges_by_gap[gap]
        if gap == 0:
            style = "line width=0.3pt, gap0, opacity=0.6"
        elif gap <= 3:
            style = "line width=0.15pt, gap1, opacity=0.25"
        elif gap <= 10:
            style = "line width=0.1pt, gapfar, opacity=0.12"
        else:
            style = "line width=0.08pt, gapfar, opacity=0.07"

        for v, w in edge_list:
            x1, y1 = coords[v]
            x2, y2 = coords[w]
            lines.append(f"\\draw[{style}] ({x1:.4f},{y1:.4f}) -- ({x2:.4f},{y2:.4f});")

    lines.append("")

    # Vertices — larger dots, colored by ring
    for v, (x, y) in sorted(coords.items()):
        ring = v_to_ring[v]
        palette_idx = min(11, int(ring * 12 / n_rings))
        color = f"ring{palette_idx}"
        lines.append(f"\\fill[{color}] ({x:.4f},{y:.4f}) circle (1.0pt);")

    lines.append("")
    lines.append(r"\end{tikzpicture}")
    lines.append(r"\end{document}")

    return "\n".join(lines)


def main():
    adj_path = os.path.join(SCRIPT_DIR, "triality_adj.txt")
    orbits_path = os.path.join(SCRIPT_DIR, "triality_orbits.txt")

    print("Loading adjacency data...")
    adj = load_adjacency(adj_path)
    print(f"  {len(adj)} vertices, valency {len(next(iter(adj.values())))}")

    print("Loading orbit data...")
    orbits = load_orbits(orbits_path)
    print(f"  {len(orbits)} orbits of size {len(orbits[0])}")

    print("Reordering orbits by connectivity...")
    orbits = reorder_orbits_by_connectivity(orbits, adj)

    print("Computing layout (117 rings x 7 vertices, optimizing offsets)...")
    coords = layout_concentric_optimized(orbits, adj, 819)

    print("Generating TikZ...")
    tikz = generate_tikz(coords, adj, orbits)

    tex_path = os.path.join(SCRIPT_DIR, "triality_hexagon.tex")
    with open(tex_path, 'w') as f:
        f.write(tikz)
    print(f"  Written: {tex_path}")

    print("Compiling PDF...")
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "triality_hexagon.tex"],
        capture_output=True, text=True, cwd=SCRIPT_DIR, timeout=120
    )
    if result.returncode == 0:
        pdf_path = os.path.join(SCRIPT_DIR, "triality_hexagon.pdf")
        print(f"  Done: {pdf_path}")
        subprocess.run(["open", pdf_path])
    else:
        print("  pdflatex failed!")
        print(result.stdout[-500:] if result.stdout else "")
        print(result.stderr[-500:] if result.stderr else "")


if __name__ == "__main__":
    main()
