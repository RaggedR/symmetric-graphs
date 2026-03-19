#!/usr/bin/env python3
"""
Draw all graphs in a single multi-page or tiled PDF.

Usage: python3 draw_all.py
"""

import subprocess
import sys
import os
import math
from draw_graph import run_gap, find_best_layout, COLORS

GRAPHVIZ_DIR = os.path.expanduser("~/git/graphviz")

ALL_GRAPHS = [
    "cube", "petersen", "heawood", "moebiuskantor",
    "pappus", "dodecahedron", "desargues", "nauru",
    "f26a", "coxeter", "tuttecoxeter",
]


def graph_to_tikz(coords, adj, n, orbits, graph_name, quotient_info,
                  cx, cy, scale=1.0, aut_order=None, aut_structure=None):
    """Generate TikZ commands for one graph, offset to (cx, cy)."""
    n_orbits = len(orbits)
    v_to_orbit = {}
    for oi, orbit in enumerate(orbits):
        for v in orbit:
            v_to_orbit[v] = oi

    lines = []

    # Title and group labels
    q_label = quotient_info.get('subgroup_structure', '?')
    lines.append(f"\\node[graphtitle] at ({cx}, {cy + 5.5 * scale}) {{{graph_name}}};")
    if aut_order and aut_structure:
        lines.append(f"\\node[graphsub] at ({cx}, {cy + 4.9 * scale}) "
                     f"{{$\\mathrm{{Aut}} = {aut_structure}$, "
                     f"$H = {q_label}$}};")
    else:
        lines.append(f"\\node[graphsub] at ({cx}, {cy + 4.9 * scale}) "
                     f"{{$H = {q_label}$, {n_orbits} orbits}};")

    # Vertices
    for v in range(1, n + 1):
        x, y = coords[v]
        oi = v_to_orbit[v]
        ci = oi % len(COLORS)
        vx = cx + x * scale
        vy = cy + y * scale
        lines.append(f"\\node[v, fill=c{ci}!30] "
                     f"(g{graph_name}{v}) at ({vx:.2f}, {vy:.2f}) "
                     f"{{\\tiny {v}}};")

    # Edges
    drawn = set()
    for v in range(1, n + 1):
        for w in adj[v - 1]:
            edge = (min(v, w), max(v, w))
            if edge not in drawn:
                drawn.add(edge)
                oi_v = v_to_orbit[v]
                oi_w = v_to_orbit[w]
                if oi_v == oi_w:
                    style = f"thick, c{oi_v % len(COLORS)}!60"
                else:
                    style = "e"
                lines.append(f"\\draw[{style}] "
                             f"(g{graph_name}{edge[0]}) -- (g{graph_name}{edge[1]});")

    return "\n".join(lines)


def main():
    print("Computing layouts for all graphs...")

    graph_data = []
    for name in ALL_GRAPHS:
        print(f"  {name}...")
        data = run_gap(name, 10)
        n = data['n']
        adj = data['adj']
        results = data['results']
        if not results:
            print(f"    Skipping {name} — no quotients")
            continue
        coords, best = find_best_layout(results, n, adj)

        graph_data.append({
            'name': data['graph'],
            'n': n,
            'adj': adj,
            'coords': coords,
            'orbits': best['orbits'],
            'quotient': best,
            'aut_order': data.get('aut_order'),
            'aut_structure': data.get('aut_structure'),
        })

    # Layout: 4 columns (smaller scale for more graphs)
    ncols = 4
    col_width = 11.0
    row_height = 13.0
    scale = 0.75

    positions = []
    for i in range(len(graph_data)):
        row = i // ncols
        col = i % ncols
        # Center the last row if incomplete
        items_in_row = min(ncols, len(graph_data) - row * ncols)
        col_offset = (ncols - items_in_row) * col_width / 2.0
        cx = col * col_width + col_offset + col_width / 2.0
        cy = -row * row_height
        positions.append((cx, cy))

    # Generate TikZ
    lines = []
    lines.append(r"\documentclass[border=20pt]{standalone}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{tikz}")
    lines.append("")

    # Colors
    for i in range(min(15, len(COLORS))):
        lines.append(f"\\definecolor{{c{i}}}{{HTML}}{{{COLORS[i]}}}")
    lines.append("")

    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}[")
    lines.append(r"  v/.style={circle, draw, thick, minimum size=5mm, inner sep=0pt},")
    lines.append(r"  e/.style={thick, gray!40},")
    lines.append(r"  graphtitle/.style={font=\large\bfseries},")
    lines.append(r"  graphsub/.style={font=\footnotesize, text=gray},")
    lines.append(r"]")
    lines.append("")

    # Main title
    total_width = ncols * col_width
    lines.append(f"\\node[font=\\Huge\\bfseries] at ({total_width/2}, 7.5) "
                 f"{{Symmetric Graph Gallery}};")
    lines.append(f"\\node[font=\\normalsize, text=gray] at ({total_width/2}, 6.5) "
                 f"{{Quotient-driven layouts --- cyclic and cylindrical symmetry}};")
    lines.append("")

    for i, gd in enumerate(graph_data):
        cx, cy = positions[i]
        lines.append(f"% === {gd['name']} ===")
        tikz = graph_to_tikz(
            gd['coords'], gd['adj'], gd['n'],
            gd['orbits'], gd['name'], gd['quotient'],
            cx, cy, scale=scale,
            aut_order=gd.get('aut_order'),
            aut_structure=gd.get('aut_structure'),
        )
        lines.append(tikz)
        lines.append("")

    lines.append(r"\end{tikzpicture}")
    lines.append(r"\end{document}")

    # Write and compile
    tex_file = os.path.join(GRAPHVIZ_DIR, "gallery.tex")
    pdf_file = os.path.join(GRAPHVIZ_DIR, "gallery.pdf")

    with open(tex_file, 'w') as f:
        f.write("\n".join(lines))
    print(f"\nWrote {tex_file}")

    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "gallery.tex"],
        capture_output=True, text=True, cwd=GRAPHVIZ_DIR
    )

    if result.returncode != 0:
        print("LaTeX failed!")
        # Show last few lines of log
        log_lines = result.stdout.split('\n')
        for line in log_lines[-20:]:
            if line.strip():
                print(f"  {line}")
        sys.exit(1)

    # Cleanup
    for ext in ['aux', 'log']:
        try:
            os.remove(os.path.join(GRAPHVIZ_DIR, f"gallery.{ext}"))
        except FileNotFoundError:
            pass

    print(f"Compiled {pdf_file}")
    subprocess.run(["open", pdf_file])


if __name__ == "__main__":
    main()
