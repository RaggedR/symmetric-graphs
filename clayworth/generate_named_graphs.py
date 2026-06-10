#!/usr/bin/env python3
"""
Generate named graph images for the Lean CSS surface codes project.

Produces 6 JPGs with consistent minimal style:
  - Black dot vertices, thin grey edges, no labels
  - Title: "Name (Foster ID)"
  - Subtitle: Aut group, |Aut|, vertex count, girth

Usage:
    python3 generate_named_graphs.py           # all 6 graphs
    python3 generate_named_graphs.py heawood   # single graph by display name
"""

import sys
import os
import math
import subprocess
import shutil
import tempfile

# Import layout engine from draw_graph.py
sys.path.insert(0, os.path.expanduser("~/git/graphviz/src"))
from draw_graph import run_gap, find_best_layout, count_all_crossings

NAMED_GRAPHS = [
    {
        "gap_name": "heawood",
        "display_name": "Heawood",
        "foster_id": "F014A",
        "girth": 6,
        "max_deg": 10,
        "output": "heawood-F014A.jpg",
    },
    {
        "gap_name": "moebiuskantor",
        "display_name": "Möbius-Kantor",
        "foster_id": "F016A",
        "girth": 6,
        "max_deg": 10,
        "output": "mobius-kantor-F016A.jpg",
    },
    {
        "gap_name": "psl2_7_c3",
        "display_name": "Klein",
        "foster_id": "F056A",
        "girth": 7,
        "max_deg": 4,
        "output": "klein-F056A.jpg",
    },
    {
        "gap_name": "psl2_13",
        "display_name": "Zhou-3",
        "foster_id": "F182A",
        "girth": 7,
        "max_deg": 4,
        "output": "zhou3-F182A.jpg",
    },
    {
        "gap_name": "tuttecage12",
        "display_name": "Tutte 12-cage",
        "foster_id": "F126A",
        "girth": 12,
        "max_deg": 4,
        "output": "tutte-12cage.jpg",
    },
    {
        "gap_name": "gh22",
        "display_name": "Langer",
        "foster_id": None,
        "girth": 6,
        "max_deg": 6,
        "output": "langer-63v.jpg",
    },
]

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "named_graphs")


def latex_safe_group(aut_structure):
    """Convert GAP StructureDescription to LaTeX-safe math string."""
    s = aut_structure
    # G_2(2) → G_{2}(2)  (proper subscript)
    s = s.replace("G_2(2)", "G_{2}(2)")
    # PSL(a,b) → \\mathrm{PSL}(a,b)
    for name in ["PSL", "PGL", "GL", "SL", "PSU"]:
        s = s.replace(name, f"\\mathrm{{{name}}}")
    # " x " → \\times
    s = s.replace(" x ", " \\times ")
    # " : " → \\rtimes
    s = s.replace(" : ", " \\rtimes ")
    return s


def generate_tikz(coords, adj, n, spec, aut_order, aut_structure):
    """Generate TikZ for a named graph with consistent minimal style."""
    # Compute bounding box
    xs = [x for x, y in coords.values()]
    ys = [y for x, y in coords.values()]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    # Node size scaled by graph size
    dot_size = max(1.2, 3.5 - n / 60.0)

    # Edge opacity scales with density
    edge_opacity = max(15, min(40, 60 - n // 3))

    # Title
    name = spec["display_name"]
    foster = spec["foster_id"]
    title = f"{name} ({foster})" if foster else name

    # Subtitle
    aut_latex = latex_safe_group(aut_structure)
    subtitle = (
        f"$\\mathrm{{Aut}} \\cong {aut_latex}$, "
        f"$|\\mathrm{{Aut}}| = {aut_order}$, "
        f"{n} vertices, girth {spec['girth']}"
    )

    title_y = max_y + 1.2
    sub_y = title_y - 0.65

    lines = []
    lines.append(r"\documentclass[border=12pt]{standalone}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{tikz}")
    lines.append("")
    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}[")
    lines.append(f"  edge/.style={{thin, black!{edge_opacity}}},")
    lines.append(r"  title/.style={font=\Large\bfseries},")
    lines.append(r"  subtitle/.style={font=\small, text=gray},")
    lines.append(r"]")
    lines.append("")

    # Title and subtitle
    lines.append(f"\\node[title] at (0, {title_y:.2f}) {{{title}}};")
    lines.append(f"\\node[subtitle] at (0, {sub_y:.2f}) {{{subtitle}}};")
    lines.append("")

    # Edges (drawn first, behind vertices)
    drawn = set()
    for v in range(1, n + 1):
        for w in adj[v - 1]:
            edge = (min(v, w), max(v, w))
            if edge not in drawn:
                drawn.add(edge)
                x1, y1 = coords[edge[0]]
                x2, y2 = coords[edge[1]]
                lines.append(f"\\draw[edge] ({x1:.3f}, {y1:.3f}) -- ({x2:.3f}, {y2:.3f});")
    lines.append("")

    # Vertices (drawn on top)
    for v in range(1, n + 1):
        x, y = coords[v]
        lines.append(f"\\fill[black] ({x:.3f}, {y:.3f}) circle ({dot_size:.1f}pt);")

    lines.append("")
    lines.append(r"\end{tikzpicture}")
    lines.append(r"\end{document}")

    return "\n".join(lines)


def build_graph(spec, build_dir):
    """Run the full pipeline for one graph: GAP → layout → TikZ → PDF → JPG."""
    name = spec["display_name"]
    gap_name = spec["gap_name"]
    output_file = spec["output"]

    print(f"\n{'='*60}")
    print(f"  {name} ({gap_name})")
    print(f"{'='*60}")

    # 1. Run GAP
    print(f"  Running GAP...")
    data = run_gap(gap_name, spec["max_deg"])

    n = data["n"]
    adj = data["adj"]
    results = data["results"]
    aut_order = data.get("aut_order", "?")
    aut_structure = data.get("aut_structure", "?")

    print(f"  {n} vertices, |Aut| = {aut_order} ({aut_structure})")
    print(f"  {len(results)} quotient(s) found")

    if not results:
        print(f"  ERROR: No quotients found for {name}!")
        return False

    # 2. Compute layout
    coords, best = find_best_layout(results, n, adj)
    crossings = count_all_crossings(coords, adj, n)
    print(f"  Layout: {best.get('label','?')}, {crossings} crossings")

    # 3. Generate TikZ
    tikz = generate_tikz(coords, adj, n, spec, aut_order, aut_structure)

    tex_path = os.path.join(build_dir, f"{gap_name}.tex")
    pdf_path = os.path.join(build_dir, f"{gap_name}.pdf")
    jpg_path = os.path.join(OUTPUT_DIR, output_file)

    with open(tex_path, "w") as f:
        f.write(tikz)

    # 4. Compile LaTeX
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", f"{gap_name}.tex"],
        capture_output=True, text=True, cwd=build_dir
    )
    if result.returncode != 0:
        print(f"  ERROR: LaTeX compilation failed!")
        print(result.stdout[-500:] if result.stdout else "")
        return False

    if not os.path.exists(pdf_path):
        print(f"  ERROR: PDF not created!")
        return False

    # 5. Convert PDF → JPG
    result = subprocess.run(
        ["sips", "-s", "format", "jpeg", "-s", "formatOptions", "95",
         "--resampleWidth", "2400", pdf_path, "--out", jpg_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  ERROR: sips conversion failed!")
        print(result.stderr)
        return False

    print(f"  → {jpg_path}")
    return True


def main():
    # Filter by name if argument given
    if len(sys.argv) > 1:
        target = sys.argv[1].lower()
        specs = [s for s in NAMED_GRAPHS
                 if target in s["display_name"].lower() or target == s["gap_name"]]
        if not specs:
            print(f"Unknown graph: {sys.argv[1]}")
            print(f"Available: {', '.join(s['display_name'] for s in NAMED_GRAPHS)}")
            sys.exit(1)
    else:
        specs = NAMED_GRAPHS

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Use a temp directory for build artifacts
    build_dir = tempfile.mkdtemp(prefix="named_graphs_")
    print(f"Build directory: {build_dir}")

    success = 0
    failed = 0
    for spec in specs:
        if build_graph(spec, build_dir):
            success += 1
        else:
            failed += 1

    # Clean up build artifacts
    shutil.rmtree(build_dir, ignore_errors=True)

    print(f"\n{'='*60}")
    print(f"  Done: {success} succeeded, {failed} failed")
    print(f"  Output: {OUTPUT_DIR}")
    print(f"{'='*60}")

    sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()
