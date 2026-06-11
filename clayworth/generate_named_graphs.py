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
        "gap_name": "k4",
        "display_name": "K4",
        "foster_id": "F004A",
        "girth": 3,
        "max_deg": 10,
        "layout": "cyclic",
        "output": "k4-F004A.jpg",
    },
    {
        "gap_name": "k33",
        "display_name": "K(3,3)",
        "foster_id": "F006A",
        "girth": 4,
        "max_deg": 10,
        "layout": "cyclic",
        "output": "k33-F006A.jpg",
    },
    {
        "gap_name": "pappus",
        "display_name": "Pappus",
        "foster_id": "F018A",
        "girth": 6,
        "max_deg": 10,
        "output": "pappus-F018A.jpg",
    },
    {
        "gap_name": "desargues",
        "display_name": "Desargues",
        "foster_id": "F020A",
        "girth": 6,
        "max_deg": 10,
        "output": "desargues-F020A.jpg",
    },
    {
        "gap_name": "nauru",
        "display_name": "Nauru",
        "foster_id": "F024A",
        "girth": 6,
        "max_deg": 10,
        "output": "nauru-F024A.jpg",
    },
    {
        "gap_name": "f26a",
        "display_name": "F26A",
        "foster_id": "F026A",
        "girth": 6,
        "max_deg": 10,
        "layout": "cyclic",
        "output": "f26a-F026A.jpg",
    },
    {
        "gap_name": "coxeter",
        "display_name": "Coxeter",
        "foster_id": "F028A",
        "girth": 7,
        "max_deg": 10,
        "output": "coxeter-F028A.jpg",
    },
    {
        "gap_name": "tuttecoxeter",
        "display_name": "Tutte-Coxeter",
        "foster_id": "F030A",
        "girth": 8,
        "max_deg": 10,
        "output": "tuttecoxeter-F030A.jpg",
    },
    {
        "gap_name": "psl2_5_c3",
        "display_name": "PSL(2,5)/C3",
        "foster_id": None,
        "girth": None,
        "max_deg": 10,
        "output": "psl2_5_c3.jpg",
    },
    {
        "gap_name": "cr19",
        "display_name": "CR(19)",
        "foster_id": None,
        "girth": None,
        "max_deg": 4,
        "layout": "cyclic",
        "output": "cr19.jpg",
    },
    {
        "gap_name": "cr31",
        "display_name": "CR(31)",
        "foster_id": None,
        "girth": None,
        "max_deg": 4,
        "output": "cr31.jpg",
    },
    {
        "gap_name": "pg23",
        "display_name": "PG(2,3)",
        "foster_id": None,
        "girth": None,
        "max_deg": 10,
        "layout": "cyclic",
        "output": "pg23.jpg",
    },
    {
        "gap_name": "meinhold72",
        "display_name": "Meinhold-72",
        "foster_id": None,
        "girth": None,
        "max_deg": 30,
        "quotient": "Cyclic Z7",  # 12 heptagonal rings (Robin, 2026-06-10)
        "layout": "rings",        # concentric rings, not one circle
        "output": "meinhold72.jpg",
    },
    {
        "gap_name": "petersen_voltage",
        "display_name": "Petersen (voltage)",
        "foster_id": None,
        "girth": 5,
        "max_deg": 10,
        "output": "petersen_voltage.jpg",
    },
    {
        "gap_name": "heawood_voltage",
        "display_name": "Heawood (voltage)",
        "foster_id": None,
        "girth": 6,
        "max_deg": 10,
        "output": "heawood_voltage.jpg",
    },
    {
        "gap_name": "mk_voltage",
        "display_name": "Möbius-Kantor (voltage)",
        "foster_id": None,
        "girth": 6,
        "max_deg": 10,
        "output": "mk_voltage.jpg",
    },
    {
        "gap_name": "cube",
        "display_name": "Cube",
        "foster_id": "F008A",
        "girth": 4,
        "max_deg": 10,
        "output": "cube-F008A.jpg",
    },
    {
        "gap_name": "petersen",
        "display_name": "Petersen",
        "foster_id": "F010A",
        "girth": 5,
        "max_deg": 10,
        "output": "petersen-F010A.jpg",
    },
    {
        "gap_name": "dodecahedron",
        "display_name": "Dodecahedron",
        "foster_id": "F020B",
        "girth": 5,
        "max_deg": 10,
        "output": "dodecahedron-F020B.jpg",
    },
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
        "layout": "cyclic",
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
        "gap_name": "biggssmith",
        "display_name": "Biggs-Smith",
        "foster_id": "F102A",
        "girth": 9,
        "max_deg": 4,
        "output": "biggssmith-F102A.jpg",
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


# Ring palette (same as surfaces/draw_triality_hexagon.py and
# clayworth/draw_langer.py): blue outermost ring -> red innermost.
RING_PALETTE = [
    "0033CC", "0066FF", "0099FF", "00CCCC", "00CC66",
    "33CC00", "99CC00", "CCCC00", "FF9900", "FF6600",
    "FF3300", "CC0000",
]


def generate_tikz(coords, adj, n, spec, aut_order, aut_structure, orbits):
    """Generate TikZ in the dmt style: ring-colored vertices, gap-colored edges."""
    # Compute bounding box
    xs = [x for x, y in coords.values()]
    ys = [y for x, y in coords.values()]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    # Rank orbits by mean radius: rank 0 = outermost ring.
    n_orbits = len(orbits)
    mean_r = []
    for orbit in orbits:
        mean_r.append(sum(math.hypot(*coords[v]) for v in orbit) / len(orbit))
    order = sorted(range(n_orbits), key=lambda oi: -mean_r[oi])
    ring_rank = {oi: rank for rank, oi in enumerate(order)}
    v_to_orbit = {}
    for oi, orbit in enumerate(orbits):
        for v in orbit:
            v_to_orbit[v] = oi

    def color_idx(oi):
        return min(11, int(ring_rank[oi] * 12 / n_orbits))

    # Node size scaled by graph size
    dot_size = max(2.0, 3.5 - n / 60.0)

    # Title
    name = spec["display_name"]
    foster = spec["foster_id"]
    title = f"{name} ({foster})" if foster else name

    # Subtitle
    aut_latex = latex_safe_group(aut_structure)
    subtitle = (
        f"$\\mathrm{{Aut}} \\cong {aut_latex}$, "
        f"$|\\mathrm{{Aut}}| = {aut_order}$, "
        f"{n} vertices"
    )
    if spec.get("girth"):
        subtitle += f", girth {spec['girth']}"

    title_y = max_y + 1.2
    sub_y = title_y - 0.65

    lines = []
    lines.append(r"\documentclass[border=12pt]{standalone}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{tikz}")
    for idx, hex_c in enumerate(RING_PALETTE):
        lines.append(f"\\definecolor{{ring{idx}}}{{HTML}}{{{hex_c}}}")
    lines.append(r"\definecolor{gap0}{HTML}{E53935}")   # same ring: red
    lines.append(r"\definecolor{gap1}{HTML}{1E88E5}")   # nearby rings: blue
    lines.append(r"\definecolor{gapfar}{HTML}{CCCCCC}") # far rings: light grey
    lines.append("")
    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}[")
    lines.append(r"  title/.style={font=\Large\bfseries},")
    lines.append(r"  subtitle/.style={font=\small, text=gray},")
    lines.append(r"]")
    lines.append("")

    # Title and subtitle
    lines.append(f"\\node[title] at (0, {title_y:.2f}) {{{title}}};")
    lines.append(f"\\node[subtitle] at (0, {sub_y:.2f}) {{{subtitle}}};")
    lines.append("")

    # On a cyclic (single-circle) layout, ring ranks are arbitrary:
    # treat every inter-orbit edge as gap 1 (blue) there.
    concentric = (max(mean_r) - min(mean_r)) > 0.1

    # Edges grouped by ring gap (far edges first, close edges on top)
    edges_by_gap = {}
    for v in range(1, n + 1):
        for w in adj[v - 1]:
            if v < w:
                if concentric:
                    gap = abs(ring_rank[v_to_orbit[v]] - ring_rank[v_to_orbit[w]])
                else:
                    gap = 0 if v_to_orbit[v] == v_to_orbit[w] else 1
                edges_by_gap.setdefault(gap, []).append((v, w))

    # Small graphs get proportionally thicker, more opaque edges
    wf = max(1.0, 2.2 - n / 40.0)
    op = 0.15 if wf > 1.0 else 0.0
    for gap in sorted(edges_by_gap.keys(), reverse=True):
        if gap == 0:
            style = f"line width={0.5*wf:.2f}pt, gap0, opacity={0.6+op:.2f}"
        elif gap <= 3:
            style = f"line width={0.4*wf:.2f}pt, gap1, opacity={0.35+op:.2f}"
        else:
            style = f"line width={0.35*wf:.2f}pt, gapfar, opacity={0.45+op:.2f}"
        for v, w in edges_by_gap[gap]:
            x1, y1 = coords[v]
            x2, y2 = coords[w]
            lines.append(f"\\draw[{style}] ({x1:.3f}, {y1:.3f}) -- ({x2:.3f}, {y2:.3f});")
    lines.append("")

    # Vertices (drawn on top), colored by ring
    for v in range(1, n + 1):
        x, y = coords[v]
        c = f"ring{color_idx(v_to_orbit[v])}"
        lines.append(f"\\fill[{c}] ({x:.3f}, {y:.3f}) circle ({dot_size:.1f}pt);")

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

    # Optional: pin the layout to a specific quotient by label
    # (e.g. "Cyclic Z7"), overriding the crossing-minimizing choice.
    want = spec.get("quotient")
    if want:
        filtered = [r for r in results if r.get("label") == want]
        if not filtered:
            print(f"  ERROR: quotient {want!r} not found "
                  f"(available: {[r.get('label') for r in results]})")
            return False
        results = filtered

    # 2. Compute layout
    coords, best = find_best_layout(results, n, adj,
                                    layout_mode=spec.get("layout", "auto"))
    crossings = count_all_crossings(coords, adj, n)
    print(f"  Layout: {best.get('label','?')}, {crossings} crossings")

    # 3. Generate TikZ
    tikz = generate_tikz(coords, adj, n, spec, aut_order, aut_structure,
                         best['orbits'])

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
