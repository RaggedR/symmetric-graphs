#!/usr/bin/env python3
"""
Draw all quotients of cubic arc-transitive graphs.

For each graph Γ, produces a multi-page PDF with 4 quotient "cards" per page.
Each card shows Γ (vertices colored by H-orbits) beside Γ/H (the quotient graph).

Usage:
    python3 draw_quotients.py petersen       # Single graph → quotients/petersen.pdf
    python3 draw_quotients.py               # All 13 graphs
"""

import subprocess
import sys
import os
import re
import math

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
GRAPHVIZ_DIR = os.path.dirname(SRC_DIR)
QUOTIENTS_DIR = os.path.join(GRAPHVIZ_DIR, "quotients")

from draw_graph import COLORS, run_gap, find_best_layout

ALL_GRAPHS = [
    "k4", "k33", "cube", "petersen", "heawood",
    "moebiuskantor", "pappus", "dodecahedron", "desargues",
    "nauru", "f26a", "coxeter", "tuttecoxeter",
]

DISPLAY_NAMES = {
    "K4": "$K_4$",
    "K33": "$K_{3,3}$",
    "Cube": "Cube ($Q_3$)",
    "Petersen": "Petersen",
    "Heawood": "Heawood",
    "Moebius-Kantor": "M\\\"obius--Kantor",
    "Pappus": "Pappus",
    "Dodecahedron": "Dodecahedron",
    "Desargues": "Desargues",
    "Nauru": "Nauru",
    "F26A": "$F_{26A}$",
    "Coxeter": "Coxeter",
    "Tutte-Coxeter": "Tutte--Coxeter",
}


# ============================================================
# GAP interface
# ============================================================

def run_gap_enumerate(graph_name):
    """Run enumerate-quotients.g in JSON mode for one graph."""
    gap_script = os.path.join(SRC_DIR, "enumerate-quotients.g")
    gap_input = (
        f'json_mode := true;\n'
        f'target_graphs := ["{graph_name}"];\n'
        f'Read("{gap_script}");\n'
    )
    result = subprocess.run(
        ["gap", "-q"], input=gap_input,
        capture_output=True, text=True, timeout=600,
        cwd=SRC_DIR,  # enumerate-quotients.g Read()s symmetric-layout.g relative to cwd
    )
    output = result.stdout + result.stderr
    match = re.search(
        r'=== JSON OUTPUT ===\n(.*?)\n=== END JSON ===', output, re.DOTALL
    )
    if not match:
        print(f"GAP output (last 2000 chars):\n{output[-2000:]}")
        raise RuntimeError(f"No JSON in GAP output for {graph_name}")
    json_text = match.group(1).replace("true", "True").replace("false", "False")
    return eval(json_text)


# ============================================================
# Layout: original graph (grouped circular by orbits)
# ============================================================

def layout_grouped_circular(n, orbits, R=2.2):
    """
    Circular layout with vertices clustered by orbit.
    Gaps between orbit groups make the partition visually clear.
    """
    n_orbits = len(orbits)
    gap = max(0.5, 1.5 * min(1.0, 6.0 / max(1, n_orbits)))
    total_slots = n + n_orbits * gap

    coords = {}
    pos = 0.0
    for orbit in orbits:
        for v in sorted(orbit):
            theta = math.radians(90 + pos * (360.0 / total_slots))
            coords[v] = (R * math.cos(theta), R * math.sin(theta))
            pos += 1
        pos += gap
    return coords


# ============================================================
# Layout: quotient graph
# ============================================================

def layout_quotient_graph(n_orbits, quot_adj, qtype, orbit_sizes):
    """
    Position quotient graph vertices based on identified type.
    quot_adj is 1-indexed (from GAP). Returns 0-indexed dict.
    """
    if n_orbits == 1:
        return {0: (0, 0)}
    if n_orbits == 2:
        return {0: (-0.7, 0), 1: (0.7, 0)}

    R = 1.0

    # Convert to 0-indexed adjacency
    adj0 = [[j - 1 for j in quot_adj[i]] for i in range(n_orbits)]
    degrees = [len(a) for a in adj0]
    endpoints = [i for i in range(n_orbits) if degrees[i] == 1]

    # Path
    if len(endpoints) == 2 and all(d <= 2 for d in degrees):
        path = _traverse_path(n_orbits, adj0, endpoints[0])
        if path and len(path) == n_orbits:
            coords = {}
            for idx, node in enumerate(path):
                x = -1.2 + idx * (2.4 / max(1, n_orbits - 1))
                coords[node] = (x, 0)
            return coords

    # Star (K1,n)
    if "star" in qtype.lower() or qtype.startswith("K1,"):
        center = max(range(n_orbits), key=lambda i: degrees[i])
        others = [i for i in range(n_orbits) if i != center]
        coords = {center: (0, 0)}
        for idx, node in enumerate(others):
            theta = math.radians(90 + idx * (360.0 / len(others)))
            coords[node] = (R * math.cos(theta), R * math.sin(theta))
        return coords

    # Default: circle
    coords = {}
    for i in range(n_orbits):
        theta = math.radians(90 + i * (360.0 / n_orbits))
        coords[i] = (R * math.cos(theta), R * math.sin(theta))
    return coords


def _traverse_path(n, adj0, start):
    """Walk a path graph from start vertex."""
    path = [start]
    visited = {start}
    while len(path) < n:
        advanced = False
        for nb in adj0[path[-1]]:
            if nb not in visited:
                path.append(nb)
                visited.add(nb)
                advanced = True
                break
        if not advanced:
            break
    return path


# ============================================================
# LaTeX helpers
# ============================================================

def latex_group(name):
    """Convert GAP StructureDescription to LaTeX-safe text."""
    s = name
    s = s.replace(" x ", r" $\times$ ")
    s = s.replace(" : ", r" $\rtimes$ ")
    return s


def loop_direction(x, y):
    """Choose self-loop direction based on vertex position."""
    if abs(x) < 0.01 and abs(y) < 0.01:
        return "above"
    if abs(y) >= abs(x):
        return "above" if y >= 0 else "below"
    return "right" if x >= 0 else "left"


# ============================================================
# TikZ: one quotient card
# ============================================================

def tikz_card(card, cx, cy, card_idx):
    """Generate TikZ commands for one quotient card at position (cx, cy)."""
    lines = []
    n = card['n']
    adj = card['adj']
    orbits = card['orbits']
    n_orbits = card['n_orbits']
    orig_coords = card['orig_coords']
    quot_coords = card['quot_coords']
    quot_adj = card['quotient_adj']
    edge_counts = card['edge_counts']
    internal_counts = card['internal_counts']
    pfx = f"c{card_idx}"

    # Vertex-to-orbit map
    v2o = {}
    for oi, orb in enumerate(orbits):
        for v in orb:
            v2o[v] = oi

    # --- Title ---
    h_tex = latex_group(card['subgroup_structure'])
    orbit_size = card['orbit_sizes'][0]  # all equal
    lines.append(
        f"\\node[cardtitle] at ({cx:.2f},{cy+5.0:.2f}) "
        f"{{{card['display_name']} / $H \\cong$ {h_tex}}};"
    )
    lines.append(
        f"\\node[cardsub] at ({cx:.2f},{cy+4.5:.2f}) "
        f"{{$|H|={card['subgroup_order']}$, "
        f"{n_orbits} $\\times$ {orbit_size}, "
        f"$\\Gamma/H \\cong$ {card['quotient_type']}}};"
    )

    # --- Original graph (left half) ---
    ox, oy = cx - 1.8, cy + 0.2
    s = 0.55

    # Edges (drawn first = behind vertices)
    drawn = set()
    for v in range(1, n + 1):
        for w in adj[v - 1]:
            if w > v:
                drawn.add((v, w))
    for v, w in sorted(drawn):
        vx = ox + orig_coords[v][0] * s
        vy = oy + orig_coords[v][1] * s
        wx = ox + orig_coords[w][0] * s
        wy = oy + orig_coords[w][1] * s
        if v2o[v] == v2o[w]:
            style = f"c{v2o[v] % len(COLORS)}!50, thin"
        else:
            style = "gray!30, thin"
        lines.append(f"\\draw[{style}] ({vx:.3f},{vy:.3f}) -- ({wx:.3f},{wy:.3f});")

    # Vertices (on top of edges)
    for v in range(1, n + 1):
        x, y = orig_coords[v]
        ci = v2o[v] % len(COLORS)
        vx = ox + x * s
        vy = oy + y * s
        if n <= 14:
            lines.append(
                f"\\node[ov, fill=c{ci}!30] at ({vx:.3f},{vy:.3f}) "
                f"{{\\scalebox{{0.45}}{{{v}}}}};"
            )
        else:
            lines.append(
                f"\\node[ovs, fill=c{ci}!40] at ({vx:.3f},{vy:.3f}) {{}};"
            )

    lines.append(
        f"\\node[glabel] at ({ox:.2f},{oy - 1.6:.2f}) {{$\\Gamma$}};"
    )

    # --- Quotient graph (right half) ---
    qx_c, qy_c = cx + 2.3, cy + 0.2
    qs = 0.8

    # Define named nodes (uniform size — all orbits equal)
    for oi in range(n_orbits):
        ci = oi % len(COLORS)
        px = qx_c + quot_coords[oi][0] * qs
        py = qy_c + quot_coords[oi][1] * qs
        lines.append(
            f"\\node[qv, fill=c{ci}!30] ({pfx}o{oi}) "
            f"at ({px:.3f},{py:.3f}) "
            f"{{\\scalebox{{0.5}}{{$O_{{{oi+1}}}$}}}};"
        )

    # Edges between orbit vertices
    drawn_q = set()
    for i in range(n_orbits):
        for j1 in quot_adj[i]:
            j = j1 - 1
            edge = (min(i, j), max(i, j))
            if edge not in drawn_q:
                drawn_q.add(edge)
                mult = edge_counts[i][j]
                label = ""
                if mult > 1:
                    label = (
                        f"node[font=\\tiny, auto, inner sep=1pt] "
                        f"{{{mult}}} "
                    )
                lines.append(
                    f"\\draw[gray!60, thick] ({pfx}o{edge[0]}) -- "
                    f"{label}({pfx}o{edge[1]});"
                )

    lines.append(
        f"\\node[glabel] at ({qx_c:.2f},{qy_c - 1.6:.2f}) {{$\\Gamma/H$}};"
    )

    # --- Legend (orbit membership) ---
    ly_start = cy - 3.8
    lx = cx - 4.0
    max_show = 8
    ncols = 2 if n_orbits > 4 else 1
    col_w = 4.2
    shown = min(n_orbits, max_show)

    for oi in range(shown):
        ci = oi % len(COLORS)
        verts = sorted(orbits[oi])
        vs = ",".join(str(v) for v in verts)
        if len(vs) > 28:
            vs = ",".join(str(v) for v in verts[:7]) + ",\\ldots"
        col = oi % ncols
        row = oi // ncols
        x = lx + col * col_w
        y = ly_start - row * 0.35
        lines.append(
            f"\\node[anchor=west, font=\\tiny] "
            f"at ({x:.2f},{y:.2f}) "
            f"{{\\textcolor{{c{ci}}}{{$\\blacksquare$}} "
            f"$O_{{{oi+1}}}$=\\{{{vs}\\}}}};"
        )

    if n_orbits > max_show:
        x = lx
        y = ly_start - (shown // ncols) * 0.35
        lines.append(
            f"\\node[anchor=west, font=\\tiny, text=gray] "
            f"at ({x:.2f},{y:.2f}) "
            f"{{\\ldots\\ and {n_orbits - max_show} more orbits}};"
        )

    return "\n".join(lines)


# ============================================================
# TikZ: complete multi-page document
# ============================================================

def generate_document(graph_name, display_name, cards):
    """Generate complete .tex document with 4 cards per page."""
    lines = []
    lines.append(r"\documentclass[a4paper]{article}")
    lines.append(r"\usepackage[margin=1cm]{geometry}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{graphicx}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\usetikzlibrary{arrows.meta}")
    lines.append(r"\pagestyle{empty}")
    lines.append("")

    for i in range(min(15, len(COLORS))):
        lines.append(f"\\definecolor{{c{i}}}{{HTML}}{{{COLORS[i]}}}")
    lines.append("")

    lines.append(r"\begin{document}")

    cards_per_page = 4
    n_pages = math.ceil(len(cards) / cards_per_page)

    # 2x2 grid positions (card centers)
    positions = [
        (-4.5,  5.5),   # top-left
        ( 4.5,  5.5),   # top-right
        (-4.5, -5.5),   # bottom-left
        ( 4.5, -5.5),   # bottom-right
    ]

    for page in range(n_pages):
        start = page * cards_per_page
        end = min(start + cards_per_page, len(cards))
        page_cards = cards[start:end]

        lines.append("")
        lines.append(r"\begin{center}")
        lines.append(r"\begin{tikzpicture}[")
        lines.append(r"  ov/.style={circle, draw, semithick, "
                     r"minimum size=3.5mm, inner sep=0pt},")
        lines.append(r"  ovs/.style={circle, draw, thin, "
                     r"minimum size=1.8mm, inner sep=0pt},")
        lines.append(r"  qv/.style={circle, draw, thick, "
                     r"minimum size=5mm, inner sep=0pt},")
        lines.append(r"  cardtitle/.style={font=\small\bfseries},")
        lines.append(r"  cardsub/.style={font=\scriptsize, text=gray},")
        lines.append(r"  glabel/.style={font=\footnotesize, text=gray},")
        lines.append(r"]")

        # Page header
        if page == 0:
            lines.append(
                f"\\node[font=\\Large\\bfseries] at (0, 11.5) "
                f"{{Quotients of {display_name}}};"
            )
            aut_info = ""
            if cards:
                aut_info = f", $|\\mathrm{{Aut}}|={cards[0].get('aut_order', '?')}$"
            lines.append(
                f"\\node[font=\\small, text=gray] at (0, 10.8) "
                f"{{{len(cards)} distinct quotients{aut_info}}};"
            )
        else:
            lines.append(
                f"\\node[font=\\normalsize, text=gray] at (0, 11.5) "
                f"{{{display_name} --- page {page+1}/{n_pages}}};"
            )

        for i, card in enumerate(page_cards):
            cx, cy = positions[i]
            lines.append(
                f"\n% Card {start+i+1}: "
                f"H={card['subgroup_structure']}"
            )
            lines.append(tikz_card(card, cx, cy, start + i))

        lines.append(r"\end{tikzpicture}")
        lines.append(r"\end{center}")

        if page < n_pages - 1:
            lines.append(r"\clearpage")

    lines.append("")
    lines.append(r"\end{document}")
    return "\n".join(lines)


# ============================================================
# PDF compilation
# ============================================================

def compile_pdf(tex_path):
    """Compile .tex → .pdf via pdflatex."""
    d = os.path.dirname(tex_path)
    f = os.path.basename(tex_path)
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", f],
        capture_output=True, text=True, cwd=d,
    )
    if result.returncode != 0:
        print(f"  LaTeX FAILED for {f}")
        for line in result.stdout.split("\n"):
            if line.startswith("!") or "Error" in line:
                print(f"    {line}")
        return False
    # cleanup
    stem = f.replace(".tex", "")
    for ext in ["aux", "log"]:
        try:
            os.remove(os.path.join(d, f"{stem}.{ext}"))
        except FileNotFoundError:
            pass
    return True


# ============================================================
# Main pipeline
# ============================================================

def collect_graph_cards(graph_name):
    """Enumerate quotients and compute layouts for one graph."""
    print(f"\n  {graph_name}...")

    data = run_gap_enumerate(graph_name)

    display_raw = data['graph']
    display_name = DISPLAY_NAMES.get(display_raw, display_raw)
    n = data['n']
    adj = data['adj']
    results = data['results']
    aut_order = data.get('aut_order', '?')
    aut_structure = data.get('aut_structure', '?')

    print(f"    n={n}, |Aut|={aut_order}, {len(results)} quotients")

    if not results:
        print(f"    (primitive)")
        return display_name, []

    # Compute the base layout for the original graph (reused across all cards)
    layout_data = run_gap(graph_name, max_deg=10)
    layout_results = layout_data['results']
    if layout_results:
        base_coords, _ = find_best_layout(layout_results, n, adj, verbose=False)
    else:
        base_coords = layout_grouped_circular(n, [list(range(1, n + 1))])

    cards = []
    for qi, r in enumerate(results):
        orbits = r['orbits']
        n_orbits = r['n_orbits']
        quot_coords = layout_quotient_graph(
            n_orbits, r['quotient_adj'], r['quotient_type'], r['orbit_sizes']
        )
        cards.append({
            'n': n,
            'adj': adj,
            'display_name': display_name,
            'orbits': orbits,
            'orbit_sizes': r['orbit_sizes'],
            'n_orbits': n_orbits,
            'subgroup_structure': r['subgroup_structure'],
            'subgroup_order': r['subgroup_order'],
            'quotient_type': r['quotient_type'],
            'quotient_adj': r['quotient_adj'],
            'edge_counts': r['edge_counts'],
            'internal_counts': r['internal_counts'],
            'orig_coords': base_coords,
            'quot_coords': quot_coords,
            'aut_order': aut_order,
        })

    return display_name, cards


def generate_combined_document(all_sections):
    """Generate one document with all graphs' quotients."""
    lines = []
    lines.append(r"\documentclass[a4paper]{article}")
    lines.append(r"\usepackage[margin=1cm]{geometry}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{graphicx}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\usetikzlibrary{arrows.meta}")
    lines.append(r"\pagestyle{empty}")
    lines.append("")

    for i in range(min(15, len(COLORS))):
        lines.append(f"\\definecolor{{c{i}}}{{HTML}}{{{COLORS[i]}}}")
    lines.append("")

    lines.append(r"\begin{document}")

    # Global title page
    lines.append(r"\begin{center}")
    lines.append(r"\vspace*{4cm}")
    lines.append(r"{\Huge\bfseries Quotients of Cubic}")
    lines.append(r"\vspace{0.5cm}")
    lines.append(r"")
    lines.append(r"{\Huge\bfseries Arc-Transitive Graphs}")
    lines.append(r"\vspace{1.5cm}")
    lines.append(r"")
    lines.append(r"{\large $n \leq 30$ (Foster census)}")
    lines.append(r"\vspace{1cm}")
    lines.append(r"")
    total_q = sum(len(cards) for _, _, cards in all_sections)
    n_graphs = sum(1 for _, _, cards in all_sections if cards)
    lines.append(
        f"{{\\normalsize {total_q} block-system quotients "
        f"across {n_graphs} imprimitive graphs}}"
    )
    lines.append(r"\end{center}")
    lines.append(r"\clearpage")

    # Card positions for 2x2 grid
    positions = [
        (-4.5,  5.5),
        ( 4.5,  5.5),
        (-4.5, -5.5),
        ( 4.5, -5.5),
    ]

    global_card_idx = 0

    for graph_name, display_name, cards in all_sections:
        if not cards:
            continue

        cards_per_page = 4
        n_pages = math.ceil(len(cards) / cards_per_page)

        for page in range(n_pages):
            start = page * cards_per_page
            end = min(start + cards_per_page, len(cards))
            page_cards = cards[start:end]

            lines.append("")
            lines.append(r"\begin{center}")
            lines.append(r"\begin{tikzpicture}[")
            lines.append(r"  ov/.style={circle, draw, semithick, "
                         r"minimum size=3.5mm, inner sep=0pt},")
            lines.append(r"  ovs/.style={circle, draw, thin, "
                         r"minimum size=1.8mm, inner sep=0pt},")
            lines.append(r"  qv/.style={circle, draw, thick, "
                         r"minimum size=5mm, inner sep=0pt},")
            lines.append(r"  cardtitle/.style={font=\small\bfseries},")
            lines.append(r"  cardsub/.style={font=\scriptsize, text=gray},")
            lines.append(r"  glabel/.style={font=\footnotesize, text=gray},")
            lines.append(r"]")

            # Page header
            if page == 0:
                aut_info = ""
                if cards:
                    aut_info = (
                        f", $|\\mathrm{{Aut}}|={cards[0].get('aut_order', '?')}$"
                    )
                lines.append(
                    f"\\node[font=\\Large\\bfseries] at (0, 11.5) "
                    f"{{Quotients of {display_name}}};"
                )
                lines.append(
                    f"\\node[font=\\small, text=gray] at (0, 10.8) "
                    f"{{{len(cards)} quotient{'s' if len(cards) != 1 else ''}"
                    f"{aut_info}}};"
                )
            else:
                lines.append(
                    f"\\node[font=\\normalsize, text=gray] at (0, 11.5) "
                    f"{{{display_name} --- page {page+1}/{n_pages}}};"
                )

            for i, card in enumerate(page_cards):
                cx, cy = positions[i]
                lines.append(
                    f"\n% Card: {card['display_name']} / "
                    f"H={card['subgroup_structure']}"
                )
                lines.append(tikz_card(card, cx, cy, global_card_idx))
                global_card_idx += 1

            lines.append(r"\end{tikzpicture}")
            lines.append(r"\end{center}")
            lines.append(r"\clearpage")

    lines.append(r"\end{document}")
    return "\n".join(lines)


def main():
    print("Collecting quotients for all graphs...")

    all_sections = []
    for g in ALL_GRAPHS:
        display_name, cards = collect_graph_cards(g)
        all_sections.append((g, display_name, cards))

    # Generate combined PDF
    tex_content = generate_combined_document(all_sections)

    os.makedirs(QUOTIENTS_DIR, exist_ok=True)
    tex_file = os.path.join(QUOTIENTS_DIR, "all_quotients.tex")
    pdf_file = os.path.join(QUOTIENTS_DIR, "all_quotients.pdf")

    with open(tex_file, 'w') as f:
        f.write(tex_content)
    print(f"\nWrote {tex_file}")

    if compile_pdf(tex_file):
        print(f"-> {pdf_file}")
        subprocess.run(["open", pdf_file])
    else:
        print("Compilation failed")


if __name__ == "__main__":
    main()
