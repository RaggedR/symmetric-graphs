#!/usr/bin/env python3
"""Zhou Z2 quotient (91v, 6-regular) coloured by its cover's rings.

The cover Zhou-3 (F182A) is drawn in named_graphs/zhou3-F182A.jpg with 14
Z13-rings coloured by the house palette (blue outer -> red inner). The 2:1
covering F182A -> 91v has fibres = blocks of size 2 under Aut = PSL(2,13).
Each fibre meets two cover rings, fusing the 14 rings into 7 pairs; the
quotient therefore has 7 rings of 13 that inherit the cover's colours:
fill = outer cover ring, border = inner cover ring.

Output: individual_graphs/zhou-z2-quotient.{tex,pdf}
"""

import json
import math
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.expanduser("~/git/graphviz/src"))
from draw_graph import run_gap, find_best_layout  # noqa: E402

ROOT = os.path.expanduser("~/git/graphviz")

RING_PALETTE = [
    "0033CC", "0066FF", "0099FF", "00CCCC", "00CC66",
    "33CC00", "99CC00", "CCCC00", "FF9900", "FF6600",
    "FF3300", "CC0000",
]


# ------------------------------------------------------------
# Step 1: reproduce the cover drawing's ring assignment
# ------------------------------------------------------------

def cover_rings():
    data = run_gap("psl2_13", 4)
    n, adj = data["n"], data["adj"]
    assert n == 182
    coords, best = find_best_layout(data["results"], n, adj,
                                    layout_mode="auto", verbose=False)
    orbits = best["orbits"]  # generator-ordered rings (1-indexed vertices)
    # Rank by mean radius, exactly as generate_named_graphs.generate_tikz
    mean_r = [sum(math.hypot(*coords[v]) for v in o) / len(o) for o in orbits]
    order = sorted(range(len(orbits)), key=lambda oi: -mean_r[oi])
    ring_rank = {oi: rank for rank, oi in enumerate(order)}

    def color_idx(oi):
        return min(11, int(ring_rank[oi] * 12 / len(orbits)))

    return adj, orbits, ring_rank, color_idx


# ------------------------------------------------------------
# Step 2: fibres via GAP block system
# ------------------------------------------------------------

def compute_blocks(adj):
    with open("/tmp/zhou3_adj.g", "w") as f:
        f.write("ZHOU3_ADJ := " + json.dumps(adj) + ";\n")
    gap_script = """
LoadPackage("grape");;
Read("/tmp/zhou3_adj.g");;
gamma := Graph(Group(()), [1..182], OnPoints,
               function(x,y) return y in ZHOU3_ADJ[x]; end, true);;
A := AutomorphismGroup(gamma);;
Print("AUTORDER ", Size(A), "\\n");
blocks := Blocks(A, [1..182]);;
Print("BLOCKS ", Length(blocks), " ", Length(blocks[1]), "\\n");
out := OutputTextFile("/tmp/zhou3_blocks.json", false);;
SetPrintFormattingStatus(out, false);;
PrintTo(out, blocks);;
CloseStream(out);;
"""
    subprocess.run(["gap", "--quitonbreak", "-q"], input=gap_script,
                   capture_output=True, text=True, timeout=300)
    with open("/tmp/zhou3_blocks.json") as f:
        blocks = json.loads(f.read().replace(" ", ""))
    assert len(blocks) == 91 and all(len(b) == 2 for b in blocks), \
        f"unexpected block system: {len(blocks)} blocks"
    print(f"covering: 91 fibres of size 2 (Aut block system)")
    return blocks


# ------------------------------------------------------------
# Step 3: quotient graph + ring fusion
# ------------------------------------------------------------

def build_quotient(adj, blocks):
    fibre_of = {}  # cover vertex (1-idx) -> quotient vertex (0-idx)
    for qi, b in enumerate(blocks):
        for v in b:
            fibre_of[v] = qi
    qadj = [set() for _ in range(91)]
    for v in range(1, 183):
        for w in adj[v - 1]:
            if fibre_of[v] != fibre_of[w]:
                qadj[fibre_of[v]].add(fibre_of[w])
    degs = sorted({len(s) for s in qadj})
    assert degs == [6], f"quotient degrees {degs}"
    edges = sorted({(min(a, b), max(a, b))
                    for a in range(91) for b in qadj[a]})
    assert len(edges) == 91 * 6 // 2
    return fibre_of, edges


def fuse_rings(orbits, ring_rank, fibre_of, blocks):
    """Pair up cover rings via fibres; return quotient rings in cover order."""
    orbit_of = {}
    for oi, orb in enumerate(orbits):
        for v in orb:
            orbit_of[v] = oi
    # which pair of cover rings does each fibre meet?
    pair_of_fibre = [tuple(sorted((orbit_of[b[0]], orbit_of[b[1]]),
                                  key=lambda oi: ring_rank[oi]))
                     for b in blocks]
    pairs = sorted(set(pair_of_fibre), key=lambda p: ring_rank[p[0]])
    assert len(pairs) == 7, f"ring fusion gave {len(pairs)} pairs: {pairs}"
    assert all(p[0] != p[1] for p in pairs), "fibre inside a single ring"
    print("cover-ring pairs (outer,inner ranks):",
          [(ring_rank[a], ring_rank[b]) for a, b in pairs])

    # quotient ring = image of the OUTER cover ring, generator order inherited
    qrings = []
    for outer, inner in pairs:
        ring = [fibre_of[v] for v in orbits[outer]]
        assert len(set(ring)) == 13
        qrings.append(ring)
    return pairs, qrings


# ------------------------------------------------------------
# Step 4: layout — 7 rings of 13, offset/reflection search
# ------------------------------------------------------------

def make_coords(qrings, offsets, mirrors, radii):
    coords = [None] * 91
    for ri, ring in enumerate(qrings):
        seq = ring[::-1] if mirrors[ri] else ring
        for j, v in enumerate(seq):
            th = 2 * math.pi * (j + offsets[ri]) / 13 + math.pi / 2
            coords[v] = (radii[ri] * math.cos(th), radii[ri] * math.sin(th))
    return coords


def edge_len_sq(coords, edges):
    return sum((coords[a][0] - coords[b][0]) ** 2 +
               (coords[a][1] - coords[b][1]) ** 2 for a, b in edges)


def count_crossings(coords, edges):
    EPS = 1e-9
    cr = 0
    m = len(edges)
    for i in range(m):
        a, b = edges[i]
        ax, ay = coords[a]
        bx, by = coords[b]
        for j in range(i + 1, m):
            c, d = edges[j]
            if c in (a, b) or d in (a, b):
                continue
            cx, cy = coords[c]
            dx, dy = coords[d]
            d1 = (dx - cx) * (ay - cy) - (dy - cy) * (ax - cx)
            d2 = (dx - cx) * (by - cy) - (dy - cy) * (bx - cx)
            d3 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            d4 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax)
            if abs(d1 - d2) < EPS or abs(d3 - d4) < EPS:
                continue
            t = d1 / (d1 - d2)
            s = d3 / (d3 - d4)
            if EPS < t < 1 - EPS and EPS < s < 1 - EPS:
                cr += 1
    return cr


def search_layout(qrings, edges, n_mirror_samples=64):
    radii = [7.5 - ri * (7.5 - 2.2) / 6 for ri in range(7)]
    rng = random.Random(42)
    best = (float("inf"), None, None)
    mirror_sets = [[False] * 7]
    for _ in range(n_mirror_samples - 1):
        mirror_sets.append([False] + [rng.random() < 0.5 for _ in range(6)])
    for mirrors in mirror_sets:
        offsets = [0.0] * 7
        for ri in range(1, 7):
            best_off, best_sc = 0.0, float("inf")
            for off in [x / 2 for x in range(26)]:  # half-steps
                offsets[ri] = off
                sc = edge_len_sq(make_coords(qrings, offsets, mirrors, radii),
                                 edges)
                if sc < best_sc:
                    best_sc, best_off = sc, off
            offsets[ri] = best_off
        proxy = edge_len_sq(make_coords(qrings, offsets, mirrors, radii),
                            edges)
        if proxy < best[0]:
            best = (proxy, list(offsets), list(mirrors))
    _, offsets, mirrors = best
    bc = count_crossings(make_coords(qrings, offsets, mirrors, radii), edges)
    improved = True
    while improved:
        improved = False
        for ri in range(7):
            for off in [x / 2 for x in range(26)]:
                if off == offsets[ri]:
                    continue
                trial = list(offsets)
                trial[ri] = off
                c = count_crossings(
                    make_coords(qrings, trial, mirrors, radii), edges)
                if c < bc:
                    bc, offsets = c, trial
                    improved = True
    print(f"crossings: {bc}")
    return make_coords(qrings, offsets, mirrors, radii), bc


# ------------------------------------------------------------
# Step 5: TikZ
# ------------------------------------------------------------

def generate_tikz(coords, edges, qrings, pairs, ring_rank, color_idx,
                  crossings):
    n_orbits = 14  # cover orbits, for colour scaling
    fill_of, border_of, qring_of = {}, {}, {}
    for ri, (ring, (outer, inner)) in enumerate(zip(qrings, pairs)):
        for v in ring:
            fill_of[v] = color_idx(outer)
            border_of[v] = color_idx(inner)
            qring_of[v] = ri

    lines = [
        r"\documentclass[border=12pt]{standalone}",
        r"\usepackage{amsmath, amssymb}",
        r"\usepackage{tikz}",
    ]
    for idx, hex_c in enumerate(RING_PALETTE):
        lines.append(f"\\definecolor{{ring{idx}}}{{HTML}}{{{hex_c}}}")
    lines += [r"\begin{document}", r"\begin{tikzpicture}"]

    y = max(c[1] for c in coords) + 1.7
    lines.append(f"\\node[font=\\Large\\bfseries] at (0, {y:.1f}) "
                 r"{Zhou $Z_2$ quotient};")
    lines.append(
        f"\\node[font=\\small, text=gray] at (0, {y - 0.7:.1f}) "
        r"{$\mathrm{Sab}(\mathrm{PSL}(2,13), D_{12})$, 91 vertices, "
        r"6-regular --- coloured by the rings of its cover F182A};")
    lines.append(f"\\node[font=\\footnotesize, text=gray!70] at "
                 f"(0, {y - 1.25:.1f}) {{fill = outer fibre ring, "
                 f"border = inner fibre ring; {crossings} crossings}};")

    # edges: intra-ring coloured by ring fill, inter-ring light grey
    for a, b in edges:
        x1, y1 = coords[a]
        x2, y2 = coords[b]
        if qring_of[a] == qring_of[b]:
            style = f"line width=0.4pt, ring{fill_of[a]}, opacity=0.35"
        else:
            style = "line width=0.3pt, gray!40, opacity=0.45"
        lines.append(f"\\draw[{style}] ({x1:.3f}, {y1:.3f}) -- "
                     f"({x2:.3f}, {y2:.3f});")

    for v in range(91):
        x, yv = coords[v]
        lines.append(
            f"\\fill[ring{fill_of[v]}] ({x:.3f}, {yv:.3f}) circle (2.8pt);"
            f"\\draw[ring{border_of[v]}, line width=1.0pt] "
            f"({x:.3f}, {yv:.3f}) circle (2.8pt);")

    lines += [r"\end{tikzpicture}", r"\end{document}"]
    return "\n".join(lines)


if __name__ == "__main__":
    adj, orbits, ring_rank, color_idx = cover_rings()
    print(f"cover: {len(orbits)} rings of {len(orbits[0])}")
    blocks = compute_blocks(adj)
    fibre_of, edges = build_quotient(adj, blocks)
    pairs, qrings = fuse_rings(orbits, ring_rank, fibre_of, blocks)
    coords, crossings = search_layout(qrings, edges)
    tex = generate_tikz(coords, edges, qrings, pairs, ring_rank, color_idx,
                        crossings)
    out = os.path.join(ROOT, "individual_graphs", "zhou-z2-quotient.tex")
    with open(out, "w") as f:
        f.write(tex)
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                        "zhou-z2-quotient.tex"],
                       capture_output=True, text=True,
                       cwd=os.path.dirname(out))
    if r.returncode != 0:
        print(r.stdout[-600:])
        raise SystemExit("pdflatex failed")
    for ext in (".aux", ".log"):
        p = out[:-4] + ext
        if os.path.exists(p):
            os.unlink(p)
    print(f"written {out} (+pdf)")
