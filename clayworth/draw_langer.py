#!/usr/bin/env python3
"""Draw the Langer graph and its Z_7 quotient graph."""

import math
import subprocess
import os

# ============================================================
# Langer graph adjacency (Zorn product formula)
# ============================================================
def to_bit_vec(n):
    return tuple((n >> i) & 1 for i in range(7))

def q62form(x):
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

q62_indices = [n for n in range(1, 128) if q62form(to_bit_vec(n)) == 0]

def langer_adj(i, j):
    if i == j: return False
    x, y = to_bit_vec(q62_indices[i]), to_bit_vec(q62_indices[j])
    alpha, gamma = x[3], y[3]
    if (alpha*gamma + x[0]*y[4] + x[1]*y[5] + x[2]*y[6]) % 2 != 0: return False
    if (alpha*gamma + x[4]*y[0] + x[5]*y[1] + x[6]*y[2]) % 2 != 0: return False
    a, b = (x[0],x[1],x[2]), (x[4],x[5],x[6])
    c, d = (y[0],y[1],y[2]), (y[4],y[5],y[6])
    bxd = ((b[1]*d[2]+b[2]*d[1])%2, (b[2]*d[0]+b[0]*d[2])%2, (b[0]*d[1]+b[1]*d[0])%2)
    axc = ((a[1]*c[2]+a[2]*c[1])%2, (a[2]*c[0]+a[0]*c[2])%2, (a[0]*c[1]+a[1]*c[0])%2)
    for k in range(3):
        if (gamma*a[k] + alpha*c[k] + bxd[k]) % 2 != 0: return False
        if (gamma*b[k] + alpha*d[k] + axc[k]) % 2 != 0: return False
    return True

adj = [[j for j in range(63) if langer_adj(i, j)] for i in range(63)]

# ============================================================
# Z_7 orbits (from companion matrix of x^3 + x + 1)
# ============================================================
orbits = [
    [0, 1, 3, 2, 5, 6, 4],
    [7, 23, 55, 47, 39, 15, 31],
    [8, 25, 56, 49, 42, 18, 32],
    [9, 24, 58, 50, 41, 16, 33],
    [10, 26, 57, 48, 40, 17, 34],
    [11, 28, 61, 52, 46, 22, 36],
    [12, 30, 62, 54, 43, 19, 35],
    [13, 27, 60, 53, 44, 21, 38],
    [14, 29, 59, 51, 45, 20, 37],
]
orbit_map = [0]*63
for oi, orb in enumerate(orbits):
    for v in orb:
        orbit_map[v] = oi

# Quotient edges
q_edges = set()
for v in range(63):
    for u in adj[v]:
        oi, oj = orbit_map[v], orbit_map[u]
        if oi < oj:
            q_edges.add((oi, oj))
q_edges = sorted(q_edges)

# ============================================================
# TikZ: Full Langer graph on 9 concentric heptagonal rings
# ============================================================

def angular_offset_search(orbits, adj, n_offsets=7):
    """Find angular offsets for each ring to minimize crossings."""
    best_offsets = [0] * 9
    # Simple greedy: try each orbit's offset
    for oi in range(1, 9):
        best_score = float('inf')
        best_off = 0
        for off in range(7):
            offsets = list(best_offsets)
            offsets[oi] = off
            # Count crossings approximately (use edge length sum as proxy)
            score = 0
            for v in range(63):
                for u in adj[v]:
                    if v < u and orbit_map[v] != orbit_map[u]:
                        ov, ou = orbit_map[v], orbit_map[u]
                        pv = orbits[ov].index(v)
                        pu = orbits[ou].index(u)
                        av = (pv + offsets[ov]) % 7
                        au = (pu + offsets[ou]) % 7
                        rv = 1.0 + ov * 0.7
                        ru = 1.0 + ou * 0.7
                        dx = rv*math.cos(2*math.pi*av/7) - ru*math.cos(2*math.pi*au/7)
                        dy = rv*math.sin(2*math.pi*av/7) - ru*math.sin(2*math.pi*au/7)
                        score += dx*dx + dy*dy
            if score < best_score:
                best_score = score
                best_off = off
        best_offsets[oi] = best_off
    return best_offsets

offsets = angular_offset_search(orbits, adj)

def vertex_pos(v, scale=1.0):
    oi = orbit_map[v]
    pos_in_orbit = orbits[oi].index(v)
    r = 1.0 + oi * 0.65
    angle = 2 * math.pi * ((pos_in_orbit + offsets[oi]) % 7) / 7 + math.pi/2
    return (r * math.cos(angle) * scale, r * math.sin(angle) * scale)

lines = []
lines.append(r"\documentclass[border=12pt]{standalone}")
lines.append(r"\usepackage{tikz}")
lines.append(r"\begin{document}")
lines.append(r"\begin{tikzpicture}[")
lines.append(r"  edge/.style={thin, black!25},")
lines.append(r"  intedge/.style={thin, black!15},")
lines.append(r"]")

# Title
lines.append(r"\node[font=\Large\bfseries] at (0, 8.2) {Langer graph (63 vertices)};")
lines.append(r"\node[font=\small, text=gray] at (0, 7.5) "
             r"{Collinearity graph of $\mathrm{GH}(2,2)$ on $Q(6,2)$, "
             r"$Z_7$ layout (9 heptagonal rings)};")

# Draw edges
drawn = set()
for v in range(63):
    for u in adj[v]:
        if v < u:
            x1, y1 = vertex_pos(v)
            x2, y2 = vertex_pos(u)
            style = "intedge" if orbit_map[v] == orbit_map[u] else "edge"
            lines.append(f"\\draw[{style}] ({x1:.3f}, {y1:.3f}) -- ({x2:.3f}, {y2:.3f});")

# Draw vertices colored by orbit
colors = ["red", "blue", "green!60!black", "orange", "purple",
          "cyan!70!black", "brown", "magenta", "olive"]
for v in range(63):
    x, y = vertex_pos(v)
    c = colors[orbit_map[v]]
    lines.append(f"\\fill[{c}] ({x:.3f}, {y:.3f}) circle (2pt);")

lines.append(r"\end{tikzpicture}")
lines.append(r"\end{document}")

with open("/tmp/langer_full.tex", "w") as f:
    f.write("\n".join(lines))

# ============================================================
# TikZ: Quotient graph (9 vertices, 19 edges)
# ============================================================

lines2 = []
lines2.append(r"\documentclass[border=12pt]{standalone}")
lines2.append(r"\usepackage{tikz}")
lines2.append(r"\begin{document}")
lines2.append(r"\begin{tikzpicture}[")
lines2.append(r"  qnode/.style={circle, draw, fill=#1!30, minimum size=18pt, "
              r"font=\footnotesize\bfseries},")
lines2.append(r"  qedge/.style={thick, black!50},")
lines2.append(r"]")

lines2.append(r"\node[font=\Large\bfseries] at (0, 4.5) {$Z_7$ quotient of the Langer graph};")
lines2.append(r"\node[font=\small, text=gray] at (0, 3.9) {9 vertices (orbits), 19 edges};")

# Place 9 vertices on a circle
qpos = {}
for i in range(9):
    angle = 2 * math.pi * i / 9 + math.pi/2
    r = 2.5
    qpos[i] = (r * math.cos(angle), r * math.sin(angle))

# Draw edges
for (a, b) in q_edges:
    x1, y1 = qpos[a]
    x2, y2 = qpos[b]
    lines2.append(f"\\draw[qedge] ({x1:.3f}, {y1:.3f}) -- ({x2:.3f}, {y2:.3f});")

# Draw vertices
for i in range(9):
    x, y = qpos[i]
    c = colors[i].split("!")[0]  # simplified color
    lines2.append(f"\\node[qnode={c}] at ({x:.3f}, {y:.3f}) {{{i}}};")

lines2.append(r"\end{tikzpicture}")
lines2.append(r"\end{document}")

with open("/tmp/langer_quotient.tex", "w") as f:
    f.write("\n".join(lines2))

# ============================================================
# Compile both
# ============================================================

for name in ["langer_full", "langer_quotient"]:
    tex = f"/tmp/{name}.tex"
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", f"{name}.tex"],
        capture_output=True, text=True, cwd="/tmp"
    )
    if result.returncode != 0:
        print(f"LaTeX error for {name}!")
        print(result.stdout[-500:])
    else:
        pdf = f"/tmp/{name}.pdf"
        jpg = f"/tmp/{name}.jpg"
        subprocess.run(
            ["sips", "-s", "format", "jpeg", "-s", "formatOptions", "95",
             "--resampleWidth", "2400", pdf, "--out", jpg],
            capture_output=True
        )
        print(f"✓ {jpg}")
