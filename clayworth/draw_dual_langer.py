#!/usr/bin/env python3
"""Draw the dual Langer graph and its Z_7 quotient graph.

The dual Langer graph is the line-collinearity graph of GH(2,2):
vertices = 63 lines of GH(2,2), two lines adjacent iff they share a point
(i.e. they are "concurrent" = collinear in the dual hexagon).

Equivalently, the vertices are the 63 triangles of the Langer graph,
with two triangles adjacent iff they share a vertex.
"""

import math
import subprocess
import os

# ============================================================
# Step 1: Build the Langer graph (Zorn product formula)
# ============================================================
def to_bit_vec(n):
    return tuple((n >> i) & 1 for i in range(7))

def q62form(x):
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

q62_indices = [n for n in range(1, 128) if q62form(to_bit_vec(n)) == 0]
assert len(q62_indices) == 63, f"Expected 63 Q(6,2) points, got {len(q62_indices)}"

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

langer_nbrs = [[j for j in range(63) if langer_adj(i, j)] for i in range(63)]

# Verify 6-regular
for i in range(63):
    assert len(langer_nbrs[i]) == 6, f"Vertex {i} has degree {len(langer_nbrs[i])}"

# ============================================================
# Step 2: Extract the 63 lines of GH(2,2) = triangles of Langer
# ============================================================
triangles = []
for v in range(63):
    nbrs = langer_nbrs[v]
    for i in range(len(nbrs)):
        for j in range(i+1, len(nbrs)):
            a, b = nbrs[i], nbrs[j]
            if b in langer_nbrs[a]:
                tri = tuple(sorted([v, a, b]))
                if tri not in triangles:
                    triangles.append(tri)

triangles.sort()
assert len(triangles) == 63, f"Expected 63 triangles (lines), got {len(triangles)}"
print(f"Recovered {len(triangles)} lines of GH(2,2)")

# Verify: each point on exactly 3 lines
for v in range(63):
    count = sum(1 for t in triangles if v in t)
    assert count == 3, f"Point {v} on {count} lines, expected 3"

# ============================================================
# Step 3: Build the dual Langer graph (concurrency graph on lines)
# Two lines (triangles) are adjacent iff they share a point.
# ============================================================
dual_nbrs = [[] for _ in range(63)]
for i in range(63):
    for j in range(i+1, 63):
        if set(triangles[i]) & set(triangles[j]):  # share a vertex
            dual_nbrs[i].append(j)
            dual_nbrs[j].append(i)

# Verify 6-regular
for i in range(63):
    assert len(dual_nbrs[i]) == 6, f"Dual vertex {i} has degree {len(dual_nbrs[i])}"

print("Dual Langer graph: 63 vertices, 6-regular")

# Verify NOT isomorphic to Langer (quick check: compare sorted degree sequences
# of distance-3 neighborhoods)
def distance_layers(nbrs, start):
    """BFS distance layers from start."""
    visited = {start}
    layers = [{start}]
    frontier = {start}
    while frontier:
        next_frontier = set()
        for v in frontier:
            for u in nbrs[v]:
                if u not in visited:
                    visited.add(u)
                    next_frontier.add(u)
        if next_frontier:
            layers.append(next_frontier)
        frontier = next_frontier
    return layers

# Check distance distribution from vertex 0
langer_layers = distance_layers(langer_nbrs, 0)
dual_layers = distance_layers(dual_nbrs, 0)
print(f"Langer distance distribution: {[len(l) for l in langer_layers]}")
print(f"Dual   distance distribution: {[len(l) for l in dual_layers]}")

# Check girth (smallest cycle)
def compute_girth(nbrs, n):
    """Compute girth by BFS from each vertex."""
    min_girth = float('inf')
    for s in range(n):
        dist = [-1]*n
        parent = [-1]*n
        dist[s] = 0
        queue = [s]
        qi = 0
        while qi < len(queue):
            v = queue[qi]; qi += 1
            for u in nbrs[v]:
                if dist[u] == -1:
                    dist[u] = dist[v] + 1
                    parent[u] = v
                    queue.append(u)
                elif parent[v] != u and parent[u] != v:
                    cycle_len = dist[v] + dist[u] + 1
                    min_girth = min(min_girth, cycle_len)
        if min_girth <= 4:
            return min_girth
    return min_girth

girth = compute_girth(dual_nbrs, 63)
print(f"Dual Langer girth: {girth}")

# ============================================================
# Step 4: Z_7 orbits on lines
#
# The companion matrix C of x^3+x+1 over F_2 acts on F_2^7 points
# and permutes Q(6,2). It acts on triangles by acting on their vertices.
# ============================================================

# Companion matrix of x^3+x+1 acting on (a0,a1,a2,b0,b1,b2,alpha) in F_2^7
# This is the same Z_7 generator used for the Langer graph.
def companion_action(v_idx):
    """Apply the Z_7 generator (companion matrix of x^3+x+1) to a Q(6,2) point index."""
    x = list(to_bit_vec(q62_indices[v_idx]))
    # The companion matrix acts as: (x0,x1,x2) -> (x1,x2,x0+x1)
    # and similarly on (x4,x5,x6), with x3 fixed under this action.
    # Actually, let's compute it properly.
    # From the Langer drawing, the orbits are known. Let me use the permutation directly.
    pass

# Instead of computing the matrix action, derive the permutation from the known orbits
langer_orbits = [
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

# Build the Z_7 permutation on the 63 points
z7_perm = [0]*63  # z7_perm[v] = image of v under the generator
for orb in langer_orbits:
    for k in range(7):
        z7_perm[orb[k]] = orb[(k+1) % 7]

# Now compute Z_7 action on lines (triangles)
# A triangle {a,b,c} maps to {z7(a), z7(b), z7(c)}
tri_set = {t: i for i, t in enumerate(triangles)}

z7_perm_lines = [0]*63
for i, tri in enumerate(triangles):
    image = tuple(sorted([z7_perm[v] for v in tri]))
    assert image in tri_set, f"Image {image} of triangle {tri} not found"
    z7_perm_lines[i] = tri_set[image]

# Compute orbits of Z_7 on lines
visited = [False]*63
line_orbits = []
for start in range(63):
    if visited[start]:
        continue
    orb = []
    v = start
    while not visited[v]:
        visited[v] = True
        orb.append(v)
        v = z7_perm_lines[v]
    line_orbits.append(orb)

print(f"\nZ_7 orbits on lines: {len(line_orbits)} orbits")
print(f"Orbit sizes: {sorted(len(o) for o in line_orbits)}")

# Verify all orbits have size 7
for orb in line_orbits:
    assert len(orb) == 7, f"Orbit has size {len(orb)}, expected 7"

orbit_map = [0]*63
for oi, orb in enumerate(line_orbits):
    for v in orb:
        orbit_map[v] = oi

# Quotient edges
q_edges = set()
for v in range(63):
    for u in dual_nbrs[v]:
        oi, oj = orbit_map[v], orbit_map[u]
        if oi < oj:
            q_edges.add((oi, oj))
q_edges = sorted(q_edges)
print(f"Quotient: {len(line_orbits)} vertices, {len(q_edges)} edges")

# Print orbit details
for oi, orb in enumerate(line_orbits):
    line_strs = [f"{{{triangles[v][0]},{triangles[v][1]},{triangles[v][2]}}}" for v in orb]
    print(f"  Orbit {oi}: lines {', '.join(line_strs)}")

# ============================================================
# Step 5: TikZ drawing — concentric heptagonal rings
# ============================================================

def angular_offset_search(orbits, nbrs, orbit_map, n_offsets=7):
    """Find angular offsets and reflections for each ring to minimize crossings."""
    n_orbits = len(orbits)
    best_offsets = [0] * n_orbits
    best_reflects = [False] * n_orbits

    for oi in range(1, n_orbits):
        best_score = float('inf')
        best_off = 0
        best_ref = False
        for off in range(7):
            for ref in [False, True]:
                offsets = list(best_offsets)
                reflects = list(best_reflects)
                offsets[oi] = off
                reflects[oi] = ref
                score = 0
                for v in range(63):
                    for u in nbrs[v]:
                        if v < u and orbit_map[v] != orbit_map[u]:
                            ov, ou = orbit_map[v], orbit_map[u]
                            pv = orbits[ov].index(v)
                            pu = orbits[ou].index(u)
                            if reflects[ov]:
                                pv = (7 - pv) % 7
                            if reflects[ou]:
                                pu = (7 - pu) % 7
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
                    best_ref = ref
        best_offsets[oi] = best_off
        best_reflects[oi] = best_ref
    return best_offsets, best_reflects

offsets, reflects = angular_offset_search(line_orbits, dual_nbrs, orbit_map)
print(f"\nOffsets:    {offsets}")
print(f"Reflects:  {reflects}")

def vertex_pos(v, scale=1.0):
    oi = orbit_map[v]
    pos_in_orbit = line_orbits[oi].index(v)
    if reflects[oi]:
        pos_in_orbit = (7 - pos_in_orbit) % 7
    r = 1.0 + oi * 0.65
    angle = 2 * math.pi * ((pos_in_orbit + offsets[oi]) % 7) / 7 + math.pi/2
    return (r * math.cos(angle) * scale, r * math.sin(angle) * scale)

# Count crossings
def segments_cross(p1, p2, p3, p4, eps=1e-9):
    """Test if segments p1-p2 and p3-p4 cross (proper intersection)."""
    def cross2d(a, b):
        return a[0]*b[1] - a[1]*b[0]
    d1 = (p2[0]-p1[0], p2[1]-p1[1])
    d2 = (p4[0]-p3[0], p4[1]-p3[1])
    denom = cross2d(d1, d2)
    if abs(denom) < eps:
        return False
    t = cross2d((p3[0]-p1[0], p3[1]-p1[1]), d2) / denom
    u = cross2d((p3[0]-p1[0], p3[1]-p1[1]), d1) / denom
    return eps < t < 1-eps and eps < u < 1-eps

edges = []
for v in range(63):
    for u in dual_nbrs[v]:
        if v < u:
            edges.append((v, u))

positions = [vertex_pos(v) for v in range(63)]
crossings = 0
for i in range(len(edges)):
    for j in range(i+1, len(edges)):
        e1, e2 = edges[i], edges[j]
        if set(e1) & set(e2):
            continue
        if segments_cross(positions[e1[0]], positions[e1[1]],
                         positions[e2[0]], positions[e2[1]]):
            crossings += 1

print(f"Edges: {len(edges)}, Crossings: {crossings}")

# ============================================================
# TikZ output
# ============================================================

lines = []
lines.append(r"\documentclass[border=12pt]{standalone}")
lines.append(r"\usepackage{tikz}")
lines.append(r"\begin{document}")
lines.append(r"\begin{tikzpicture}[")
lines.append(r"  edge/.style={thin, black!25},")
lines.append(r"  intedge/.style={thin, black!15},")
lines.append(r"]")

# Title
lines.append(r"\node[font=\Large\bfseries] at (0, 8.2) {Dual Langer graph (63 vertices)};")
lines.append(r"\node[font=\small, text=gray] at (0, 7.5) "
             r"{Line-collinearity graph of $\mathrm{GH}(2,2)$, "
             r"$Z_7$ layout (9 heptagonal rings)};")

# Draw edges
for v, u in edges:
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

outdir = "/tmp"
tex_path = os.path.join(outdir, "dual_langer_full.tex")
with open(tex_path, "w") as f:
    f.write("\n".join(lines))

# ============================================================
# TikZ: Quotient graph
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

lines2.append(r"\node[font=\Large\bfseries] at (0, 4.5) {$Z_7$ quotient of the dual Langer graph};")
lines2.append(r"\node[font=\small, text=gray] at (0, 3.9) {"
              + f"{len(line_orbits)} vertices (orbits), {len(q_edges)} edges"
              + r"};")

qpos = {}
for i in range(len(line_orbits)):
    angle = 2 * math.pi * i / len(line_orbits) + math.pi/2
    r = 2.5
    qpos[i] = (r * math.cos(angle), r * math.sin(angle))

for (a, b) in q_edges:
    x1, y1 = qpos[a]
    x2, y2 = qpos[b]
    lines2.append(f"\\draw[qedge] ({x1:.3f}, {y1:.3f}) -- ({x2:.3f}, {y2:.3f});")

for i in range(len(line_orbits)):
    x, y = qpos[i]
    c = colors[i].split("!")[0]
    lines2.append(f"\\node[qnode={c}] at ({x:.3f}, {y:.3f}) {{{i}}};")

lines2.append(r"\end{tikzpicture}")
lines2.append(r"\end{document}")

tex_path2 = os.path.join(outdir, "dual_langer_quotient.tex")
with open(tex_path2, "w") as f:
    f.write("\n".join(lines2))

# ============================================================
# Compile both
# ============================================================

for name in ["dual_langer_full", "dual_langer_quotient"]:
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
        print(f"  {jpg}")
