#!/usr/bin/env python3
"""Draw the Zhou Z₂ quotient — Sab(PSL(2,13), D₁₂), 91 vertices, 6-regular.

Two layouts:
  Z₇:  13 concentric heptagonal rings
  Z₁₃: 7 concentric 13-gonal rings

Uses precomputed trig tables and numpy-vectorized crossing counts.
"""

import math
import subprocess
import random
import numpy as np

N = 91

# ============================================================
# Adjacency list (1-indexed, from GAP)
# ============================================================
ADJ = [
    [4, 23, 30, 39, 54, 85], [5, 8, 12, 36, 57, 78],
    [6, 17, 19, 39, 43, 64], [1, 12, 25, 27, 50, 85],
    [2, 17, 30, 33, 56, 78], [3, 15, 25, 36, 46, 64],
    [8, 16, 21, 22, 48, 89], [2, 7, 31, 36, 77, 89],
    [10, 27, 31, 32, 60, 65], [9, 21, 29, 34, 40, 60],
    [12, 14, 16, 34, 51, 76], [2, 4, 11, 27, 51, 57],
    [14, 28, 38, 39, 41, 91], [11, 13, 33, 34, 83, 91],
    [6, 16, 25, 28, 69, 80], [7, 11, 15, 22, 69, 76],
    [3, 5, 22, 33, 43, 44], [19, 24, 26, 29, 55, 75],
    [3, 18, 35, 39, 67, 75], [21, 33, 35, 37, 42, 81],
    [7, 10, 20, 29, 42, 48], [7, 16, 17, 24, 44, 79],
    [1, 24, 30, 34, 74, 87], [18, 22, 23, 26, 79, 87],
    [4, 6, 15, 26, 50, 53], [18, 24, 25, 31, 53, 66],
    [4, 9, 12, 38, 65, 68], [13, 15, 29, 38, 49, 80],
    [10, 18, 21, 28, 49, 55], [1, 5, 23, 32, 56, 59],
    [8, 9, 26, 32, 66, 77], [9, 30, 31, 35, 59, 82],
    [5, 14, 17, 20, 81, 83], [10, 11, 14, 23, 40, 74],
    [19, 20, 32, 37, 67, 82], [2, 6, 8, 37, 46, 62],
    [20, 35, 36, 38, 62, 90], [13, 27, 28, 37, 68, 90],
    [1, 3, 13, 19, 41, 54], [10, 34, 43, 62, 63, 86],
    [13, 39, 45, 48, 63, 66], [20, 21, 53, 54, 58, 73],
    [3, 17, 40, 58, 68, 86], [17, 22, 47, 52, 60, 90],
    [41, 46, 48, 51, 82, 87], [6, 36, 45, 60, 84, 87],
    [44, 51, 54, 55, 77, 90], [7, 21, 41, 45, 50, 72],
    [28, 29, 57, 59, 61, 70], [4, 25, 48, 61, 72, 83],
    [11, 12, 45, 47, 55, 82], [44, 53, 57, 60, 67, 91],
    [25, 26, 42, 52, 58, 91], [1, 39, 42, 47, 73, 77],
    [18, 29, 47, 51, 56, 88], [5, 30, 55, 63, 69, 88],
    [2, 12, 49, 52, 67, 70], [42, 43, 53, 59, 68, 76],
    [30, 32, 49, 58, 61, 76], [9, 10, 44, 46, 52, 84],
    [49, 50, 59, 62, 79, 83], [36, 37, 40, 61, 63, 79],
    [40, 41, 56, 62, 66, 69], [3, 6, 65, 70, 74, 88],
    [9, 27, 64, 73, 79, 88], [26, 31, 41, 63, 70, 81],
    [19, 35, 52, 57, 69, 71], [27, 38, 43, 58, 71, 87],
    [15, 16, 56, 63, 67, 71], [49, 57, 64, 66, 74, 81],
    [67, 68, 69, 77, 83, 87], [48, 50, 74, 75, 78, 90],
    [42, 54, 65, 78, 79, 80], [23, 34, 64, 70, 72, 90],
    [18, 19, 72, 76, 78, 84], [11, 16, 58, 59, 75, 84],
    [8, 31, 47, 54, 71, 83], [2, 5, 72, 73, 75, 80],
    [22, 24, 61, 62, 65, 73], [15, 28, 73, 78, 82, 86],
    [20, 33, 66, 70, 84, 85], [32, 35, 45, 51, 80, 86],
    [14, 33, 50, 61, 71, 77], [46, 60, 75, 76, 81, 85],
    [1, 4, 81, 84, 86, 89], [40, 43, 80, 82, 85, 89],
    [23, 24, 45, 46, 68, 71], [55, 56, 64, 65, 89, 91],
    [7, 8, 85, 86, 88, 91], [37, 38, 44, 47, 72, 74],
    [13, 14, 52, 53, 88, 89],
]

# Edge list (0-indexed) as numpy array for vectorized operations
EDGES = []
for v in range(N):
    for u in ADJ[v]:
        if v < u - 1:
            EDGES.append((v, u - 1))
EDGES = np.array(EDGES, dtype=np.int32)  # shape (273, 2)

# ============================================================
# Z₇ orbits: 13 orbits of 7 (ordered under generator)
# ============================================================
Z7_ORBITS = [
    [1, 2, 86, 85, 30, 77, 87], [3, 43, 66, 82, 10, 6, 89],
    [4, 90, 5, 25, 68, 78, 59], [7, 17, 67, 20, 57, 28, 32],
    [8, 71, 49, 35, 23, 50, 61], [9, 52, 73, 11, 69, 34, 16],
    [12, 81, 37, 51, 62, 54, 76], [13, 53, 19, 84, 47, 91, 40],
    [14, 44, 79, 56, 60, 63, 58], [15, 72, 41, 26, 22, 24, 88],
    [18, 29, 45, 80, 75, 39, 55], [21, 36, 48, 65, 70, 64, 38],
    [27, 83, 31, 74, 42, 46, 33],
]

# ============================================================
# Z₁₃ orbits: 7 orbits of 13 (ordered under generator)
# ============================================================
Z13_ORBITS = [
    [1, 2, 20, 45, 34, 76, 61, 29, 57, 24, 70, 83, 63],
    [3, 22, 28, 65, 88, 90, 46, 25, 72, 50, 12, 10, 62],
    [4, 15, 68, 87, 58, 35, 74, 52, 27, 64, 79, 56, 21],
    [5, 33, 73, 48, 32, 69, 80, 84, 55, 13, 30, 9, 42],
    [6, 39, 47, 37, 40, 18, 8, 67, 86, 82, 81, 91, 51],
    [7, 66, 85, 53, 26, 43, 17, 71, 49, 14, 59, 36, 11],
    [16, 77, 78, 44, 38, 60, 31, 23, 41, 19, 75, 89, 54],
]


# ============================================================
# Precomputed trig tables and coordinate builder
# ============================================================

def precompute_trig(k):
    """Precompute cos/sin for all k angular positions (with pi/2 offset)."""
    cos_table = np.zeros(k)
    sin_table = np.zeros(k)
    for i in range(k):
        theta = 2 * math.pi * i / k + math.pi / 2
        cos_table[i] = math.cos(theta)
        sin_table[i] = math.sin(theta)
    return cos_table, sin_table


def compute_radii(n_rings):
    """Radii for concentric rings."""
    R_max, R_min = 7.5, 2.2
    if n_rings == 1:
        return np.array([R_max])
    return np.array([R_max - i * (R_max - R_min) / (n_rings - 1)
                      for i in range(n_rings)])


def build_orbit_arrays(orbits):
    """Build arrays: for each vertex (0-indexed), its orbit and position."""
    orbit_of = np.zeros(N, dtype=np.int32)
    pos_of = np.zeros(N, dtype=np.int32)
    for oi, orb in enumerate(orbits):
        for pi, v in enumerate(orb):
            orbit_of[v - 1] = oi
            pos_of[v - 1] = pi
    return orbit_of, pos_of


def make_coords_np(orbit_of, pos_of, offsets, radii, cos_table, sin_table, k):
    """Vectorized coordinate computation using precomputed trig."""
    angular_pos = (pos_of + offsets[orbit_of]) % k
    r = radii[orbit_of]
    x = r * cos_table[angular_pos]
    y = r * sin_table[angular_pos]
    return np.column_stack([x, y])  # shape (N, 2)


# ============================================================
# Vectorized crossing count
# ============================================================

def count_crossings_np(coords, edges):
    """Count edge crossings using numpy vectorization."""
    n_edges = len(edges)
    p1 = coords[edges[:, 0]]  # (E, 2)
    p2 = coords[edges[:, 1]]  # (E, 2)

    crossings = 0
    for i in range(n_edges):
        ei = edges[i]
        rest = edges[i + 1:]
        mask = ((rest[:, 0] != ei[0]) & (rest[:, 0] != ei[1]) &
                (rest[:, 1] != ei[0]) & (rest[:, 1] != ei[1]))
        if not mask.any():
            continue

        # Edge i endpoints (broadcast as row vectors)
        ax, ay = p1[i, 0], p1[i, 1]
        bx, by = p2[i, 0], p2[i, 1]

        # Candidate edges endpoints
        cx = p1[i + 1:, 0][mask]
        cy = p1[i + 1:, 1][mask]
        dx = p2[i + 1:, 0][mask]
        dy = p2[i + 1:, 1][mask]

        # cross(cd, ca) and cross(cd, cb)
        d1 = (dx - cx) * (ay - cy) - (dy - cy) * (ax - cx)
        d2 = (dx - cx) * (by - cy) - (dy - cy) * (bx - cx)
        # cross(ab, ac) and cross(ab, ad)
        d3 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        d4 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax)

        denom1 = d1 - d2
        denom2 = d3 - d4
        valid = (np.abs(denom1) > 1e-12) & (np.abs(denom2) > 1e-12)
        if not valid.any():
            continue

        t = np.where(valid, d1 / np.where(valid, denom1, 1.0), 0.0)
        s = np.where(valid, d3 / np.where(valid, denom2, 1.0), 0.0)

        EPS = 1e-9
        cross_mask = valid & (t > EPS) & (t < 1 - EPS) & (s > EPS) & (s < 1 - EPS)
        crossings += int(cross_mask.sum())

    return crossings


# ============================================================
# Edge-length proxy (very fast)
# ============================================================

def edge_length_sum(coords, edges):
    """Total edge length — fast proxy for layout quality."""
    diff = coords[edges[:, 0]] - coords[edges[:, 1]]
    return np.sqrt((diff ** 2).sum(axis=1)).sum()


# ============================================================
# Greedy offset search with precomputed trig
# ============================================================

def greedy_offsets(orbits, orbit_of, pos_of, radii, cos_table, sin_table, k):
    """Greedy per-ring offset search using edge-length-squared proxy."""
    n_rings = len(orbits)
    offsets = np.zeros(n_rings, dtype=np.int32)

    for ring in range(1, n_rings):
        best_score = float('inf')
        best_off = 0
        for off in range(k):
            offsets[ring] = off
            coords = make_coords_np(orbit_of, pos_of, offsets, radii,
                                     cos_table, sin_table, k)
            # Only count edges involving rings <= current
            mask = (orbit_of[EDGES[:, 0]] <= ring) & (orbit_of[EDGES[:, 1]] <= ring)
            if not mask.any():
                continue
            sub_edges = EDGES[mask]
            diff = coords[sub_edges[:, 0]] - coords[sub_edges[:, 1]]
            score = (diff ** 2).sum()
            if score < best_score:
                best_score = score
                best_off = off
        offsets[ring] = best_off
    return offsets


# ============================================================
# Search: greedy + sampled reflections
# ============================================================

def search_layout(orbits, n_mirror_samples=None):
    """Find best offsets + reflections for a set of orbits."""
    k = len(orbits[0])
    n_rings = len(orbits)
    radii = compute_radii(n_rings)
    cos_table, sin_table = precompute_trig(k)

    total_mirrors = 1 << (n_rings - 1)
    if n_mirror_samples is None:
        # Exhaustive if ≤ 256, else sample 256
        n_mirror_samples = min(total_mirrors, 256)

    if n_mirror_samples >= total_mirrors:
        mirror_bits_list = range(total_mirrors)
        print(f"  Exhaustive reflection sweep: {total_mirrors} combos")
    else:
        # Always include 0 (no reflections), then random sample
        mirror_bits_list = [0] + [random.randint(1, total_mirrors - 1)
                                   for _ in range(n_mirror_samples - 1)]
        print(f"  Sampled reflection sweep: {n_mirror_samples}/{total_mirrors} combos")

    best_proxy = float('inf')
    best_offsets = None
    best_seqs = None

    for mirror_bits in mirror_bits_list:
        # Build trial sequences with reflections
        trial_seqs = [orbits[0]]
        for i in range(n_rings - 1):
            if mirror_bits & (1 << i):
                trial_seqs.append(orbits[i + 1][::-1])
            else:
                trial_seqs.append(orbits[i + 1])

        orbit_of, pos_of = build_orbit_arrays(trial_seqs)
        offsets = greedy_offsets(trial_seqs, orbit_of, pos_of, radii,
                                 cos_table, sin_table, k)
        coords = make_coords_np(orbit_of, pos_of, offsets, radii,
                                 cos_table, sin_table, k)
        proxy = edge_length_sum(coords, EDGES)

        if proxy < best_proxy:
            best_proxy = proxy
            best_offsets = offsets.copy()
            best_seqs = list(trial_seqs)

    # Compute actual crossings for the winner
    orbit_of, pos_of = build_orbit_arrays(best_seqs)
    coords = make_coords_np(orbit_of, pos_of, best_offsets, radii,
                             cos_table, sin_table, k)
    crossings = count_crossings_np(coords, EDGES)
    print(f"  Best proxy: {best_proxy:.1f}, crossings: {crossings}")

    return best_offsets, best_seqs, crossings, coords


# ============================================================
# Exhaustive search for Z₁₃ (13^6 ≈ 4.8M with rotations+reflections)
# ============================================================

def exhaustive_z13(orbits):
    """Full exhaustive search for Z₁₃: (2×13)^6 is too big, but
    greedy+reflections with actual crossing verification on top-N candidates."""
    from itertools import product
    k = len(orbits[0])  # 13
    n_rings = len(orbits)  # 7
    radii = compute_radii(n_rings)
    cos_table, sin_table = precompute_trig(k)

    # Phase 1: find top-50 (offset, reflection) combos by proxy
    total_mirrors = 1 << (n_rings - 1)  # 64
    candidates = []

    print(f"  Phase 1: greedy search over {total_mirrors} reflection combos")
    for mirror_bits in range(total_mirrors):
        trial_seqs = [orbits[0]]
        for i in range(n_rings - 1):
            if mirror_bits & (1 << i):
                trial_seqs.append(orbits[i + 1][::-1])
            else:
                trial_seqs.append(orbits[i + 1])

        orbit_of, pos_of = build_orbit_arrays(trial_seqs)
        offsets = greedy_offsets(trial_seqs, orbit_of, pos_of, radii,
                                 cos_table, sin_table, k)
        coords = make_coords_np(orbit_of, pos_of, offsets, radii,
                                 cos_table, sin_table, k)
        proxy = edge_length_sum(coords, EDGES)
        candidates.append((proxy, offsets.copy(), list(trial_seqs)))

    # Phase 2: verify top 20 by actual crossing count
    candidates.sort()
    print(f"  Phase 2: crossing count on top 20 candidates")
    best_cross = float('inf')
    best_offsets = None
    best_seqs = None
    best_coords = None

    for proxy, offsets, seqs in candidates[:20]:
        orbit_of, pos_of = build_orbit_arrays(seqs)
        coords = make_coords_np(orbit_of, pos_of, offsets, radii,
                                 cos_table, sin_table, k)
        cross = count_crossings_np(coords, EDGES)
        if cross < best_cross:
            best_cross = cross
            best_offsets = offsets
            best_seqs = seqs
            best_coords = coords

    # Phase 3: local search — try all k offsets for each ring of the winner
    print(f"  Phase 3: local offset refinement from {best_cross} crossings")
    improved = True
    while improved:
        improved = False
        orbit_of, pos_of = build_orbit_arrays(best_seqs)
        for ring in range(1, n_rings):
            for off in range(k):
                if off == best_offsets[ring]:
                    continue
                trial_offsets = best_offsets.copy()
                trial_offsets[ring] = off
                coords = make_coords_np(orbit_of, pos_of, trial_offsets, radii,
                                         cos_table, sin_table, k)
                cross = count_crossings_np(coords, EDGES)
                if cross < best_cross:
                    best_cross = cross
                    best_offsets = trial_offsets
                    best_coords = coords
                    improved = True
                    print(f"    Ring {ring} offset {off}: {cross} crossings")

    print(f"  Final: {best_cross} crossings")
    return best_offsets, best_seqs, best_cross, best_coords


# ============================================================
# TikZ generation
# ============================================================

# Ring palette (same as draw_langer.py / draw_triality_hexagon.py):
# blue outer -> red inner. Orbit 0 is the innermost ring, so reverse.
RING_PALETTE = [
    "0033CC", "0066FF", "0099FF", "00CCCC", "00CC66",
    "33CC00", "99CC00", "CCCC00", "FF9900", "FF6600",
    "FF3300", "CC0000",
]


def generate_tikz(coords_np, orbits, title, subtitle, crossings):
    n_rings = len(orbits)
    orbit_of, _ = build_orbit_arrays(orbits)

    # Rank orbits by mean radius: outermost ring -> ring0 (blue), innermost -> red
    mean_r = []
    for orb in orbits:
        mean_r.append(sum(math.hypot(*coords_np[v - 1]) for v in orb) / len(orb))
    rank = {oi: r for r, oi in enumerate(sorted(range(n_rings), key=lambda oi: -mean_r[oi]))}

    def ring_color_idx(oi):
        return min(11, int(rank[oi] * 12 / n_rings))

    lines = []
    lines.append(r"\documentclass[border=12pt]{standalone}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    for idx, hex_c in enumerate(RING_PALETTE):
        lines.append(f"\\definecolor{{ring{idx}}}{{HTML}}{{{hex_c}}}")
    lines.append(r"\definecolor{gap0}{HTML}{E53935}")   # same ring: red
    lines.append(r"\definecolor{gap1}{HTML}{1E88E5}")   # nearby rings: blue
    lines.append(r"\definecolor{gapfar}{HTML}{CCCCCC}") # far: light grey
    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}")

    y_title = coords_np[:, 1].max() + 1.5
    lines.append(f"\\node[font=\\Large\\bfseries] at (0, {y_title:.1f}) {{{title}}};")
    lines.append(f"\\node[font=\\small, text=gray] at (0, {y_title - 0.7:.1f}) {{{subtitle}}};")
    lines.append(f"\\node[font=\\footnotesize, text=gray!70] at (0, {y_title - 1.2:.1f}) {{{crossings} crossings}};")

    # Edges by ring gap (far edges first, close edges on top — dmt style,
    # same convention as clayworth/draw_langer.py): same ring red, nearby
    # rings blue, far rings light grey.
    edges_by_gap = {}
    for v, u in EDGES:
        gap = abs(rank[int(orbit_of[v])] - rank[int(orbit_of[u])])
        edges_by_gap.setdefault(gap, []).append((v, u))

    for gap in sorted(edges_by_gap.keys(), reverse=True):
        if gap == 0:
            style = "line width=0.5pt, gap0, opacity=0.6"
        elif gap <= 3:
            style = "line width=0.4pt, gap1, opacity=0.3"
        else:
            style = "line width=0.35pt, gapfar, opacity=0.45"
        for v, u in edges_by_gap[gap]:
            x1, y1 = coords_np[v]
            x2, y2 = coords_np[u]
            lines.append(f"\\draw[{style}] ({x1:.3f}, {y1:.3f}) -- ({x2:.3f}, {y2:.3f});")

    # Vertices colored by ring (blue outer -> red inner)
    for v in range(N):
        x, y = coords_np[v]
        c = f"ring{ring_color_idx(int(orbit_of[v]))}"
        lines.append(f"\\fill[{c}] ({x:.3f}, {y:.3f}) circle (2.2pt);")

    lines.append(r"\end{tikzpicture}")
    lines.append(r"\end{document}")
    return "\n".join(lines)


def compile_and_convert(name, tikz_code):
    tex_path = f"/tmp/{name}.tex"
    with open(tex_path, "w") as f:
        f.write(tikz_code)
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", f"{name}.tex"],
        capture_output=True, text=True, cwd="/tmp"
    )
    if result.returncode != 0:
        print(f"LaTeX error for {name}!")
        print(result.stdout[-500:])
        return None
    pdf = f"/tmp/{name}.pdf"
    jpg = f"/tmp/{name}.jpg"
    subprocess.run(
        ["sips", "-s", "format", "jpeg", "-s", "formatOptions", "95",
         "--resampleWidth", "2400", pdf, "--out", jpg],
        capture_output=True
    )
    print(f"  -> {jpg}")
    return jpg


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    import sys
    random.seed(42)

    print("=== Z₇ layout: 13 heptagonal rings ===")
    z7_offs, z7_seqs, z7_cross, z7_coords = exhaustive_z13(Z7_ORBITS)
    z7_tikz = generate_tikz(
        z7_coords, z7_seqs,
        r"Zhou $Z_2$ quotient --- $Z_7$ layout",
        r"$\mathrm{Sab}(\mathrm{PSL}(2,13), D_{12})$, 91 vertices, 6-regular, 13 heptagonal rings",
        z7_cross,
    )
    compile_and_convert("zhou_quotient_z7", z7_tikz)

    print()
    print("=== Z₁₃ layout: 7 tridecagonal rings ===")
    z13_offs, z13_seqs, z13_cross, z13_coords = exhaustive_z13(Z13_ORBITS)
    z13_tikz = generate_tikz(
        z13_coords, z13_seqs,
        r"Zhou $Z_2$ quotient --- $Z_{13}$ layout",
        r"$\mathrm{Sab}(\mathrm{PSL}(2,13), D_{12})$, 91 vertices, 6-regular, 7 tridecagonal rings",
        z13_cross,
    )
    compile_and_convert("zhou_quotient_z13", z13_tikz)

    print()
    print(f"Z₇:  {z7_cross} crossings")
    print(f"Z₁₃: {z13_cross} crossings")
    print(f"Better: {'Z₇' if z7_cross <= z13_cross else 'Z₁₃'}")
