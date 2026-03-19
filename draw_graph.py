#!/usr/bin/env python3
"""
Symmetric graph drawing pipeline.

Uses Eades-Hong theory to detect whether a graph should be drawn with:
  - Cyclic symmetry (single circle, interleaved orbits)
  - Cylindrical symmetry (concentric rings, one per orbit)

Then outputs TikZ, compiles PDF, opens it.

Usage:
    python3 draw_graph.py petersen
    python3 draw_graph.py dodecahedron --layout rings
    python3 draw_graph.py heawood --layout cyclic
    python3 draw_graph.py desargues
"""

import subprocess
import sys
import os
import re
import math

GRAPHVIZ_DIR = os.path.expanduser("~/git/graphviz")
GAP_SCRIPT = os.path.join(GRAPHVIZ_DIR, "symmetric-layout.g")

COLORS = [
    "E53935", "1E88E5", "43A047", "FB8C00", "8E24AA",
    "00ACC1", "F4511E", "3949AB", "7CB342", "C2185B",
    "FF7043", "26A69A", "5C6BC0", "FFCA28", "AB47BC",
]


# ============================================================
# GAP interface
# ============================================================

def run_gap(graph_name, max_deg=4):
    """Run GAP and capture the JSON output."""
    gap_input = f'graph_name := "{graph_name}";\nmax_deg := {max_deg};\nRead("{GAP_SCRIPT}");\n'
    result = subprocess.run(
        ["gap", "-q"], input=gap_input,
        capture_output=True, text=True, timeout=120
    )
    output = result.stdout + result.stderr
    match = re.search(r'=== JSON OUTPUT ===\n(.*?)\n=== END JSON ===', output, re.DOTALL)
    if not match:
        raise RuntimeError(f"No JSON in GAP output for {graph_name}")
    json_text = match.group(1).replace("true", "True").replace("false", "False")
    return eval(json_text)


# ============================================================
# Quotient selection
# ============================================================

def score_quotient(q):
    score = 0
    if q.get('is_cycle'):
        score += 100
    score -= q['max_valency'] * 10
    n_orbs = q['n_orbits']
    if 3 <= n_orbs <= 10:
        score += 20
    elif n_orbs == 2:
        score += 15  # 2 equal orbits = generalized Petersen
    sizes = q['orbit_sizes']
    if len(set(sizes)) == 1:
        score += 30
    elif len(set(sizes)) <= 2:
        score += 15
    if q.get('has_internal_cycle'):
        score += sum(q['has_internal_cycle']) * 5
    # Bonus for 2 equal orbits — generalized Petersen pattern
    # In GP(n,k): outer ring forms a cycle, inner ring forms a star (not a cycle)
    # This is the canonical drawing for GP graphs, so score it highest
    if n_orbs == 2 and len(set(sizes)) == 1:
        has_cycle = q.get('has_internal_cycle', [])
        if len(has_cycle) >= 2 and (has_cycle[0] or has_cycle[1]):
            score += 150  # GP pattern — outscores everything else
    return score


def select_best(results):
    if not results:
        return None
    scored = [(score_quotient(r), i, r) for i, r in enumerate(results)]
    scored.sort(reverse=True)
    print(f"\nTop quotients:")
    for score, idx, r in scored[:3]:
        print(f"  score={score:4d}  {r.get('label','?')}: "
              f"{r['n_orbits']} orbits {r['orbit_sizes']}, "
              f"max_deg={r['max_valency']}, cycle={r.get('is_cycle',False)}")
    return scored[0][2]


# ============================================================
# Layout: detect cyclic vs cylindrical
# ============================================================

def detect_symmetry_type(quotient, adj):
    """
    Detect whether the graph should be drawn with cyclic or cylindrical symmetry.

    Cylindrical (concentric rings) if:
      - The quotient is a path (each orbit connects only to adjacent orbits)
      - Or there are exactly 2 equal-sized orbits (generalized Petersen pattern)

    Cyclic (single circle, interleaved) otherwise.
    """
    n_orbits = quotient['n_orbits']
    quot_adj = quotient.get('quotient_adj', [])
    valencies = quotient.get('quotient_valencies', [])
    sizes = quotient['orbit_sizes']

    # Check if quotient is a path
    if n_orbits >= 3:
        deg_1_count = sum(1 for v in valencies if v == 1)
        deg_2_count = sum(1 for v in valencies if v == 2)
        if deg_1_count == 2 and deg_2_count == n_orbits - 2:
            return "cylindrical"

    # 2 equal orbits with edge between them → generalized Petersen style
    if n_orbits == 2 and len(set(sizes)) == 1:
        return "cylindrical"

    # 2-4 orbits, all equal size, max degree ≤ 2 → cylindrical often looks better
    if n_orbits <= 4 and len(set(sizes)) == 1 and max(valencies) <= 2:
        return "cylindrical"

    return "cyclic"


def find_path_ordering(quotient):
    """
    If quotient graph is a path, return orbit indices in path order.
    Prefers endpoint with internal cycle as outer ring (Schlegel aesthetics).
    Returns None if quotient is not a path.
    """
    n_orbits = quotient['n_orbits']
    if n_orbits <= 2:
        return None

    valencies = quotient.get('quotient_valencies', [])
    quot_adj = quotient.get('quotient_adj', [])

    if not valencies or not quot_adj:
        return None

    endpoints = [i for i in range(n_orbits) if valencies[i] == 1]
    if len(endpoints) != 2:
        return None

    # Prefer starting from an endpoint with internal cycle (pentagon → outer ring)
    has_cycle = quotient.get('has_internal_cycle', [])
    start = endpoints[0]
    if has_cycle and len(has_cycle) > max(endpoints):
        if has_cycle[endpoints[1]] and not has_cycle[endpoints[0]]:
            start = endpoints[1]

    # Traverse the path
    path = [start]
    visited = {start}
    while len(path) < n_orbits:
        current = path[-1]
        neighbors = [nb - 1 for nb in quot_adj[current]]  # 1-indexed → 0-indexed
        advanced = False
        for nb in neighbors:
            if nb not in visited:
                path.append(nb)
                visited.add(nb)
                advanced = True
                break
        if not advanced:
            break

    if len(path) != n_orbits:
        return None

    return path


def get_orbit_sequence(quotient, orbit_idx, orbit):
    """Get the ordered sequence of vertices in an orbit under the generator."""
    data = quotient.get('orbit_order', None)
    if data and orbit_idx < len(data):
        seq = data[orbit_idx]
        if set(seq) == set(orbit):
            return seq
    return sorted(orbit)


# ============================================================
# Layout: Eades-Hong rotation (cyclic)
# ============================================================

def find_hamiltonian_cycle(adj, n):
    """Check if 1-2-...-n-1 is a Hamiltonian cycle, or find one by DFS."""
    # Check natural ordering
    is_ham = True
    for i in range(n):
        next_v = (i % n) + 1  # next vertex in 1-indexed
        curr_v = i + 1
        if next_v not in adj[i]:
            is_ham = False
            break
    if is_ham:
        return list(range(1, n + 1))

    if n > 20:
        return None

    # DFS for small graphs
    neighbors = [set() for _ in range(n)]
    for i in range(n):
        for j in adj[i]:
            neighbors[i].add(j - 1)

    path = [0]
    visited = {0}

    def dfs():
        if len(path) == n:
            return 0 in neighbors[path[-1]]
        for nb in sorted(neighbors[path[-1]]):
            if nb not in visited:
                path.append(nb)
                visited.add(nb)
                if dfs():
                    return True
                path.pop()
                visited.discard(nb)
        return False

    if dfs():
        return [v + 1 for v in path]
    return None


def layout_cyclic(quotient, n, adj):
    """
    Cyclic layout with optimized vertex ordering.

    For equal-sized orbits: Eades-Hong interleaving with rotation+reflection
    optimization across all orbits to minimize crossings.
    For unequal orbits: Hamiltonian cycle fallback.
    """
    orbits = quotient['orbits']
    m = len(orbits)
    R = 4.5

    # If orbits are unequal sizes, use Hamiltonian cycle
    sizes = [len(o) for o in orbits]
    if len(set(sizes)) > 1:
        ham = find_hamiltonian_cycle(adj, n)
        if ham:
            print(f"  Using Hamiltonian cycle (unequal orbits)")
            coords = {}
            for i, v in enumerate(ham):
                theta = math.radians(90 + i * (360.0 / n))
                coords[v] = (R * math.cos(theta), R * math.sin(theta))
            return coords

    # Get orbit sequences
    sequences = []
    for i in range(m):
        sequences.append(get_orbit_sequence(quotient, i, orbits[i]))

    k = len(sequences[0])  # orbit size (all equal)

    def build_cyclic_coords(seqs, rotations):
        """Place interleaved orbits on one circle with given rotations."""
        coords = {}
        for i in range(m):
            seq = seqs[i]
            ki = len(seq)
            rot = rotations[i] % ki
            for j in range(ki):
                v = seq[(j + rot) % ki]
                position = i + j * m
                theta = math.radians(90 + position * (360.0 / n))
                coords[v] = (R * math.cos(theta), R * math.sin(theta))
        return coords

    # Optimize: try all rotation+reflection combinations
    # Orbit 0 fixed, orbits 1..m-1 get k rotations × 2 reflections each
    from itertools import product

    total = (2 * k) ** (m - 1)
    if total <= 20000:
        best_coords = None
        best_crossings = float('inf')

        options = list(range(2 * k))  # 0..k-1 = rotations, k..2k-1 = reflected+rotations

        for combo in product(*([options] * (m - 1))):
            trial_seqs = [sequences[0]]
            trial_rots = [0]
            for c in combo:
                orbit_idx = len(trial_seqs)
                if c < k:
                    trial_seqs.append(sequences[orbit_idx])
                    trial_rots.append(c)
                else:
                    trial_seqs.append(sequences[orbit_idx][::-1])
                    trial_rots.append(c - k)

            trial_coords = build_cyclic_coords(trial_seqs, trial_rots)
            crossings = count_all_crossings(trial_coords, adj, n)

            if crossings < best_crossings:
                best_crossings = crossings
                best_coords = trial_coords

        print(f"  Cyclic optimization: {total} combos → {best_crossings} crossings")
        return best_coords
    else:
        # Standard Eades-Hong interleaving
        coords = {}
        for i in range(m):
            seq = sequences[i]
            for j, v in enumerate(seq):
                position = i + j * m
                theta = math.radians(90 + position * (360.0 / n))
                coords[v] = (R * math.cos(theta), R * math.sin(theta))
        return coords


# ============================================================
# Geometric crossing count
# ============================================================

def segments_cross(p1, p2, p3, p4):
    """
    Test if line segments p1-p2 and p3-p4 cross.

    Counts both proper crossings (interior of both segments) and
    degenerate crossings (vertex on other edge). This is important
    for the optimizer: degenerate layouts where vertices lie on
    other edges should be penalized, not rewarded.
    """
    def cross2d(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    d1 = cross2d(p3, p4, p1)
    d2 = cross2d(p3, p4, p2)
    d3 = cross2d(p1, p2, p3)
    d4 = cross2d(p1, p2, p4)

    denom1 = d1 - d2
    denom2 = d3 - d4

    if abs(denom1) < 1e-12 or abs(denom2) < 1e-12:
        return False  # Parallel or degenerate

    t = d1 / denom1  # Parameter along p1-p2
    s = d3 / denom2  # Parameter along p3-p4

    EPS = 1e-9
    return EPS < t < 1 - EPS and EPS < s < 1 - EPS


def count_all_crossings(coords, adj, n, perturb=True):
    """
    Count total edge crossings using geometric segment intersection.

    perturb: apply deterministic micro-perturbation to break vertex-on-edge
    degeneracies. Without this, layouts where vertices lie on non-incident
    edges can report artificially low crossing counts.
    """
    if perturb:
        # Deterministic perturbation: unique tiny offset per vertex
        eps = 1e-6
        coords = {v: (x + eps * math.sin(v * 17.3),
                       y + eps * math.cos(v * 23.7))
                  for v, (x, y) in coords.items()}

    edges = []
    for v in range(1, n + 1):
        for w in adj[v - 1]:
            if v < w:
                edges.append((v, w))

    crossings = 0
    for i in range(len(edges)):
        u1, v1 = edges[i]
        for j in range(i + 1, len(edges)):
            u2, v2 = edges[j]
            if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                continue
            if segments_cross(coords[u1], coords[v1], coords[u2], coords[v2]):
                crossings += 1
    return crossings


# ============================================================
# Layout: cylindrical (concentric rings)
# ============================================================

def count_crossings_between_rings(seq_outer, seq_inner, adj, offset_steps):
    """
    Count edge crossings between two concentric orbit rings.

    Buchheim & Hong insight: for edges between two concentric circles,
    two edges (u1→v1) and (u2→v2) cross iff the angular intervals
    they span on the annulus overlap in a crossing pattern.

    offset_steps: how many positions to rotate the inner ring relative to outer.
    """
    k = len(seq_outer)
    if k == 0:
        return 0

    # Build position maps
    pos_outer = {v: i for i, v in enumerate(seq_outer)}
    pos_inner = {v: (i + offset_steps) % k for i, v in enumerate(seq_inner)}

    # Find inter-orbit edges
    inter_edges = []
    outer_set = set(seq_outer)
    inner_set = set(seq_inner)
    for u in seq_outer:
        for v in adj[u - 1]:
            if v in inner_set:
                inter_edges.append((pos_outer[u], pos_inner[v]))

    # Count crossings: edges (a1, b1) and (a2, b2) cross on concentric circles
    # iff the arcs (a1→b1) and (a2→b2) interleave when mapped to the circle
    crossings = 0
    for i in range(len(inter_edges)):
        a1, b1 = inter_edges[i]
        for j in range(i + 1, len(inter_edges)):
            a2, b2 = inter_edges[j]
            # Two edges between concentric circles cross iff
            # (a1 - a2) and (b1 - b2) have opposite signs modulo k
            # More precisely: they cross iff one separates the other
            da = (a2 - a1) % k
            db = (b2 - b1) % k
            if da != 0 and db != 0 and da != db:
                # Check if they actually interleave
                if (da < k // 2 + 1) != (db < k // 2 + 1):
                    crossings += 1
    return crossings


def spoke_deviation(seq_outer, seq_inner, adj, offset_steps):
    """
    Measure how far spokes deviate from radial.
    A radial spoke connects outer pos i to inner pos i (deviation = 0).
    Returns sum of squared angular deviations.
    """
    k = len(seq_outer)
    outer_set = set(seq_outer)
    inner_set = set(seq_inner)
    pos_outer = {v: i for i, v in enumerate(seq_outer)}
    pos_inner = {v: (i + offset_steps) % k for i, v in enumerate(seq_inner)}

    total_dev = 0
    n_spokes = 0
    for u in seq_outer:
        for v in adj[u - 1]:
            if v in inner_set:
                po = pos_outer[u]
                pi = pos_inner[v]
                # Angular deviation (minimum arc distance on circle)
                dev = min(abs(pi - po), k - abs(pi - po))
                total_dev += dev * dev
                n_spokes += 1

    return total_dev


def find_best_exponent(seq_outer, seq_inner, adj, k):
    """
    Buchheim & Hong Corollary 1: try all k possible angular offsets
    for the inner ring.

    Primary: minimize inter-orbit crossings.
    Secondary (tie-break): minimize spoke deviation from radial.
    """
    best_offset = 0
    best_score = (float('inf'), float('inf'))

    for offset in range(k):
        crossings = count_crossings_between_rings(seq_outer, seq_inner, adj, offset)
        deviation = spoke_deviation(seq_outer, seq_inner, adj, offset)
        score = (crossings, deviation)
        if score < best_score:
            best_score = score
            best_offset = offset

    return best_offset, best_score[0]


def layout_cylindrical(quotient, n, adj):
    """
    Concentric rings with global crossing minimization.

    1. Reorder orbits to follow quotient path (if applicable)
    2. Assign orbits to concentric rings
    3. Global angular offset optimization minimizing total crossings
    """
    orbits = list(quotient['orbits'])
    n_orbits = quotient['n_orbits']
    has_cycle = quotient.get('has_internal_cycle', [False] * n_orbits)

    # For 2-orbit GP pattern: put the cycle orbit on the outside
    if n_orbits == 2 and len(has_cycle) == 2:
        if has_cycle[1] and not has_cycle[0]:
            orbits = [orbits[1], orbits[0]]
            quotient = dict(quotient)
            oo = quotient.get('orbit_order', None)
            if oo and len(oo) == 2:
                quotient['orbit_order'] = [oo[1], oo[0]]

    # For 3+ orbits: reorder to follow quotient graph path
    elif n_orbits >= 3:
        path_order = find_path_ordering(quotient)
        if path_order is not None and path_order != list(range(n_orbits)):
            print(f"  Reordering orbits to follow quotient path: "
                  f"{[i+1 for i in path_order]}")
            orbits = [orbits[i] for i in path_order]
            quotient = dict(quotient)
            oo = quotient.get('orbit_order', None)
            if oo and len(oo) == n_orbits:
                quotient['orbit_order'] = [oo[i] for i in path_order]
            hc = quotient.get('has_internal_cycle', None)
            if hc and len(hc) == n_orbits:
                quotient['has_internal_cycle'] = [hc[i] for i in path_order]

    max_size = max(len(o) for o in orbits)

    # Get orbit sequences
    sequences = []
    for oi in range(n_orbits):
        sequences.append(get_orbit_sequence(quotient, oi, orbits[oi]))

    # Radii
    R_max = 4.5
    R_min = max(1.0, R_max * 0.35)
    radii = []
    for i in range(n_orbits):
        if n_orbits == 1:
            radii.append(R_max)
        else:
            radii.append(R_max - i * (R_max - R_min) / (n_orbits - 1))

    def build_coords(offsets):
        """Build coordinate dict from angular offsets (degrees per ring)."""
        coords = {}
        for oi in range(n_orbits):
            seq = sequences[oi]
            k = len(seq)
            r = radii[oi]
            for vi, v in enumerate(seq):
                theta = math.radians(90 + vi * (360.0 / k) + offsets[oi])
                coords[v] = (r * math.cos(theta), r * math.sin(theta))
        return coords

    # Global angular offset optimization
    # Ring 0 is fixed at offset 0; try all offset combinations for rings 1..n-1
    orbit_sizes = [len(sequences[i]) for i in range(n_orbits)]
    total_combos = 1
    for i in range(1, n_orbits):
        total_combos *= orbit_sizes[i]

    # Also try reflected (reversed) sequences for each ring
    # Reflection = traversing orbit in opposite direction under generator
    total_combos_with_mirror = total_combos * (2 ** (n_orbits - 1))

    if total_combos_with_mirror <= 20000:
        # Exhaustive global search with rotations AND reflections
        from itertools import product

        best_angular_offsets = [0.0] * n_orbits
        best_sequences = list(sequences)
        best_crossings = float('inf')
        best_edge_len = float('inf')

        offset_ranges = [range(orbit_sizes[i]) for i in range(1, n_orbits)]
        mirror_ranges = [[False, True]] * (n_orbits - 1)

        for mirrors in product(*mirror_ranges):
            # Build trial sequences with optional reversal
            trial_seqs = [sequences[0]]  # Ring 0 fixed
            for i, m in enumerate(mirrors):
                seq = sequences[i + 1]
                trial_seqs.append(seq[::-1] if m else seq)

            def build_trial_coords(offsets):
                coords = {}
                for oi in range(n_orbits):
                    seq = trial_seqs[oi]
                    k = len(seq)
                    r = radii[oi]
                    for vi, v in enumerate(seq):
                        theta = math.radians(90 + vi * (360.0 / k) + offsets[oi])
                        coords[v] = (r * math.cos(theta), r * math.sin(theta))
                return coords

            for combo in product(*offset_ranges):
                offsets = [0.0]
                for idx, c in enumerate(combo):
                    k_i = orbit_sizes[idx + 1]
                    offsets.append(c * (360.0 / k_i))

                trial_coords = build_trial_coords(offsets)
                crossings = count_all_crossings(trial_coords, adj, n)

                if crossings < best_crossings:
                    best_crossings = crossings
                    best_edge_len = sum(
                        math.hypot(trial_coords[v][0] - trial_coords[w][0],
                                   trial_coords[v][1] - trial_coords[w][1])
                        for v in range(1, n + 1)
                        for w in adj[v - 1] if v < w
                    )
                    best_angular_offsets = offsets[:]
                    best_sequences = list(trial_seqs)
                elif crossings == best_crossings:
                    edge_len = sum(
                        math.hypot(trial_coords[v][0] - trial_coords[w][0],
                                   trial_coords[v][1] - trial_coords[w][1])
                        for v in range(1, n + 1)
                        for w in adj[v - 1] if v < w
                    )
                    if edge_len < best_edge_len:
                        best_edge_len = edge_len
                        best_angular_offsets = offsets[:]
                        best_sequences = list(trial_seqs)

        print(f"  Global optimization: {total_combos_with_mirror} combos "
              f"(rotations + reflections) → {best_crossings} crossings")
        angular_offsets = best_angular_offsets
        sequences = best_sequences
    else:
        # Fall back to greedy pairwise for large search spaces
        angular_offsets = [0.0] * n_orbits
        step = 360.0 / max_size
        for i in range(1, n_orbits):
            k = len(sequences[i])
            best_offset, best_cross = find_best_exponent(
                sequences[i - 1], sequences[i], adj, k
            )
            angular_offsets[i] = angular_offsets[i - 1] + best_offset * step
            print(f"    Ring {i}: offset={best_offset} steps "
                  f"({best_offset * step:.0f}°), crossings={best_cross}")

    # Place vertices
    coords = build_coords(angular_offsets)
    total_cross = count_all_crossings(coords, adj, n)
    print(f"  Total crossings: {total_cross}")

    return coords


def compute_layout(quotient, n, adj, mode="auto"):
    """
    Layout modes:
      auto      - detect cyclic vs cylindrical from quotient structure
      cyclic    - Eades-Hong rotation (single circle, interleaved)
      rings     - cylindrical (concentric rings)
    """
    if mode == "auto":
        sym_type = detect_symmetry_type(quotient, adj)
        print(f"  Detected: {sym_type} symmetry")
    else:
        sym_type = mode

    if sym_type in ("cylindrical", "rings"):
        print(f"  Layout: concentric rings ({quotient['n_orbits']} rings)")
        coords = layout_cylindrical(quotient, n, adj)
    else:
        print(f"  Layout: cyclic (single circle)")
        coords = layout_cyclic(quotient, n, adj)
        total_cross = count_all_crossings(coords, adj, n)
        print(f"  Total crossings: {total_cross}")

    return coords


# ============================================================
# TikZ output
# ============================================================

def generate_tikz(coords, adj, n, orbits, graph_name, quotient_info,
                   aut_order=None, aut_structure=None):
    n_orbits = len(orbits)
    v_to_orbit = {}
    for oi, orbit in enumerate(orbits):
        for v in orbit:
            v_to_orbit[v] = oi

    lines = []
    lines.append(r"\documentclass[border=15pt]{standalone}")
    lines.append(r"\usepackage{amsmath, amssymb}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\usetikzlibrary{arrows.meta}")
    lines.append("")

    for i in range(min(n_orbits, len(COLORS))):
        lines.append(f"\\definecolor{{block{i}}}{{HTML}}{{{COLORS[i % len(COLORS)]}}}")
    lines.append("")

    lines.append(r"\begin{document}")
    lines.append(r"\begin{tikzpicture}[")
    lines.append(r"  vertex/.style={circle, draw, thick, minimum size=6.5mm, inner sep=0pt,")
    lines.append(r"                 font=\scriptsize\bfseries},")
    lines.append(r"  edge/.style={thick, gray!45},")
    lines.append(r"  title/.style={font=\Large\bfseries},")
    lines.append(r"  subtitle/.style={font=\small, text=gray},")
    lines.append(r"]")
    lines.append("")

    q_label = quotient_info.get('subgroup_structure', '?')
    lines.append(f"\\node[title] at (0, 6.5) {{{graph_name}}};")
    if aut_structure and aut_order:
        lines.append(f"\\node[subtitle] at (0, 5.8) "
                     f"{{$\\mathrm{{Aut}} = {aut_structure}$, "
                     f"$|\\mathrm{{Aut}}| = {aut_order}$}};")
        lines.append(f"\\node[subtitle] at (0, 5.2) "
                     f"{{$H = {q_label}$ $\\leq$ $\\mathrm{{Aut}}$, "
                     f"{n_orbits} orbits of {quotient_info['orbit_sizes']}}};")
    else:
        lines.append(f"\\node[subtitle] at (0, 5.8) "
                     f"{{Orbits of ${q_label}$ ({n_orbits} orbits of "
                     f"{quotient_info['orbit_sizes']})}};")
    lines.append("")

    # Vertices
    for v in range(1, n + 1):
        x, y = coords[v]
        oi = v_to_orbit[v]
        ci = oi % len(COLORS)
        lines.append(f"\\node[vertex, fill=block{ci}!30] "
                     f"(v{v}) at ({x:.3f}, {y:.3f}) {{{v}}};")

    lines.append("")

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
                    style = f"thick, block{oi_v % len(COLORS)}!60"
                else:
                    style = "edge"
                lines.append(f"\\draw[{style}] (v{edge[0]}) -- (v{edge[1]});")

    lines.append("")

    # Legend
    orbit_items = []
    for oi, orbit in enumerate(orbits):
        ci = oi % len(COLORS)
        verts = ",".join(str(v) for v in sorted(orbit))
        has_cycle = quotient_info.get('has_internal_cycle', [False] * n_orbits)
        cycle_note = " (cycle)" if oi < len(has_cycle) and has_cycle[oi] else ""
        orbit_items.append(
            f"\\textcolor{{block{ci}}}{{$\\blacksquare$}} "
            f"$O_{{{oi+1}}}$ = \\{{{verts}\\}}{cycle_note}"
        )

    legend = " \\\\\n  ".join(orbit_items)
    lines.append(f"\\node[draw, rounded corners, fill=gray!5, "
                 f"text width=12cm, align=left, font=\\small]")
    lines.append(f"  at (0, -6.5) {{")
    lines.append(f"  \\textbf{{Orbit decomposition}} "
                 f"(${q_label}$, order {quotient_info.get('subgroup_order','?')}) \\\\[4pt]")
    lines.append(f"  {legend}")
    lines.append(f"}};")
    lines.append("")
    lines.append(r"\end{tikzpicture}")
    lines.append(r"\end{document}")

    return "\n".join(lines)


# ============================================================
# Best layout selection
# ============================================================

def find_best_layout(results, n, adj, layout_mode="auto", verbose=True):
    """
    Find the best layout:
    1. Use top-scored quotient (preserves visual idiom)
    2. Apply improved algorithms (path reordering, global offset optimization)
    3. Try both auto and rings modes for that quotient
    4. Also try sibling quotients with same orbit structure for better offsets
    """
    scored = [(score_quotient(r), i, r) for i, r in enumerate(results)]
    scored.sort(reverse=True)

    top_q = scored[0][2]
    top_structure = (top_q['n_orbits'], tuple(sorted(top_q['orbit_sizes'])))

    modes = [layout_mode] if layout_mode != "auto" else ["auto"]

    best = None
    best_coords = None
    best_crossings = float('inf')

    # Try top-scored quotient and siblings with same orbit structure
    for score, idx, q in scored:
        structure = (q['n_orbits'], tuple(sorted(q['orbit_sizes'])))
        if structure != top_structure:
            continue
        for mode in modes:
            try:
                coords = compute_layout(q, n, adj, mode=mode)
                crossings = count_all_crossings(coords, adj, n)
            except Exception as e:
                if verbose:
                    print(f"    → ERROR: {e}")
                continue
            if verbose:
                print(f"    → {crossings} crossings")
            if crossings < best_crossings or \
               (crossings == best_crossings and score > score_quotient(best)):
                best_crossings = crossings
                best_coords = coords
                best = q

    if verbose:
        print(f"\n  Selected: {best.get('label','?')} with {best_crossings} crossings")

    return best_coords, best


# ============================================================
# Main
# ============================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Draw symmetric graphs as PDF")
    parser.add_argument("graph", help="Graph name (petersen, dodecahedron, heawood, ...)")
    parser.add_argument("--max-deg", type=int, default=4, help="Max quotient degree")
    parser.add_argument("--layout", choices=["auto", "cyclic", "rings"],
                        default="auto",
                        help="auto (detect), cyclic (single circle), rings (concentric)")
    parser.add_argument("--no-open", action="store_true")
    args = parser.parse_args()

    print(f"Running GAP for {args.graph}...")
    data = run_gap(args.graph, args.max_deg)

    graph_name = data['graph']
    n = data['n']
    adj = data['adj']
    results = data['results']

    print(f"Graph: {graph_name} ({n} vertices)")

    if not results:
        print("No suitable quotients found!")
        sys.exit(1)

    coords, best = find_best_layout(results, n, adj, layout_mode=args.layout)
    tikz = generate_tikz(coords, adj, n, best['orbits'], graph_name, best,
                         aut_order=data.get('aut_order'),
                         aut_structure=data.get('aut_structure'))

    tex_file = os.path.join(GRAPHVIZ_DIR, f"auto-{args.graph}.tex")
    pdf_file = os.path.join(GRAPHVIZ_DIR, f"auto-{args.graph}.pdf")

    with open(tex_file, 'w') as f:
        f.write(tikz)

    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", f"auto-{args.graph}.tex"],
        capture_output=True, text=True, cwd=GRAPHVIZ_DIR
    )

    if result.returncode != 0:
        print("LaTeX compilation failed!")
        print(result.stdout[-500:])
        sys.exit(1)

    for ext in ['aux', 'log']:
        try:
            os.remove(os.path.join(GRAPHVIZ_DIR, f"auto-{args.graph}.{ext}"))
        except FileNotFoundError:
            pass

    print(f"  {pdf_file}")

    if not args.no_open:
        subprocess.run(["open", pdf_file])


if __name__ == "__main__":
    main()
