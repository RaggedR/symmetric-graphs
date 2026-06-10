#!/usr/bin/env python3
"""Extract the cyclic quotient orbit data used for each named graph layout."""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/git/graphviz"))
from draw_graph import run_gap, score_quotient, detect_symmetry_type

GRAPHS = [
    ("heawood", 10),
    ("moebiuskantor", 10),
    ("dodecahedron", 10),
    ("psl2_7_c3", 4),      # Klein
    ("tuttecage12", 4),     # Tutte 12-cage
    ("psl2_13", 4),         # Zhou
    ("gh22", 6),            # Langer
]

for gap_name, max_deg in GRAPHS:
    print(f"\n{'='*60}")
    print(f"  {gap_name}")
    print(f"{'='*60}")

    data = run_gap(gap_name, max_deg)
    n = data["n"]
    results = data["results"]
    aut_order = data.get("aut_order", "?")

    print(f"  {n} vertices, |Aut| = {aut_order}")
    print(f"  {len(results)} quotient(s) found")

    if not results:
        continue

    # Score and pick the best quotient (same logic as find_best_layout)
    best = max(results, key=lambda r: score_quotient(r, n))

    orbits = best["orbits"]  # list of lists: each sublist is an orbit
    orbit_sizes = best["orbit_sizes"]
    n_orbits = len(orbits)
    is_cycle = best.get("is_cycle", False)
    sym_type = detect_symmetry_type(best, n)

    print(f"  Best quotient: {n_orbits} orbits, sizes {orbit_sizes}")
    print(f"  Is cycle: {is_cycle}, layout type: {sym_type}")

    # Build the orbit map: vertex (1-indexed) → orbit index (0-indexed)
    orbit_map = [0] * (n + 1)  # 1-indexed
    for oi, orbit in enumerate(orbits):
        for v in orbit:
            orbit_map[v] = oi

    # Print as 0-indexed array for Lean (Fin n → Fin n_orbits)
    lean_map = [orbit_map[v+1] for v in range(n)]
    print(f"  Orbit map (0-indexed, for Lean Fin {n} → Fin {n_orbits}):")

    # Print in rows of 10
    for start in range(0, n, 10):
        chunk = lean_map[start:start+10]
        print(f"    {chunk}")

    # Print orbits
    for oi, orbit in enumerate(orbits):
        orbit_0 = sorted([v-1 for v in orbit])  # 0-indexed
        print(f"  Orbit {oi} (size {len(orbit)}): {orbit_0}")

    # Quotient graph adjacency
    q_adj = best.get("quotient_adj", [])
    print(f"  Quotient adjacency: {q_adj}")

    # Check cyclic group order: orbit sizes should all be equal for cyclic quotient
    if len(set(orbit_sizes)) == 1:
        cyclic_order = orbit_sizes[0]
        print(f"  → Cyclic group: Z_{cyclic_order}")
    else:
        print(f"  → Mixed orbit sizes (not a free cyclic action)")
