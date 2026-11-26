#!/usr/bin/env python3
"""
Generalized graph drawing for non-cubic symmetric graphs.
Bypasses GAP and accepts adjacency lists + orbit structure directly.

Usage:
    python3 draw_noncubic.py --graph e8_mandala
    python3 draw_noncubic.py --graph schlafli
"""

import sys
import os
import argparse
import subprocess
import tempfile

# Import the drawing functions from draw_graph.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from draw_graph import (
    score_quotient, detect_symmetry_type, compute_layout,
    generate_tikz, COLORS
)

def create_mock_quotient(orbits, orbit_adj, edge_counts, has_internal_cycle=None):
    """
    Create a quotient structure compatible with draw_graph.py

    Args:
        orbits: list of lists, e.g. [[0,1,2], [3,4,5], ...]
        orbit_adj: adjacency list for orbit quotient graph (1-indexed)
        edge_counts: matrix of edge counts between orbits
        has_internal_cycle: list of booleans, one per orbit
    """
    n_orbits = len(orbits)
    orbit_sizes = [len(o) for o in orbits]

    # Convert to 1-indexed orbits (for compatibility with GAP output)
    orbits_1indexed = [[v+1 for v in orb] for orb in orbits]

    # Compute valencies
    valencies = [len(adj) for adj in orbit_adj]
    max_valency = max(valencies) if valencies else 0

    # Check if quotient is a cycle
    is_cycle = all(v == 2 for v in valencies)

    # Create orbit_order (the vertices in order within each orbit)
    orbit_order = orbits_1indexed.copy()

    if has_internal_cycle is None:
        has_internal_cycle = [False] * n_orbits

    return {
        'n_orbits': n_orbits,
        'orbits': orbits_1indexed,
        'orbit_sizes': orbit_sizes,
        'quotient_adj': orbit_adj,
        'quotient_valencies': valencies,
        'max_valency': max_valency,
        'edge_counts': edge_counts,
        'internal_counts': [0] * n_orbits,  # TODO: compute properly
        'is_cycle': is_cycle,
        'has_internal_cycle': has_internal_cycle,
        'orbit_order': orbit_order,
        'offsets': [[0]*n_orbits for _ in range(n_orbits)]  # TODO: compute
    }


def build_e8_mandala():
    """
    Build the E8 mandala graph (ip=+1, 56-regular, 240 vertices).
    Uses pre-computed root structure with 8 rings of 30.
    """
    import numpy as np
    import itertools

    # Build 240 E8 roots
    roots = []
    # Type A: (±1,±1,0^6)
    for i, j in itertools.combinations(range(8), 2):
        for si, sj in itertools.product([1, -1], repeat=2):
            v = [0]*8; v[i] = si; v[j] = sj
            roots.append(v)
    # Type B: (±1/2)^8 with even number of minus signs
    for signs in itertools.product([1, -1], repeat=8):
        if signs.count(-1) % 2 == 0:
            roots.append([s*0.5 for s in signs])

    R = np.array(roots, dtype=float)

    # Build adjacency list (ip = +1 means 60-degree angle)
    adj = [[] for _ in range(240)]
    for i in range(240):
        for j in range(i+1, 240):
            ip = np.dot(R[i], R[j])
            if abs(ip - 1.0) < 1e-6:  # inner product = +1
                adj[i].append(j)
                adj[j].append(i)

    # Verify degree
    degrees = [len(a) for a in adj]
    print(f"E8 mandala: {len(R)} vertices, degrees: min={min(degrees)}, max={max(degrees)}")

    # Create 8 orbits of 30 vertices each based on Coxeter plane projection
    # For now, simple grouping by radius
    e = np.eye(8)
    a1 = 0.5*(e[0] - e[1] - e[2] - e[3] - e[4] - e[5] - e[6] + e[7])
    a2 = e[0] + e[1]
    a3 = e[1] - e[0]
    a4 = e[2] - e[1]
    a5 = e[3] - e[2]
    a6 = e[4] - e[3]
    a7 = e[5] - e[4]
    a8 = e[6] - e[5]
    simple = [a1, a2, a3, a4, a5, a6, a7, a8]

    def refl(alpha):
        return np.eye(8) - np.outer(alpha, alpha)

    C = np.eye(8)
    for a in simple:
        C = C @ refl(a)

    w, V = np.linalg.eig(C)
    theta = 2*np.pi/30
    idx = np.argmin(np.abs(np.angle(w) - theta))
    v = V[:, idx]
    u1, u2 = np.real(v), np.imag(v)
    u1 = u1/np.linalg.norm(u1)
    u2 = u2 - (u2@u1)*u1
    u2 = u2/np.linalg.norm(u2)

    P = np.column_stack([R @ u1, R @ u2])
    radii = np.round(np.hypot(P[:,0], P[:,1]), 4)

    # Group by radius
    radius_values = sorted(set(radii))
    orbits = []
    for r in radius_values:
        orbit = [i for i, rad in enumerate(radii) if rad == r]
        # Sort by angle for proper ordering
        angles = np.arctan2(P[orbit, 1], P[orbit, 0])
        orbit = [orbit[i] for i in np.argsort(angles)]
        orbits.append(orbit)

    print(f"Created {len(orbits)} orbits with sizes: {[len(o) for o in orbits]}")

    # Build orbit quotient
    n_orbits = len(orbits)
    orbit_map = {}
    for i, orb in enumerate(orbits):
        for v in orb:
            orbit_map[v] = i

    orbit_adj = [set() for _ in range(n_orbits)]
    edge_counts = [[0]*n_orbits for _ in range(n_orbits)]

    for i in range(240):
        orb_i = orbit_map[i]
        for j in adj[i]:
            orb_j = orbit_map[j]
            if orb_i != orb_j:
                orbit_adj[orb_i].add(orb_j)
                edge_counts[orb_i][orb_j] += 1

    orbit_adj_list = [sorted(list(s)) for s in orbit_adj]
    # Convert to 1-indexed
    orbit_adj_1indexed = [[x+1 for x in adj] for adj in orbit_adj_list]

    # All orbits form 30-fold cycles
    has_internal_cycle = [True] * n_orbits

    return {
        'name': 'E8-Mandala',
        'n': 240,
        'adj': [[v+1 for v in a] for a in adj],  # 1-indexed
        'quotient': create_mock_quotient(orbits, orbit_adj_1indexed,
                                         edge_counts, has_internal_cycle)
    }


def build_schlafli():
    """
    Build the Schläfli graph SRG(27,10,1,5).
    27 lines on cubic surface, E6 representation.
    """
    # TODO: Import from E6 structure or hardcode adjacency
    # For now, placeholder
    print("Schläfli graph construction not yet implemented")
    return None


GRAPH_BUILDERS = {
    'e8_mandala': build_e8_mandala,
    'schlafli': build_schlafli,
}


def main():
    parser = argparse.ArgumentParser(description='Draw non-cubic symmetric graphs')
    parser.add_argument('--graph', required=True, choices=GRAPH_BUILDERS.keys(),
                       help='Graph to draw')
    parser.add_argument('--layout', choices=['rings', 'cyclic', 'auto'], default='auto',
                       help='Layout mode')
    parser.add_argument('--output', help='Output directory',
                       default='individual_graphs')

    args = parser.parse_args()

    print(f"Building {args.graph}...")
    builder = GRAPH_BUILDERS[args.graph]
    data = builder()

    if data is None:
        print(f"Failed to build {args.graph}")
        return 1

    # Layout
    quotient = data['quotient']
    n = data['n']
    adj = data['adj']

    mode = args.layout if args.layout != 'auto' else 'auto'
    coords = compute_layout(quotient, n, adj, mode=mode)

    # Extract orbits from quotient
    orbits = quotient['orbits']

    # Write TikZ
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    base_name = data['name'].lower().replace(' ', '-').replace('_', '-')
    tex_path = os.path.join(output_dir, f"auto-{base_name}.tex")

    quotient_info = {
        'label': 'Direct construction',
        'n_orbits': quotient['n_orbits'],
        'orbit_sizes': quotient['orbit_sizes'],
        'max_valency': quotient['max_valency']
    }

    tikz_code = generate_tikz(
        coords, adj, n, orbits,
        data['name'], quotient_info,
        aut_structure='Unknown',
        aut_order='Unknown'
    )

    with open(tex_path, 'w') as f:
        f.write(tikz_code)

    print(f"Wrote {tex_path}")

    # Compile
    try:
        result = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', tex_path],
            cwd=output_dir,
            capture_output=True,
            timeout=30
        )
        pdf_path = tex_path.replace('.tex', '.pdf')
        if os.path.exists(pdf_path):
            print(f"Generated {pdf_path}")
        else:
            print(f"pdflatex failed, see {tex_path}")
    except Exception as e:
        print(f"Could not compile PDF: {e}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
