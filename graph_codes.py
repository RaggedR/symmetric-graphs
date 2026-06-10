#!/usr/bin/env python3
"""
Error-correcting codes from symmetric graphs.

Constructs three families of classical codes from cubic arc-transitive graphs:

1. Cycle codes: null space of the vertex-edge incidence matrix over F_2
   Parameters: [m, m-n+1, g] where m = edges, n = vertices, g = girth

2. Incidence codes: for bipartite graphs, the biadjacency matrix gives
   point codes (row space) and dual codes (null space) over F_2

3. Quotient codes: project cycle codes through orbit quotient maps,
   producing shorter codes that inherit symmetry

Usage:
    python3 graph_codes.py petersen          # Single graph
    python3 graph_codes.py --all             # Full census table
    python3 graph_codes.py heawood --verbose # Detailed output with generator matrix
"""

import sys
import os
import argparse
from copy import deepcopy

# Reuse GAP interface from the drawing pipeline
GRAPHVIZ_DIR = os.path.expanduser("~/git/graphviz")
sys.path.insert(0, GRAPHVIZ_DIR)
from draw_graph import run_gap


# ============================================================
# F_2 linear algebra
# ============================================================

def mat_mod2(rows):
    """Ensure all entries are mod 2."""
    return [[x % 2 for x in row] for row in rows]


def rref_f2(matrix):
    """
    Row-reduce a matrix over F_2 (in-place).
    Returns (rref_matrix, pivot_columns, rank).
    """
    if not matrix:
        return [], [], 0

    M = [row[:] for row in matrix]  # copy
    nrows = len(M)
    ncols = len(M[0])
    pivots = []
    row_idx = 0

    for col in range(ncols):
        # Find pivot in this column
        pivot = None
        for r in range(row_idx, nrows):
            if M[r][col] % 2 == 1:
                pivot = r
                break
        if pivot is None:
            continue

        # Swap to row_idx
        M[row_idx], M[pivot] = M[pivot], M[row_idx]
        pivots.append(col)

        # Eliminate all other 1s in this column
        for r in range(nrows):
            if r != row_idx and M[r][col] % 2 == 1:
                M[r] = [(M[r][j] + M[row_idx][j]) % 2 for j in range(ncols)]

        row_idx += 1

    return M, pivots, len(pivots)


def null_space_f2(matrix):
    """
    Compute the null space of a matrix over F_2.
    Returns a list of basis vectors (as lists of 0/1).
    """
    if not matrix:
        return []

    nrows = len(matrix)
    ncols = len(matrix[0])

    _, pivots, rank = rref_f2(matrix)
    pivot_set = set(pivots)
    free_vars = [j for j in range(ncols) if j not in pivot_set]

    # For each free variable, construct a null space vector
    # First, get the RREF form
    M, _, _ = rref_f2(matrix)

    basis = []
    for fv in free_vars:
        vec = [0] * ncols
        vec[fv] = 1
        # For each pivot column, solve for the pivot variable
        for i, pc in enumerate(pivots):
            if i < len(M):
                vec[pc] = M[i][fv] % 2  # since we're in F_2, this is the value
        basis.append(vec)

    return basis


def row_space_f2(matrix):
    """
    Compute a basis for the row space over F_2.
    Returns the non-zero rows of the RREF.
    """
    M, pivots, rank = rref_f2(matrix)
    return [M[i] for i in range(rank)]


def hamming_weight(vec):
    """Hamming weight of a binary vector."""
    return sum(v % 2 for v in vec)


def minimum_distance_bruteforce(basis):
    """
    Compute minimum distance by exhaustive search over all non-zero codewords.
    Only feasible for small codes (dimension <= ~20).
    """
    k = len(basis)
    if k == 0:
        return 0
    n = len(basis[0])

    if k > 20:
        return None  # too large for brute force

    min_wt = n + 1
    # Iterate over all 2^k - 1 non-zero linear combinations
    for mask in range(1, 1 << k):
        codeword = [0] * n
        for i in range(k):
            if mask & (1 << i):
                codeword = [(codeword[j] + basis[i][j]) % 2 for j in range(n)]
        wt = hamming_weight(codeword)
        if wt < min_wt:
            min_wt = wt
    return min_wt


# ============================================================
# Graph utilities
# ============================================================

def edges_from_adj(adj, n):
    """Extract edge list from adjacency list. Returns list of (u, v) with u < v."""
    edges = []
    for v in range(n):
        for w in adj[v]:
            w_idx = w - 1  # 1-indexed to 0-indexed
            if v < w_idx:
                edges.append((v, w_idx))
    return edges


def incidence_matrix_f2(adj, n):
    """
    Build the vertex-edge incidence matrix B over F_2.
    B[v][e] = 1 iff vertex v is an endpoint of edge e.
    Returns (matrix, edge_list).
    """
    edges = edges_from_adj(adj, n)
    m = len(edges)
    B = [[0] * m for _ in range(n)]
    for e_idx, (u, v) in enumerate(edges):
        B[u][e_idx] = 1
        B[v][e_idx] = 1
    return B, edges


def is_bipartite(adj, n):
    """
    Check if graph is bipartite. Returns (True, part_A, part_B) or (False, None, None).
    Uses BFS 2-coloring.
    """
    color = [-1] * n
    color[0] = 0
    queue = [0]
    qi = 0
    while qi < len(queue):
        u = queue[qi]
        qi += 1
        for w in adj[u]:
            w_idx = w - 1
            if color[w_idx] == -1:
                color[w_idx] = 1 - color[u]
                queue.append(w_idx)
            elif color[w_idx] == color[u]:
                return False, None, None
    part_a = [i for i in range(n) if color[i] == 0]
    part_b = [i for i in range(n) if color[i] == 1]
    return True, part_a, part_b


def biadjacency_matrix_f2(adj, n, part_a, part_b):
    """
    Build the biadjacency matrix M for a bipartite graph.
    M[i][j] = 1 iff vertex part_a[i] is adjacent to vertex part_b[j].
    """
    b_index = {v: j for j, v in enumerate(part_b)}
    M = [[0] * len(part_b) for _ in range(len(part_a))]
    for i, u in enumerate(part_a):
        for w in adj[u]:
            w_idx = w - 1
            if w_idx in b_index:
                M[i][b_index[w_idx]] = 1
    return M


def compute_girth(adj, n):
    """
    Compute girth (shortest cycle length) by BFS from each vertex.
    """
    girth = float('inf')
    for start in range(n):
        dist = [-1] * n
        dist[start] = 0
        parent = [-1] * n
        queue = [start]
        qi = 0
        while qi < len(queue):
            u = queue[qi]
            qi += 1
            for w in adj[u]:
                w_idx = w - 1
                if dist[w_idx] == -1:
                    dist[w_idx] = dist[u] + 1
                    parent[w_idx] = u
                    queue.append(w_idx)
                elif parent[u] != w_idx and parent[w_idx] != u:
                    # Found a cycle
                    cycle_len = dist[u] + dist[w_idx] + 1
                    girth = min(girth, cycle_len)
            if dist[u] >= girth // 2:
                break  # can't find shorter cycle from here
    return girth if girth < float('inf') else 0


# ============================================================
# Code constructions
# ============================================================

def cycle_code(adj, n, verbose=False):
    """
    Construct the cycle code of a graph.

    The cycle space over F_2 is the null space of the vertex-edge incidence
    matrix. For a connected graph with n vertices and m edges:
      - Length: m
      - Dimension: m - n + 1  (by Euler's formula)
      - Minimum distance: girth g

    Returns dict with code parameters and optionally the generator matrix.
    """
    B, edges = incidence_matrix_f2(adj, n)
    m = len(edges)

    # Cycle space = null space of B over F_2
    basis = null_space_f2(B)
    k = len(basis)

    # Expected: k = m - n + 1 for connected graphs
    expected_k = m - n + 1

    # Girth = minimum distance
    girth = compute_girth(adj, n)

    # Compute exact minimum distance for small codes
    exact_d = None
    if k <= 20:
        exact_d = minimum_distance_bruteforce(basis)

    result = {
        'type': 'cycle',
        'n': m,       # code length = number of edges
        'k': k,       # dimension
        'd_girth': girth,  # minimum distance from girth
        'd_exact': exact_d,
        'edges': edges,
    }

    if verbose:
        result['generator'] = basis
        print(f"  Cycle code: [{m}, {k}, {girth}]")
        print(f"    Vertices: {n}, Edges: {m}")
        print(f"    Expected dimension: {expected_k}, Actual: {k}")
        if exact_d is not None:
            print(f"    Exact minimum distance: {exact_d}")
            if exact_d != girth:
                print(f"    NOTE: exact d={exact_d} differs from girth={girth}")

    return result


def incidence_code(adj, n, verbose=False):
    """
    Construct incidence codes from a bipartite graph.

    For a bipartite graph with parts A, B and biadjacency matrix M:
      - Point code C_p = row_space(M) over F_2: length |B|, dim = rank(M)
      - Dual code C_p^perp = null_space(M^T) over F_2

    Returns dict with both codes, or None if graph is not bipartite.
    """
    bip, part_a, part_b = is_bipartite(adj, n)
    if not bip:
        if verbose:
            print("  Not bipartite — no incidence code")
        return None

    M = biadjacency_matrix_f2(adj, n, part_a, part_b)
    na, nb = len(part_a), len(part_b)

    # Point code: row space of M (length nb)
    point_basis = row_space_f2(M)
    point_k = len(point_basis)

    # Dual code: null space of M^T (length na)
    MT = [[M[i][j] for i in range(na)] for j in range(nb)]
    dual_basis = null_space_f2(MT)
    dual_k = len(dual_basis)

    # Also compute the "full" incidence code: null space of M (length nb)
    null_basis = null_space_f2(M)
    null_k = len(null_basis)

    # Minimum distances
    point_d = minimum_distance_bruteforce(point_basis) if point_k <= 20 else None
    dual_d = minimum_distance_bruteforce(dual_basis) if dual_k <= 20 else None
    null_d = minimum_distance_bruteforce(null_basis) if null_k <= 20 else None

    result = {
        'type': 'incidence',
        'bipartition': (na, nb),
        'point_code': {'n': nb, 'k': point_k, 'd': point_d},
        'dual_code': {'n': na, 'k': dual_k, 'd': dual_d},
        'null_code': {'n': nb, 'k': null_k, 'd': null_d},
    }

    if verbose:
        print(f"  Bipartite: ({na}, {nb})")
        print(f"  Point code (row space of M): [{nb}, {point_k}, {point_d}]")
        print(f"  Null code (null space of M): [{nb}, {null_k}, {null_d}]")
        print(f"  Dual code (null space of M^T): [{na}, {dual_k}, {dual_d}]")

    if verbose:
        result['point_generator'] = point_basis
        result['null_generator'] = null_basis

    return result


def quotient_code(adj, n, orbits, verbose=False):
    """
    Construct the quotient code by projecting the cycle code through the orbit map.

    Given orbits O_1, ..., O_k, the projection pi: G -> G/H collapses each orbit
    to a single vertex. The image of the cycle code under pi is the cycle code of
    the quotient graph.

    Returns dict with quotient code parameters.
    """
    # Build quotient adjacency
    n_orbs = len(orbits)
    orbit_of = {}
    for oi, orb in enumerate(orbits):
        for v in orb:
            orbit_of[v - 1] = oi  # convert to 0-indexed

    # Quotient adjacency (may have multi-edges, but we reduce to simple)
    quot_adj = [set() for _ in range(n_orbs)]
    self_loops = [False] * n_orbs
    for v in range(n):
        ov = orbit_of[v]
        for w in adj[v]:
            ow = orbit_of[w - 1]
            if ov == ow:
                self_loops[ov] = True
            else:
                quot_adj[ov].add(ow + 1)  # 1-indexed for consistency

    # Convert to list format matching draw_graph conventions
    quot_adj_list = [sorted(list(s)) for s in quot_adj]

    # Compute cycle code of quotient
    # Number of edges in quotient (simple graph)
    quot_edges = edges_from_adj(quot_adj_list, n_orbs)
    m_quot = len(quot_edges)

    if m_quot == 0:
        if verbose:
            print("  Quotient has no edges")
        return None

    B_quot, _ = incidence_matrix_f2(quot_adj_list, n_orbs)
    basis_quot = null_space_f2(B_quot)
    k_quot = len(basis_quot)

    girth_quot = compute_girth(quot_adj_list, n_orbs)

    d_quot = None
    if k_quot <= 20 and k_quot > 0:
        d_quot = minimum_distance_bruteforce(basis_quot)

    result = {
        'type': 'quotient_cycle',
        'n_orbits': n_orbs,
        'orbit_sizes': [len(o) for o in orbits],
        'n': m_quot,
        'k': k_quot,
        'd_girth': girth_quot,
        'd_exact': d_quot,
        'self_loops': sum(self_loops),
    }

    if verbose:
        print(f"  Quotient: {n_orbs} orbits, sizes {[len(o) for o in orbits]}")
        print(f"  Quotient edges: {m_quot}")
        print(f"  Quotient cycle code: [{m_quot}, {k_quot}, {girth_quot}]")
        if d_quot is not None and d_quot != girth_quot:
            print(f"    Exact minimum distance: {d_quot}")
        if sum(self_loops) > 0:
            print(f"    {sum(self_loops)} orbits with internal edges (self-loops in quotient)")

    return result


# ============================================================
# Analysis pipeline
# ============================================================

def analyse_graph(graph_name, verbose=False, max_deg=4):
    """
    Run the full code analysis pipeline on a graph.
    """
    print(f"\n{'='*60}")
    print(f"  {graph_name}")
    print(f"{'='*60}")

    try:
        data = run_gap(graph_name, max_deg)
    except Exception as e:
        print(f"  GAP error: {e}")
        return None

    n = data['n']
    adj = data['adj']
    results = data.get('results', [])
    aut_order = data.get('aut_order', '?')
    aut_structure = data.get('aut_structure', '?')

    print(f"  |V| = {n}, |Aut| = {aut_order} ({aut_structure})")
    print()

    # 1. Cycle code
    print("  --- Cycle Code ---")
    cc = cycle_code(adj, n, verbose=verbose)
    d_str = cc['d_exact'] if cc['d_exact'] is not None else cc['d_girth']
    print(f"  C_cycle = [{cc['n']}, {cc['k']}, {d_str}]")
    print()

    # 2. Incidence code (bipartite graphs only)
    print("  --- Incidence Code ---")
    ic = incidence_code(adj, n, verbose=verbose)
    if ic is None:
        print("  (not bipartite)")
    print()

    # 3. Quotient codes from best quotients
    print("  --- Quotient Codes ---")
    qc_list = []
    if results:
        # Score and sort quotients
        from draw_graph import score_quotient
        scored = [(score_quotient(r, n), i, r) for i, r in enumerate(results)]
        scored.sort(reverse=True)

        seen_structures = set()
        for score, idx, q in scored[:5]:
            orbit_key = tuple(sorted(q['orbit_sizes']))
            if orbit_key in seen_structures:
                continue
            seen_structures.add(orbit_key)

            label = q.get('subgroup_structure', '?')
            print(f"\n  H = {label} (order {q.get('subgroup_order', '?')})")
            qc = quotient_code(adj, n, q['orbits'], verbose=verbose)
            if qc:
                qc['label'] = label
                qc_list.append(qc)
    else:
        print("  No quotients available")

    return {
        'name': graph_name,
        'n': n,
        'aut_order': aut_order,
        'aut_structure': aut_structure,
        'cycle_code': cc,
        'incidence_code': ic,
        'quotient_codes': qc_list,
    }


def print_summary_table(all_results):
    """Print a summary table of all codes found."""
    print(f"\n{'='*80}")
    print(f"  SUMMARY: Error-Correcting Codes from Symmetric Graphs")
    print(f"{'='*80}")

    # Cycle codes
    print(f"\n  {'Graph':<20} {'|V|':>4} {'|Aut|':>8}  {'Cycle Code':>18}  {'Rate':>6}  {'Bipartite'}")
    print(f"  {'-'*20} {'-'*4} {'-'*8}  {'-'*18}  {'-'*6}  {'-'*9}")

    for r in all_results:
        if r is None:
            continue
        cc = r['cycle_code']
        d = cc['d_exact'] if cc['d_exact'] is not None else cc['d_girth']
        code_str = f"[{cc['n']}, {cc['k']}, {d}]"
        rate = cc['k'] / cc['n'] if cc['n'] > 0 else 0
        bip = "yes" if r['incidence_code'] is not None else "no"
        print(f"  {r['name']:<20} {r['n']:>4} {r['aut_order']:>8}  {code_str:>18}  {rate:>6.3f}  {bip}")

    # Incidence codes
    inc_results = [r for r in all_results if r and r['incidence_code'] is not None]
    if inc_results:
        print(f"\n  Incidence Codes (bipartite graphs only):")
        print(f"  {'Graph':<20} {'Partition':>10}  {'Point Code':>18}  {'Null Code':>18}")
        print(f"  {'-'*20} {'-'*10}  {'-'*18}  {'-'*18}")
        for r in inc_results:
            ic = r['incidence_code']
            pc = ic['point_code']
            nc = ic['null_code']
            part_str = f"({ic['bipartition'][0]},{ic['bipartition'][1]})"
            pc_str = f"[{pc['n']}, {pc['k']}, {pc['d']}]"
            nc_str = f"[{nc['n']}, {nc['k']}, {nc['d']}]"
            print(f"  {r['name']:<20} {part_str:>10}  {pc_str:>18}  {nc_str:>18}")

    # Quotient codes
    quot_results = [r for r in all_results if r and r['quotient_codes']]
    if quot_results:
        print(f"\n  Quotient Codes (cycle code of G/H):")
        print(f"  {'Graph':<20} {'H':>12} {'Orbits':>8}  {'Code':>18}")
        print(f"  {'-'*20} {'-'*12} {'-'*8}  {'-'*18}")
        for r in quot_results:
            for qc in r['quotient_codes']:
                d = qc['d_exact'] if qc['d_exact'] is not None else qc['d_girth']
                if qc['k'] > 0:
                    code_str = f"[{qc['n']}, {qc['k']}, {d}]"
                else:
                    code_str = "(trivial)"
                print(f"  {r['name']:<20} {qc['label']:>12} {qc['n_orbits']:>8}  {code_str:>18}")

    print()


# ============================================================
# Main
# ============================================================

# Graphs from the Foster census (and beyond) available in symmetric-layout.g
FOSTER_GRAPHS = [
    'k4', 'k33', 'cube', 'petersen', 'heawood', 'moebiuskantor',
    'pappus', 'dodecahedron', 'desargues', 'nauru', 'f26a',
    'coxeter', 'tuttecoxeter',
]


def main():
    parser = argparse.ArgumentParser(
        description="Construct error-correcting codes from symmetric graphs"
    )
    parser.add_argument("graph", nargs='?', default=None,
                        help="Graph name (petersen, heawood, ...)")
    parser.add_argument("--all", action="store_true",
                        help="Analyse all Foster census graphs")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show detailed output including generator matrices")
    parser.add_argument("--max-deg", type=int, default=4,
                        help="Max quotient degree for GAP search")
    args = parser.parse_args()

    if args.all:
        graphs = FOSTER_GRAPHS
    elif args.graph:
        graphs = [args.graph]
    else:
        parser.print_help()
        sys.exit(1)

    all_results = []
    for g in graphs:
        r = analyse_graph(g, verbose=args.verbose, max_deg=args.max_deg)
        all_results.append(r)

    if len(all_results) > 1:
        print_summary_table(all_results)


if __name__ == "__main__":
    main()
