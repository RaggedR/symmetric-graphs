"""
foster_census_spectrum.py

Computes adjacency spectrum, Laplacian spectrum, Fiedler eigenvalue,
spectral gap, Ramanujan status, and related invariants for all Foster census
cubic arc-transitive graphs up to 30 vertices.
"""

import numpy as np
from collections import deque
import json

# ---------------------------------------------------------------------------
# Graph constructors
# ---------------------------------------------------------------------------

def gp_graph(n, k):
    """Generalized Petersen graph GP(n, k).
    Outer cycle: 0..n-1, inner star: n..2n-1.
    Spokes: i — n+i. Outer cycle: i — (i+1)%n. Inner star: n+i — n+(i+k)%n.
    """
    adj = [set() for _ in range(2 * n)]
    for i in range(n):
        ip = (i + 1) % n
        adj[i].add(ip); adj[ip].add(i)
        adj[i].add(n + i); adj[n + i].add(i)
        j = n + (i + k) % n
        adj[n + i].add(j); adj[j].add(n + i)
    return [sorted(s) for s in adj]


def lcf_graph(n, shifts):
    """LCF graph: Hamiltonian cycle 0..n-1 plus chords from repeating shift list."""
    adj = [set() for _ in range(n)]
    period = len(shifts)
    for i in range(n):
        j = (i + 1) % n
        adj[i].add(j); adj[j].add(i)
        s = shifts[i % period]
        j = (i + s) % n
        adj[i].add(j); adj[j].add(i)
    return [sorted(s) for s in adj]


def complete_graph(n):
    return [list(range(j for j in range(n) if j != i)) for i in range(n)]


def k4():
    return [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


def k33():
    """K_{3,3}: parts {0,1,2} and {3,4,5}."""
    adj = []
    for i in range(3):
        adj.append([3, 4, 5])
    for i in range(3):
        adj.append([0, 1, 2])
    return adj


def dodecahedron():
    """Explicit adjacency (1-indexed in comments, 0-indexed here)."""
    adj_1 = [
        [2, 5, 6], [1, 3, 8], [2, 4, 10], [3, 5, 12], [1, 4, 14],
        [1, 7, 15], [6, 8, 16], [2, 7, 9], [8, 10, 17], [3, 9, 11],
        [10, 12, 18], [4, 11, 13], [12, 14, 19], [5, 13, 15], [6, 14, 20],
        [7, 17, 20], [9, 16, 18], [11, 17, 19], [13, 18, 20], [15, 16, 19]
    ]
    return [[v - 1 for v in row] for row in adj_1]


def coxeter():
    """Coxeter graph: hub 0-6, three rings of 7 with steps 1, 2, 3."""
    adj = [set() for _ in range(28)]
    steps = [1, 2, 3]
    for j in range(7):
        for ring_idx, step in enumerate(steps):
            hub = j
            ring_v = 7 * (ring_idx + 1) + j
            adj[hub].add(ring_v); adj[ring_v].add(hub)
            next_in_ring = 7 * (ring_idx + 1) + (j + step) % 7
            adj[ring_v].add(next_in_ring); adj[next_in_ring].add(ring_v)
    return [sorted(s) for s in adj]


# ---------------------------------------------------------------------------
# Spectral computations
# ---------------------------------------------------------------------------

def adjacency_matrix(adj):
    n = len(adj)
    A = np.zeros((n, n))
    for i, nbrs in enumerate(adj):
        for j in nbrs:
            A[i, j] = 1.0
    return A


def laplacian_matrix(adj):
    n = len(adj)
    A = adjacency_matrix(adj)
    D = np.diag([len(nbrs) for nbrs in adj])
    return D - A


def compute_spectrum(A):
    """Return eigenvalues sorted descending (real, rounded)."""
    evals = np.linalg.eigvalsh(A)
    evals = np.sort(evals)[::-1]
    return evals


def group_multiplicities(evals, tol=1e-8):
    """Group eigenvalues into (value, multiplicity) pairs."""
    groups = []
    i = 0
    while i < len(evals):
        val = evals[i]
        count = 1
        while i + count < len(evals) and abs(evals[i + count] - val) < tol:
            count += 1
        groups.append((round(val, 6), count))
        i += count
    return groups


def bfs_distances(adj, src):
    n = len(adj)
    dist = [-1] * n
    dist[src] = 0
    q = deque([src])
    while q:
        v = q.popleft()
        for w in adj[v]:
            if dist[w] == -1:
                dist[w] = dist[v] + 1
                q.append(w)
    return dist


def compute_diameter(adj):
    n = len(adj)
    return max(max(bfs_distances(adj, v)) for v in range(n))


def compute_girth(adj):
    """BFS from each vertex to find shortest cycle."""
    n = len(adj)
    best = float('inf')
    for src in range(n):
        dist = [-1] * n
        parent = [-1] * n
        dist[src] = 0
        q = deque([src])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if dist[w] == -1:
                    dist[w] = dist[v] + 1
                    parent[w] = v
                    q.append(w)
                elif parent[v] != w:
                    # Found a cycle of length dist[v] + dist[w] + 1
                    cycle_len = dist[v] + dist[w] + 1
                    if cycle_len < best:
                        best = cycle_len
    return best


def analyze_graph(name, adj):
    n = len(adj)
    A = adjacency_matrix(adj)
    L = laplacian_matrix(adj)

    adj_evals = compute_spectrum(A)
    lap_evals = np.sort(np.linalg.eigvalsh(L))

    lambda1 = adj_evals[0]
    lambda2 = adj_evals[1] if len(adj_evals) > 1 else None
    lambda_min = adj_evals[-1]

    spectral_gap = float(lambda1 - lambda2) if lambda2 is not None else None

    # Fiedler = smallest nonzero Laplacian eigenvalue
    lap_nonzero = [v for v in lap_evals if v > 1e-8]
    fiedler = float(min(lap_nonzero)) if lap_nonzero else 0.0

    groups = group_multiplicities(adj_evals)
    n_distinct = len(groups)

    ramanujan_bound = 2 * np.sqrt(3 - 1)  # 2*sqrt(d-1) for d=3
    is_ramanujan = bool(float(lambda2) <= ramanujan_bound + 1e-8)

    is_bipartite = bool(abs(float(lambda_min) + float(lambda1)) < 1e-8)

    diameter = compute_diameter(adj)
    girth = compute_girth(adj)

    # Distance-transitive iff #distinct eigenvalues = diameter + 1
    is_distance_transitive = (n_distinct == diameter + 1)

    return {
        "name": name,
        "n": n,
        "adj_eigenvalues": [float(x) for x in adj_evals],
        "eigenvalue_groups": [(float(v), m) for v, m in groups],
        "lambda1": float(lambda1),
        "lambda2": float(lambda2),
        "lambda_min": float(lambda_min),
        "spectral_gap": float(spectral_gap),
        "fiedler": fiedler,
        "n_distinct": n_distinct,
        "ramanujan_bound": float(ramanujan_bound),
        "is_ramanujan": is_ramanujan,
        "is_bipartite": is_bipartite,
        "diameter": diameter,
        "girth": girth,
        "is_distance_transitive": is_distance_transitive,
    }


# ---------------------------------------------------------------------------
# All Foster census graphs
# ---------------------------------------------------------------------------

GRAPHS = [
    ("K_4",           k4()),
    ("K_{3,3}",       k33()),
    ("Cube",          gp_graph(4, 1)),
    ("Petersen",      gp_graph(5, 2)),
    ("Heawood",       lcf_graph(14, [5, -5])),
    ("Moebius-Kantor", gp_graph(8, 3)),
    ("Pappus",        lcf_graph(18, [5, 7, -7, 7, -7, -5])),
    ("Dodecahedron",  dodecahedron()),
    ("Desargues",     gp_graph(10, 3)),
    ("Nauru",         gp_graph(12, 5)),
    ("F26A",          lcf_graph(26, [-7, 7])),
    ("Coxeter",       coxeter()),
    ("Tutte-Coxeter", lcf_graph(30, [-13, -9, 7, -7, 9, 13])),
]


def main():
    results = []
    for name, adj in GRAPHS:
        r = analyze_graph(name, adj)
        results.append(r)

    # Print table
    print(f"{'Graph':<20} {'n':>4} {'λ₂':>8} {'λ_min':>8} {'Gap':>8} "
          f"{'Fiedler':>8} {'#dist':>6} {'Gram':>4} {'Diam':>4} {'Girth':>5} "
          f"{'Bip':>4} {'Ram':>4} {'DistTrans':>10}")
    print("-" * 110)
    for r in results:
        print(f"{r['name']:<20} {r['n']:>4} {r['lambda2']:>8.4f} {r['lambda_min']:>8.4f} "
              f"{r['spectral_gap']:>8.4f} {r['fiedler']:>8.4f} {r['n_distinct']:>6} "
              f"{3:>4} {r['diameter']:>4} {r['girth']:>5} "
              f"{'Y' if r['is_bipartite'] else 'N':>4} "
              f"{'Y' if r['is_ramanujan'] else 'N':>4} "
              f"{'Y' if r['is_distance_transitive'] else 'N':>10}")

    print("\n\nEigenvalue spectra:")
    for r in results:
        print(f"\n{r['name']} (n={r['n']}):")
        for val, mult in r['eigenvalue_groups']:
            print(f"  {val:+.6f}  (mult {mult})")

    # Save JSON for LaTeX
    with open("/Users/robin/git/graphviz/spectral/foster_spectrum_data.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nData saved to foster_spectrum_data.json")

    return results


if __name__ == "__main__":
    main()
