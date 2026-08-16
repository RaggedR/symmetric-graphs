"""
g2_spectrum.py
Spectral analysis of the three graphs arising from G_2(2),
the automorphism group of the split Cayley hexagon GH(2,2).

Graphs:
  1. Tutte 12-cage: incidence graph of GH(2,2), 126 vertices, cubic, girth 12
  2. GH(2,2) collinearity graph: point graph of GH(2,2), 63 vertices, 6-regular
  3. Meinhold-72: quotient of Clayworth via PSL(3,2), 72 vertices, 7-regular

Author: Robin Langer
"""

import numpy as np
from collections import deque
from fractions import Fraction

# ============================================================
# Adjacency lists (1-indexed in source; converted to 0-indexed here)
# ============================================================

# --- Tutte 12-cage: 126 vertices, cubic ---
_tuttecage12_1indexed = [
    [ 64, 65, 66 ], [ 67, 68, 69 ], [ 70, 71, 72 ], [ 73, 74, 75 ],
    [ 76, 77, 78 ], [ 79, 80, 81 ], [ 82, 83, 84 ], [ 67, 73, 79 ],
    [ 67, 85, 86 ], [ 73, 87, 88 ], [ 79, 89, 90 ], [ 91, 92, 93 ],
    [ 94, 95, 96 ], [ 97, 98, 99 ], [ 100, 101, 102 ], [ 64, 74, 76 ],
    [ 64, 103, 104 ], [ 74, 91, 97 ], [ 76, 105, 106 ], [ 87, 107, 108 ],
    [ 100, 109, 110 ], [ 88, 111, 112 ], [ 94, 113, 114 ], [ 70, 75, 82 ],
    [ 70, 115, 116 ], [ 75, 94, 100 ], [ 82, 117, 118 ], [ 97, 119, 120 ],
    [ 88, 121, 122 ], [ 91, 123, 124 ], [ 87, 125, 126 ], [ 65, 68, 71 ],
    [ 65, 107, 109 ], [ 68, 92, 95 ], [ 71, 119, 121 ], [ 85, 103, 115 ],
    [ 101, 104, 125 ], [ 86, 113, 123 ], [ 98, 111, 116 ], [ 69, 77, 83 ],
    [ 69, 98, 101 ], [ 77, 108, 114 ], [ 83, 122, 124 ], [ 95, 105, 126 ],
    [ 92, 112, 117 ], [ 86, 106, 118 ], [ 85, 110, 120 ], [ 66, 80, 84 ],
    [ 66, 111, 113 ], [ 80, 93, 102 ], [ 84, 120, 126 ], [ 89, 109, 124 ],
    [ 99, 107, 118 ], [ 90, 104, 117 ], [ 96, 103, 122 ], [ 72, 78, 81 ],
    [ 72, 123, 125 ], [ 78, 110, 112 ], [ 81, 96, 99 ], [ 102, 106, 121 ],
    [ 90, 114, 119 ], [ 89, 105, 116 ], [ 93, 108, 115 ],
    [ 1, 16, 17 ], [ 1, 32, 33 ], [ 1, 48, 49 ], [ 2, 8, 9 ],
    [ 2, 32, 34 ], [ 2, 40, 41 ], [ 3, 24, 25 ], [ 3, 32, 35 ],
    [ 3, 56, 57 ], [ 4, 8, 10 ], [ 4, 16, 18 ], [ 4, 24, 26 ],
    [ 5, 16, 19 ], [ 5, 40, 42 ], [ 5, 56, 58 ], [ 6, 8, 11 ],
    [ 6, 48, 50 ], [ 6, 56, 59 ], [ 7, 24, 27 ], [ 7, 40, 43 ],
    [ 7, 48, 51 ], [ 9, 36, 47 ], [ 9, 38, 46 ], [ 10, 20, 31 ],
    [ 10, 22, 29 ], [ 11, 52, 62 ], [ 11, 54, 61 ], [ 12, 18, 30 ],
    [ 12, 34, 45 ], [ 12, 50, 63 ], [ 13, 23, 26 ], [ 13, 34, 44 ],
    [ 13, 55, 59 ], [ 14, 18, 28 ], [ 14, 39, 41 ], [ 14, 53, 59 ],
    [ 15, 21, 26 ], [ 15, 37, 41 ], [ 15, 50, 60 ], [ 17, 36, 55 ],
    [ 17, 37, 54 ], [ 19, 44, 62 ], [ 19, 46, 60 ], [ 20, 33, 53 ],
    [ 20, 42, 63 ], [ 21, 33, 52 ], [ 21, 47, 58 ], [ 22, 39, 49 ],
    [ 22, 45, 58 ], [ 23, 38, 49 ], [ 23, 42, 61 ], [ 25, 36, 63 ],
    [ 25, 39, 62 ], [ 27, 45, 54 ], [ 27, 46, 53 ], [ 28, 35, 61 ],
    [ 28, 47, 51 ], [ 29, 35, 60 ], [ 29, 43, 55 ], [ 30, 38, 57 ],
    [ 30, 43, 52 ], [ 31, 37, 57 ], [ 31, 44, 51 ]
]

# --- GH(2,2) collinearity graph: 63 vertices, 6-regular ---
_gh22_1indexed = [
    [ 19, 20, 31, 41, 44, 60 ], [ 15, 16, 35, 38, 49, 51 ],
    [ 23, 24, 27, 36, 51, 60 ], [ 5, 21, 40, 44, 45, 57 ],
    [ 4, 25, 38, 43, 47, 57 ], [ 7, 17, 48, 51, 52, 63 ],
    [ 6, 28, 29, 39, 45, 63 ], [ 10, 13, 55, 60, 61, 62 ],
    [ 11, 12, 33, 44, 56, 58 ], [ 8, 26, 32, 34, 43, 62 ],
    [ 9, 25, 32, 46, 52, 58 ], [ 9, 14, 30, 38, 50, 56 ],
    [ 8, 22, 29, 30, 37, 55 ], [ 12, 21, 28, 50, 53, 61 ],
    [ 2, 28, 33, 37, 49, 54 ], [ 2, 18, 34, 35, 45, 59 ],
    [ 6, 18, 26, 33, 42, 48 ], [ 16, 17, 25, 42, 59, 61 ],
    [ 1, 30, 32, 41, 42, 54 ], [ 1, 22, 31, 39, 43, 59 ],
    [ 4, 14, 22, 40, 52, 53 ], [ 13, 20, 21, 37, 52, 59 ],
    [ 3, 27, 34, 37, 47, 50 ], [ 3, 26, 29, 36, 46, 53 ],
    [ 5, 11, 18, 46, 47, 61 ], [ 10, 17, 24, 33, 43, 53 ],
    [ 3, 23, 39, 40, 42, 58 ], [ 7, 14, 15, 39, 54, 61 ],
    [ 7, 13, 24, 30, 45, 46 ], [ 12, 13, 19, 29, 38, 42 ],
    [ 1, 20, 46, 48, 49, 50 ], [ 10, 11, 19, 34, 52, 54 ],
    [ 9, 15, 17, 26, 37, 44 ], [ 10, 16, 23, 32, 45, 50 ],
    [ 2, 16, 41, 53, 55, 58 ], [ 3, 24, 54, 56, 57, 59 ],
    [ 13, 15, 22, 23, 33, 47 ], [ 2, 5, 12, 30, 43, 51 ],
    [ 7, 20, 27, 28, 43, 58 ], [ 4, 21, 27, 42, 49, 62 ],
    [ 1, 19, 35, 47, 53, 63 ], [ 17, 18, 19, 27, 30, 40 ],
    [ 5, 10, 20, 26, 38, 39 ], [ 1, 4, 9, 33, 45, 60 ],
    [ 4, 7, 16, 29, 34, 44 ], [ 11, 24, 25, 29, 31, 49 ],
    [ 5, 23, 25, 37, 41, 63 ], [ 6, 17, 31, 50, 55, 57 ],
    [ 2, 15, 31, 40, 46, 62 ], [ 12, 14, 23, 31, 34, 48 ],
    [ 2, 3, 6, 38, 52, 60 ], [ 6, 11, 21, 22, 32, 51 ],
    [ 14, 21, 24, 26, 35, 41 ], [ 15, 19, 28, 32, 36, 57 ],
    [ 8, 13, 35, 48, 57, 58 ], [ 9, 12, 36, 59, 62, 63 ],
    [ 4, 5, 36, 48, 54, 55 ], [ 9, 11, 27, 35, 39, 55 ],
    [ 16, 18, 20, 22, 36, 56 ], [ 1, 3, 8, 44, 51, 61 ],
    [ 8, 14, 18, 25, 28, 60 ], [ 8, 10, 40, 49, 56, 63 ],
    [ 6, 7, 41, 47, 56, 62 ]
]

# --- Meinhold-72: 72 vertices, 7-regular ---
_meinhold72_1indexed = [
    [ 4, 24, 32, 38, 43, 52, 59 ],
    [ 10, 15, 17, 20, 33, 44, 69 ],
    [ 10, 11, 20, 29, 34, 37, 47 ],
    [ 1, 14, 18, 28, 48, 65, 71 ],
    [ 11, 22, 25, 35, 37, 56, 66 ],
    [ 7, 15, 22, 33, 38, 59, 66 ],
    [ 6, 8, 23, 39, 41, 54, 57 ],
    [ 7, 9, 22, 24, 25, 32, 53 ],
    [ 8, 12, 23, 26, 58, 61, 64 ],
    [ 2, 3, 31, 39, 40, 51, 57 ],
    [ 3, 5, 14, 36, 40, 65, 72 ],
    [ 9, 25, 30, 44, 56, 60, 69 ],
    [ 26, 45, 49, 50, 61, 67, 70 ],
    [ 4, 11, 16, 20, 43, 56, 69 ],
    [ 2, 6, 18, 45, 57, 70, 71 ],
    [ 14, 18, 26, 36, 45, 46, 64 ],
    [ 2, 18, 19, 23, 39, 48, 64 ],
    [ 4, 15, 16, 17, 21, 38, 63 ],
    [ 17, 21, 30, 34, 44, 47, 68 ],
    [ 2, 3, 14, 45, 48, 49, 62 ],
    [ 18, 19, 27, 28, 46, 50, 70 ],
    [ 5, 6, 8, 36, 45, 61, 62 ],
    [ 7, 9, 17, 30, 33, 55, 63 ],
    [ 1, 8, 31, 39, 48, 58, 62 ],
    [ 5, 8, 12, 31, 57, 65, 71 ],
    [ 9, 13, 16, 32, 43, 60, 63 ],
    [ 21, 35, 38, 42, 52, 66, 68 ],
    [ 4, 21, 30, 55, 56, 59, 66 ],
    [ 3, 41, 49, 51, 54, 67, 72 ],
    [ 12, 19, 23, 28, 54, 65, 72 ],
    [ 10, 24, 25, 34, 35, 44, 52 ],
    [ 1, 8, 26, 36, 54, 65, 67 ],
    [ 2, 6, 23, 40, 49, 61, 72 ],
    [ 3, 19, 31, 48, 50, 65, 67 ],
    [ 5, 27, 31, 40, 50, 58, 61 ],
    [ 11, 16, 22, 32, 38, 47, 68 ],
    [ 3, 5, 41, 45, 46, 50, 57 ],
    [ 1, 6, 18, 27, 36, 39, 40 ],
    [ 7, 10, 17, 24, 38, 42, 47 ],
    [ 10, 11, 33, 35, 38, 43, 63 ],
    [ 7, 29, 37, 42, 53, 55, 66 ],
    [ 27, 39, 41, 46, 51, 58, 64 ],
    [ 1, 14, 26, 40, 49, 51, 58 ],
    [ 2, 12, 19, 31, 61, 62, 70 ],
    [ 13, 15, 16, 20, 22, 37, 53 ],
    [ 16, 21, 37, 42, 47, 56, 60 ],
    [ 3, 19, 36, 39, 46, 54, 62 ],
    [ 4, 17, 20, 24, 34, 53, 55 ],
    [ 13, 20, 29, 33, 43, 55, 59 ],
    [ 13, 21, 34, 35, 37, 55, 63 ],
    [ 10, 29, 42, 43, 52, 60, 69 ],
    [ 1, 27, 31, 51, 67, 70, 71 ],
    [ 8, 41, 45, 48, 64, 67, 71 ],
    [ 7, 29, 30, 32, 47, 59, 60 ],
    [ 23, 28, 41, 48, 49, 50, 58 ],
    [ 5, 12, 14, 28, 46, 58, 62 ],
    [ 7, 10, 15, 25, 37, 60, 63 ],
    [ 9, 24, 35, 42, 43, 55, 56 ],
    [ 1, 6, 28, 49, 54, 62, 70 ],
    [ 12, 26, 46, 51, 54, 57, 70 ],
    [ 9, 13, 22, 33, 35, 44, 68 ],
    [ 20, 22, 24, 44, 47, 56, 59 ],
    [ 18, 23, 26, 40, 50, 57, 65 ],
    [ 9, 16, 17, 42, 53, 68, 69 ],
    [ 4, 11, 25, 30, 32, 34, 63 ],
    [ 5, 6, 27, 28, 41, 71, 72 ],
    [ 13, 29, 32, 34, 52, 53, 68 ],
    [ 19, 27, 36, 61, 64, 67, 72 ],
    [ 2, 12, 14, 51, 64, 71, 72 ],
    [ 13, 15, 21, 44, 52, 59, 60 ],
    [ 4, 15, 25, 52, 53, 66, 69 ],
    [ 11, 29, 30, 33, 66, 68, 69 ]
]


def to_zero_indexed(adj_1indexed):
    """Convert 1-indexed adjacency list to 0-indexed."""
    return [[v - 1 for v in neighbors] for neighbors in adj_1indexed]


def build_adjacency_matrix(adj, n):
    """Build adjacency matrix from adjacency list."""
    A = np.zeros((n, n), dtype=float)
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            A[u, v] = 1.0
    return A


def compute_spectrum(adj, n):
    """Compute eigenvalues of the adjacency matrix, sorted descending."""
    A = build_adjacency_matrix(adj, n)
    eigenvalues = np.linalg.eigvalsh(A)
    eigenvalues = np.sort(eigenvalues)[::-1]  # descending
    return eigenvalues


def round_nice(x, tol=1e-9):
    """Round eigenvalue to nearest integer or half-integer if close."""
    r = round(x)
    if abs(x - r) < tol:
        return float(r)
    h = round(2 * x) / 2
    if abs(x - h) < tol:
        return h
    return x


def group_eigenvalues(eigenvalues, tol=1e-6):
    """Group eigenvalues by value (with tolerance), return (value, multiplicity) pairs."""
    groups = []
    current = eigenvalues[0]
    count = 1
    for ev in eigenvalues[1:]:
        if abs(ev - current) < tol:
            count += 1
        else:
            groups.append((round_nice(current), count))
            current = ev
            count = 1
    groups.append((round_nice(current), count))
    return groups


def bfs_girth(adj, n):
    """Compute girth using BFS from each vertex."""
    girth = float('inf')
    for start in range(n):
        dist = [-1] * n
        dist[start] = 0
        parent = [-1] * n
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    queue.append(v)
                elif parent[u] != v:
                    # Found a cycle: length dist[u] + dist[v] + 1
                    cycle_len = dist[u] + dist[v] + 1
                    if cycle_len < girth:
                        girth = cycle_len
        if girth == 3:
            break  # can't do better
    return girth


def check_bipartite(adj, n):
    """Check if graph is bipartite using 2-coloring."""
    color = [-1] * n
    for start in range(n):
        if color[start] != -1:
            continue
        color[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if color[v] == -1:
                    color[v] = 1 - color[u]
                    queue.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def analyze_graph(name, n, regularity, adj_1indexed, ramanujan_bound_latex):
    """Full spectral analysis of a regular graph."""
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"  n={n}, {regularity}-regular")
    print(f"{'='*60}")

    adj = to_zero_indexed(adj_1indexed)

    # Verify regularity
    degrees = [len(neighbors) for neighbors in adj]
    assert all(d == regularity for d in degrees), \
        f"Degree check failed: {set(degrees)}"
    print(f"  Regularity check: OK (all degrees = {regularity})")

    # Spectrum
    eigenvalues = compute_spectrum(adj, n)
    spectrum = group_eigenvalues(eigenvalues)

    print(f"  Spectrum (descending, {len(spectrum)} distinct values):")
    for val, mult in spectrum:
        print(f"    λ = {val:+.6f}  (multiplicity {mult})")

    # Ramanujan check
    lambda1 = eigenvalues[0]
    lambda2 = eigenvalues[1]  # second-largest
    lambda_min = eigenvalues[-1]
    ramanujan_threshold = 2 * np.sqrt(regularity - 1)
    is_ramanujan = max(abs(lambda2), abs(lambda_min)) <= ramanujan_threshold + 1e-9

    print(f"\n  λ₁ = {lambda1:.6f}  (should be {regularity})")
    print(f"  λ₂ = {lambda2:.6f}")
    print(f"  λ_min = {lambda_min:.6f}")
    print(f"  Ramanujan threshold = 2√{regularity-1} ≈ {ramanujan_threshold:.6f}")
    print(f"  Ramanujan: {'YES' if is_ramanujan else 'NO'}")

    # Fiedler value and spectral gap
    # For connected graph, eigenvalues sorted descending
    # Algebraic connectivity = second smallest = index n-2
    fiedler = eigenvalues[-2]
    spectral_gap = lambda1 - lambda2  # gap from top

    print(f"\n  Fiedler eigenvalue (algebraic connectivity) = {fiedler:.6f}")
    print(f"  Spectral gap (λ₁ - λ₂) = {spectral_gap:.6f}")
    print(f"  Number of distinct eigenvalues: {len(spectrum)}")

    # Bipartiteness
    is_bipartite = check_bipartite(adj, n)
    spectrum_symmetric = abs(lambda1 + lambda_min) < 1e-6
    print(f"  Bipartite: {is_bipartite}")
    print(f"  Spectrum symmetric about 0: {spectrum_symmetric}")

    # Girth
    g = bfs_girth(adj, n)
    print(f"  Girth: {g}")

    return {
        "name": name,
        "n": n,
        "regularity": regularity,
        "spectrum": spectrum,
        "eigenvalues_raw": eigenvalues,
        "lambda1": lambda1,
        "lambda2": lambda2,
        "lambda_min": lambda_min,
        "ramanujan_threshold": ramanujan_threshold,
        "is_ramanujan": is_ramanujan,
        "fiedler": fiedler,
        "spectral_gap": spectral_gap,
        "num_distinct_eigenvalues": len(spectrum),
        "is_bipartite": is_bipartite,
        "spectrum_symmetric": spectrum_symmetric,
        "girth": g,
    }


def print_latex_table(result):
    """Print a LaTeX-friendly spectrum table."""
    print(f"\n% Spectrum table for {result['name']}")
    print(r"\begin{tabular}{rr}")
    print(r"  \hline")
    print(r"  Eigenvalue & Multiplicity \\")
    print(r"  \hline")
    for val, mult in result["spectrum"]:
        if val == int(val):
            print(f"  ${int(val)}$ & ${mult}$ \\\\")
        else:
            # Try to represent as sqrt fraction
            print(f"  ${val:.4f}$ & ${mult}$ \\\\")
    print(r"  \hline")
    print(r"\end{tabular}")


def main():
    print("Spectral Analysis: G₂(2) Graphs")
    print("=" * 60)

    results = []

    # 1. Tutte 12-cage
    r1 = analyze_graph(
        name="Tutte 12-cage (incidence graph of GH(2,2))",
        n=126,
        regularity=3,
        adj_1indexed=_tuttecage12_1indexed,
        ramanujan_bound_latex=r"2\sqrt{2}",
    )
    results.append(r1)

    # 2. GH(2,2) collinearity graph
    r2 = analyze_graph(
        name="GH(2,2) collinearity graph (point graph)",
        n=63,
        regularity=6,
        adj_1indexed=_gh22_1indexed,
        ramanujan_bound_latex=r"2\sqrt{5}",
    )
    results.append(r2)

    # 3. Meinhold-72
    r3 = analyze_graph(
        name="Meinhold-72 (Clayworth/PSL(3,2) quotient)",
        n=72,
        regularity=7,
        adj_1indexed=_meinhold72_1indexed,
        ramanujan_bound_latex=r"2\sqrt{6}",
    )
    results.append(r3)

    # Summary comparison table
    print("\n\n" + "=" * 60)
    print("SUMMARY COMPARISON")
    print("=" * 60)
    print(f"{'Graph':<35} {'n':>4} {'k':>3} {'λ₂':>8} {'λ_min':>8} {'Ramanujan':>10} {'Girth':>6} {'#ev':>4}")
    print("-" * 80)
    for r in results:
        print(
            f"{r['name'][:35]:<35} {r['n']:>4} {r['regularity']:>3} "
            f"{r['lambda2']:>8.4f} {r['lambda_min']:>8.4f} "
            f"{'YES' if r['is_ramanujan'] else 'NO':>10} "
            f"{r['girth']:>6} {r['num_distinct_eigenvalues']:>4}"
        )

    # LaTeX tables
    print("\n\n% ==================== LaTeX tables ====================")
    for r in results:
        print_latex_table(r)


if __name__ == "__main__":
    main()
