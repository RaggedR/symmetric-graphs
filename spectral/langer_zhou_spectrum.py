"""
langer_zhou_spectrum.py

Spectral analysis of four graphs from exceptional/classical groups over small fields:
  1. Langer graph: collinearity graph of GH(2,2), 63 vertices, 6-regular
  2. Dual Langer graph: line-collinearity graph of GH(2,2), 63 vertices, 6-regular
  3. Zhou graph: Sab(PSL(2,13), D_12), 91 vertices, 6-regular
  4. Clayworth graph: Sab(G_2(2), C_3), 4032 vertices, cubic

Key question: Are the Langer and Dual Langer cospectral?

Author: Robin Langer
"""

import numpy as np
import json
from collections import deque

# ============================================================
# Build the Langer graph (Zorn product formula on Q(6,2))
# ============================================================

def to_bit_vec(n):
    return tuple((n >> i) & 1 for i in range(7))

def q62form(x):
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

q62_indices = [n for n in range(1, 128) if q62form(to_bit_vec(n)) == 0]
assert len(q62_indices) == 63

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

print("Building Langer graph...")
langer_nbrs = [[j for j in range(63) if langer_adj(i, j)] for i in range(63)]
for i in range(63):
    assert len(langer_nbrs[i]) == 6, f"Vertex {i} has degree {len(langer_nbrs[i])}"
print(f"  Langer: 63 vertices, {sum(len(n) for n in langer_nbrs)//2} edges")

# ============================================================
# Build the Dual Langer graph (line-collinearity of GH(2,2))
# ============================================================

print("Building Dual Langer graph...")
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
assert len(triangles) == 63, f"Expected 63 triangles, got {len(triangles)}"

dual_nbrs = [[] for _ in range(63)]
for i in range(63):
    for j in range(i+1, 63):
        if set(triangles[i]) & set(triangles[j]):
            dual_nbrs[i].append(j)
            dual_nbrs[j].append(i)

for i in range(63):
    assert len(dual_nbrs[i]) == 6, f"Dual vertex {i} has degree {len(dual_nbrs[i])}"
print(f"  Dual Langer: 63 vertices, {sum(len(n) for n in dual_nbrs)//2} edges")

# ============================================================
# Zhou graph adjacency list (1-indexed from GAP)
# ============================================================

ADJ_ZHOU_1 = [
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

zhou_nbrs = [[u-1 for u in nbrs] for nbrs in ADJ_ZHOU_1]
for i in range(91):
    assert len(zhou_nbrs[i]) == 6
print(f"  Zhou: 91 vertices, {sum(len(n) for n in zhou_nbrs)//2} edges")


# ============================================================
# Utility functions
# ============================================================

def adjacency_matrix(nbrs, n):
    A = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in nbrs[i]:
            A[i, j] = 1.0
    return A

def compute_spectrum(A):
    """Compute sorted eigenvalues (descending) of a symmetric matrix."""
    evals = np.linalg.eigvalsh(A)
    evals = np.sort(evals)[::-1]
    return evals

def round_to_nice(x, tol=1e-8):
    """Round eigenvalue to a nice form if close to a recognisable algebraic number."""
    # Check integers first
    r = round(x)
    if abs(x - r) < tol:
        return float(r)
    # Check half-integers
    r2 = round(2*x) / 2.0
    if abs(x - r2) < tol:
        return r2
    return x

def group_eigenvalues(evals, tol=1e-6):
    """Group eigenvalues by value, returning list of (value, multiplicity)."""
    groups = []
    for ev in evals:
        ev_nice = round_to_nice(ev)
        if groups and abs(ev_nice - groups[-1][0]) < tol:
            groups[-1][1] += 1
        else:
            groups.append([ev_nice, 1])
    return [(v, m) for v, m in groups]

def bfs_distances(nbrs, n, source):
    """BFS from source; returns distance array."""
    dist = [-1] * n
    dist[source] = 0
    q = deque([source])
    while q:
        v = q.popleft()
        for u in nbrs[v]:
            if dist[u] == -1:
                dist[u] = dist[v] + 1
                q.append(u)
    return dist

def compute_diameter(nbrs, n):
    diam = 0
    for s in range(n):
        d = bfs_distances(nbrs, n, s)
        diam = max(diam, max(d))
    return diam

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
                elif parent[v] != u:
                    cycle_len = dist[v] + dist[u] + 1
                    min_girth = min(min_girth, cycle_len)
        if min_girth <= 4:
            return min_girth
    return min_girth

def compute_intersection_array(nbrs, n, k):
    """
    For a distance-regular graph: compute intersection numbers b_i, c_i.
    Returns (b_list, c_list) where b_list[i] = b_i, c_list[i] = c_i.
    Uses vertex 0 as reference; checks consistency across a sample.
    """
    # BFS layers from vertex 0
    dist = bfs_distances(nbrs, n, 0)
    diameter = max(dist)

    b = []  # b[i] = |Gamma_{i+1}(u) cap Gamma(v)| for v in Gamma_i(0)
    c = []  # c[i] = |Gamma_{i-1}(u) cap Gamma(v)| for v in Gamma_i(0)

    for i in range(diameter + 1):
        layer_i = [v for v in range(n) if dist[v] == i]
        # Pick first vertex u in layer i, pick any v in Gamma(u) at dist i
        u = layer_i[0]
        # b_i: neighbours of u at dist i+1 from source
        b_count = sum(1 for w in nbrs[u] if dist[w] == i+1)
        # c_i: neighbours of u at dist i-1 from source
        c_count = sum(1 for w in nbrs[u] if dist[w] == i-1)
        b.append(b_count)
        c.append(c_count)

    return b[:-1], c[1:]  # b_0..b_{d-1}, c_1..c_d

def is_distance_regular(nbrs, n, k, tol=1e-6):
    """Check distance-regularity by verifying intersection numbers are vertex-independent."""
    dist0 = bfs_distances(nbrs, n, 0)
    diameter = max(dist0)

    # Compute intersection numbers from vertex 0
    b_ref = []
    c_ref = []
    for i in range(diameter + 1):
        layer_i = [v for v in range(n) if dist0[v] == i]
        u = layer_i[0]
        b_ref.append(sum(1 for w in nbrs[u] if dist0[w] == i+1))
        c_ref.append(sum(1 for w in nbrs[u] if dist0[w] == i-1))

    # Verify from a sample of other vertices
    import random
    sample = random.sample(range(n), min(10, n))
    for s in sample:
        dist_s = bfs_distances(nbrs, n, s)
        for i in range(diameter + 1):
            layer_i = [v for v in range(n) if dist_s[v] == i]
            if not layer_i:
                continue
            u = layer_i[0]
            b_i = sum(1 for w in nbrs[u] if dist_s[w] == i+1)
            c_i = sum(1 for w in nbrs[u] if dist_s[w] == i-1)
            if b_i != b_ref[i] or c_i != c_ref[i]:
                return False, None, None
    return True, b_ref[:-1], c_ref[1:]

def fiedler_value(A, n):
    """Compute algebraic connectivity (Fiedler value) = 2nd smallest Laplacian eigenvalue."""
    k = A.sum(axis=1)[0]  # degree
    L = np.diag(A.sum(axis=1)) - A
    evals_L = np.linalg.eigvalsh(L)
    evals_L = np.sort(evals_L)
    return float(evals_L[1])

def check_ramanujan(evals_sorted, k):
    """
    For a k-regular graph, check if |lambda_2| <= 2*sqrt(k-1).
    Returns (is_ramanujan, lambda_2, bound).
    """
    # evals_sorted is descending, so evals_sorted[0] = k
    bound = 2 * np.sqrt(k - 1)
    # second largest in absolute value (excluding k and possibly -k)
    abs_vals = sorted([abs(e) for e in evals_sorted], reverse=True)
    # abs_vals[0] = k (the largest), lambda_2 is the next one
    lambda_2 = abs_vals[1]
    return lambda_2 <= bound + 1e-6, float(lambda_2), float(bound)

def d3_connectivity(nbrs, n):
    """
    Compute the number of edges in the distance-3 subgraph.
    That is, count pairs (u,v) with d(u,v)=3 that are adjacent.
    Actually: count edges in the induced subgraph on distance-3 layer from vertex 0,
    then check consistency. Better: compute the total number of edges in G^{(3)}
    (the graph where u~v iff d(u,v)=3).
    """
    count = 0
    for s in range(n):
        dist = bfs_distances(nbrs, n, s)
        d3 = [v for v in range(n) if dist[v] == 3]
        # edges among d3 vertices
        for v in d3:
            for w in nbrs[v]:
                if dist[w] == 3 and w > v:
                    count += 1
    return count


# ============================================================
# Run spectral analysis
# ============================================================

import random
random.seed(42)

results = {}

print("\n" + "="*60)
print("LANGER GRAPH (63 vertices, 6-regular)")
print("="*60)

A_L = adjacency_matrix(langer_nbrs, 63)
evals_L = compute_spectrum(A_L)
groups_L = group_eigenvalues(evals_L)
diam_L = compute_diameter(langer_nbrs, 63)
girth_L = compute_girth(langer_nbrs, 63)
is_dr_L, b_L, c_L = is_distance_regular(langer_nbrs, 63, 6)
fiedler_L = fiedler_value(A_L, 63)
ram_L, lam2_L, bound_L = check_ramanujan(evals_L, 6)
d3_L = d3_connectivity(langer_nbrs, 63)

print(f"Diameter: {diam_L}")
print(f"Girth: {girth_L}")
print(f"Distance-regular: {is_dr_L}")
if is_dr_L:
    print(f"  b = {b_L}")
    print(f"  c = {c_L}")
    print(f"  Intersection array: {{{','.join(map(str,b_L))};{','.join(map(str,c_L))}}}")
print(f"Fiedler value (lambda_2 of Laplacian): {fiedler_L:.6f}")
print(f"Ramanujan: {ram_L} (|lambda_2| = {lam2_L:.6f}, bound = {bound_L:.6f})")
print(f"D3-connectivity count: {d3_L}")
print(f"Eigenvalues (value, multiplicity):")
for v, m in groups_L:
    print(f"  lambda = {v:>10.6f}  (multiplicity {m})")

results['langer'] = {
    'n': 63, 'k': 6,
    'diameter': diam_L,
    'girth': girth_L,
    'distance_regular': is_dr_L,
    'b_array': b_L,
    'c_array': c_L,
    'fiedler': fiedler_L,
    'ramanujan': ram_L,
    'lambda2_abs': lam2_L,
    'ramanujan_bound': bound_L,
    'd3_connectivity': d3_L,
    'eigenvalues': [(float(v), int(m)) for v,m in groups_L],
    'raw_evals': [float(e) for e in evals_L],
}

print("\n" + "="*60)
print("DUAL LANGER GRAPH (63 vertices, 6-regular)")
print("="*60)

A_DL = adjacency_matrix(dual_nbrs, 63)
evals_DL = compute_spectrum(A_DL)
groups_DL = group_eigenvalues(evals_DL)
diam_DL = compute_diameter(dual_nbrs, 63)
girth_DL = compute_girth(dual_nbrs, 63)
is_dr_DL, b_DL, c_DL = is_distance_regular(dual_nbrs, 63, 6)
fiedler_DL = fiedler_value(A_DL, 63)
ram_DL, lam2_DL, bound_DL = check_ramanujan(evals_DL, 6)
d3_DL = d3_connectivity(dual_nbrs, 63)

print(f"Diameter: {diam_DL}")
print(f"Girth: {girth_DL}")
print(f"Distance-regular: {is_dr_DL}")
if is_dr_DL:
    print(f"  b = {b_DL}")
    print(f"  c = {c_DL}")
    print(f"  Intersection array: {{{','.join(map(str,b_DL))};{','.join(map(str,c_DL))}}}")
print(f"Fiedler value: {fiedler_DL:.6f}")
print(f"Ramanujan: {ram_DL} (|lambda_2| = {lam2_DL:.6f}, bound = {bound_DL:.6f})")
print(f"D3-connectivity count: {d3_DL}")
print(f"Eigenvalues (value, multiplicity):")
for v, m in groups_DL:
    print(f"  lambda = {v:>10.6f}  (multiplicity {m})")

results['dual_langer'] = {
    'n': 63, 'k': 6,
    'diameter': diam_DL,
    'girth': girth_DL,
    'distance_regular': is_dr_DL,
    'b_array': b_DL,
    'c_array': c_DL,
    'fiedler': fiedler_DL,
    'ramanujan': ram_DL,
    'lambda2_abs': lam2_DL,
    'ramanujan_bound': bound_DL,
    'd3_connectivity': d3_DL,
    'eigenvalues': [(float(v), int(m)) for v,m in groups_DL],
    'raw_evals': [float(e) for e in evals_DL],
}

# Cospectrality check
evals_L_rounded = [round(e, 8) for e in sorted(evals_L)]
evals_DL_rounded = [round(e, 8) for e in sorted(evals_DL)]
cospectral = all(abs(a-b) < 1e-6 for a,b in zip(evals_L_rounded, evals_DL_rounded))

print(f"\nCOSPECTRAL (Langer vs Dual Langer): {cospectral}")
max_diff = max(abs(a-b) for a,b in zip(evals_L_rounded, evals_DL_rounded))
print(f"Max eigenvalue difference: {max_diff:.2e}")
results['cospectral_langer_dual'] = cospectral
results['max_eval_diff'] = float(max_diff)

print("\n" + "="*60)
print("ZHOU GRAPH (91 vertices, 6-regular)")
print("="*60)

A_Z = adjacency_matrix(zhou_nbrs, 91)
evals_Z = compute_spectrum(A_Z)
groups_Z = group_eigenvalues(evals_Z)
diam_Z = compute_diameter(zhou_nbrs, 91)
girth_Z = compute_girth(zhou_nbrs, 91)
is_dr_Z, b_Z, c_Z = is_distance_regular(zhou_nbrs, 91, 6)
fiedler_Z = fiedler_value(A_Z, 91)
ram_Z, lam2_Z, bound_Z = check_ramanujan(evals_Z, 6)
d3_Z = d3_connectivity(zhou_nbrs, 91)

print(f"Diameter: {diam_Z}")
print(f"Girth: {girth_Z}")
print(f"Distance-regular: {is_dr_Z}")
if is_dr_Z:
    print(f"  b = {b_Z}")
    print(f"  c = {c_Z}")
    print(f"  Intersection array: {{{','.join(map(str,b_Z))};{','.join(map(str,c_Z))}}}")
print(f"Fiedler value: {fiedler_Z:.6f}")
print(f"Ramanujan: {ram_Z} (|lambda_2| = {lam2_Z:.6f}, bound = {bound_Z:.6f})")
print(f"D3-connectivity count: {d3_Z}")
print(f"Eigenvalues (value, multiplicity):")
for v, m in groups_Z:
    print(f"  lambda = {v:>10.6f}  (multiplicity {m})")

# Check for sqrt(13) structure in Zhou
sqrt13 = np.sqrt(13)
for v, m in groups_Z:
    ratio = abs(v) / sqrt13 if abs(v) > 0.01 else 0
    if abs(ratio - round(ratio)) < 1e-5:
        print(f"  -> lambda = {v:.6f} ~ {round(ratio):.0f}*sqrt(13) = {round(ratio)*sqrt13:.6f}")

results['zhou'] = {
    'n': 91, 'k': 6,
    'diameter': diam_Z,
    'girth': girth_Z,
    'distance_regular': is_dr_Z,
    'b_array': b_Z,
    'c_array': c_Z,
    'fiedler': fiedler_Z,
    'ramanujan': ram_Z,
    'lambda2_abs': lam2_Z,
    'ramanujan_bound': bound_Z,
    'd3_connectivity': d3_Z,
    'eigenvalues': [(float(v), int(m)) for v,m in groups_Z],
    'raw_evals': [float(e) for e in evals_Z],
}

# ============================================================
# Identify algebraic forms for eigenvalues
# ============================================================

def identify_algebraic(x, tol=1e-7):
    """Try to identify x as a simple algebraic number."""
    # Integers
    r = round(x)
    if abs(x - r) < tol:
        return str(int(r))
    # sqrt-free quadratics: try x = a + b*sqrt(c) for small a,b,c
    for c in [2, 3, 5, 7, 13, 17, 21]:
        sq = np.sqrt(c)
        # x = a + b*sq => b = (x-a)/sq, check if b is rational
        for a in np.arange(-6, 7, 0.5):
            b = (x - a) / sq
            b_r = round(2*b) / 2.0
            if abs(b - b_r) < tol and abs(b_r) <= 4:
                if abs(a) < 1e-9:
                    if b_r == 1.0: return f"sqrt({c})"
                    if b_r == -1.0: return f"-sqrt({c})"
                    if b_r > 0: return f"{b_r:.4g}*sqrt({c})"
                    return f"{b_r:.4g}*sqrt({c})"
                else:
                    a_r = round(2*a) / 2.0
                    if abs(a_r) < 1e-9:
                        continue
                    sign = "+" if b_r > 0 else ""
                    return f"{a_r:.4g}{sign}{b_r:.4g}*sqrt({c})"
    return f"{x:.6f}"

print("\n" + "="*60)
print("ALGEBRAIC FORMS")
print("="*60)

print("\nLanger eigenvalues:")
for v, m in groups_L:
    alg = identify_algebraic(v)
    print(f"  {v:>10.6f} (x{m}) ~ {alg}")

print("\nDual Langer eigenvalues:")
for v, m in groups_DL:
    alg = identify_algebraic(v)
    print(f"  {v:>10.6f} (x{m}) ~ {alg}")

print("\nZhou eigenvalues:")
for v, m in groups_Z:
    alg = identify_algebraic(v)
    print(f"  {v:>10.6f} (x{m}) ~ {alg}")

# ============================================================
# Non-isomorphism check via d3-subgraph
# ============================================================

print("\n" + "="*60)
print("NON-ISOMORPHISM: d3-SUBGRAPH ANALYSIS")
print("="*60)

# Count degree distribution in d3-subgraph
def d3_degree_sequence(nbrs, n):
    """For each vertex, count how many of its distance-3 vertices are adjacent to it."""
    # This doesn't make sense — let me compute the degree sequence of the d3-subgraph
    # The d3-subgraph has vertices = V(G), edges = {u,v : d(u,v)=3 in G}
    # But this is NOT a subset of edges of G. Let me compute the full d3-subgraph.
    # Actually, the invariant used in the memory is "d3-subgraph connectivity"
    # meaning we look at the induced subgraph on Gamma_3(v) for each vertex v.

    dist0 = bfs_distances(nbrs, n, 0)
    d3_verts = [v for v in range(n) if dist0[v] == 3]

    # Edge count in induced subgraph on d3_verts
    edge_count = 0
    for v in d3_verts:
        for w in nbrs[v]:
            if dist0[w] == 3 and w > v:
                edge_count += 1
    return len(d3_verts), edge_count

# Use vertex 0 for both
n_d3_L, edges_d3_L = d3_degree_sequence(langer_nbrs, 63)
n_d3_DL, edges_d3_DL = d3_degree_sequence(dual_nbrs, 63)

print(f"Langer:      Gamma_3(0) has {n_d3_L} vertices, {edges_d3_L} internal edges")
print(f"Dual Langer: Gamma_3(0) has {n_d3_DL} vertices, {edges_d3_DL} internal edges")

# Try more vertices
print("\nD3 subgraph edge counts across vertices 0..4:")
for v in range(5):
    dist_v = bfs_distances(langer_nbrs, 63, v)
    d3v = [u for u in range(63) if dist_v[u] == 3]
    ec = sum(1 for u in d3v for w in langer_nbrs[u] if dist_v[w] == 3 and w > u)

    dist_v2 = bfs_distances(dual_nbrs, 63, v)
    d3v2 = [u for u in range(63) if dist_v2[u] == 3]
    ec2 = sum(1 for u in d3v2 for w in dual_nbrs[u] if dist_v2[w] == 3 and w > u)

    print(f"  v={v}: Langer d3-edges={ec}, Dual d3-edges={ec2}")

results['d3_analysis'] = {
    'langer_d3_from_0': {'verts': n_d3_L, 'edges': edges_d3_L},
    'dual_langer_d3_from_0': {'verts': n_d3_DL, 'edges': edges_d3_DL},
}

# ============================================================
# Characteristic polynomials (for cospectrality proof)
# ============================================================

print("\n" + "="*60)
print("CHARACTERISTIC POLYNOMIAL COEFFICIENTS (Langer vs Dual Langer)")
print("="*60)

# Compute charpoly from eigenvalues (product form)
def charpoly_coeffs(evals, tol=1e-6):
    """Compute coefficients of characteristic polynomial from eigenvalues.
    Uses Newton's identities / power sums."""
    n = len(evals)
    # p_k = sum of k-th powers of eigenvalues
    p = [sum(e**k for e in evals) for k in range(n+1)]
    # coefficients via Newton's identities
    c = [0.0] * (n+1)
    c[0] = 1.0
    for k in range(1, n+1):
        c[k] = -sum(c[k-j] * p[j] for j in range(1, k+1)) / k
    return [round(x) for x in c]

# This is too slow for n=63; use trace formula for first few coefficients
def trace_invariants(A, max_k=8):
    """Compute tr(A^k) for k=1..max_k."""
    Ak = np.eye(A.shape[0])
    traces = []
    for k in range(1, max_k+1):
        Ak = Ak @ A
        traces.append(float(np.trace(Ak)))
    return traces

traces_L = trace_invariants(A_L)
traces_DL = trace_invariants(A_DL)

print("tr(A^k) for k=1..8:")
print(f"{'k':>3}  {'Langer':>14}  {'Dual Langer':>14}  {'Match':>6}")
for k, (tL, tDL) in enumerate(zip(traces_L, traces_DL), 1):
    match = abs(tL - tDL) < 1e-4
    print(f"{k:>3}  {tL:>14.2f}  {tDL:>14.2f}  {'YES' if match else 'NO!':>6}")

results['trace_comparison'] = {
    'traces_langer': [float(t) for t in traces_L],
    'traces_dual_langer': [float(t) for t in traces_DL],
    'all_match': all(abs(a-b) < 1e-4 for a,b in zip(traces_L, traces_DL)),
}

# ============================================================
# Save JSON
# ============================================================

def make_json_safe(obj):
    if isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [make_json_safe(v) for v in obj]
    elif isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    elif isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    return obj

outfile = "/Users/robin/git/graphviz/spectral/langer_zhou_data.json"
with open(outfile, "w") as f:
    json.dump(make_json_safe(results), f, indent=2)
print(f"\nData saved to {outfile}")

# ============================================================
# Clayworth graph: Sab(G_2(2), C_3), 4032 vertices, cubic
# ============================================================
# G_2(2) has order 12096. The stabiliser C_3 has order 3.
# This gives 12096/3 = 4032 vertices, each vertex has degree 3.
# The graph is a {12,3} tiling of a genus-505 surface.
# With 4032 vertices, the adjacency matrix is 4032x4032 — feasible with numpy.

print("\n" + "="*60)
print("CLAYWORTH GRAPH (4032 vertices, cubic, Sab(G_2(2), C_3))")
print("="*60)
print("Building Clayworth graph via Langer adjacency lift...")

# The Clayworth graph is a 2-fold cover of the Langer graph (the 63-vertex
# 6-regular Langer graph is related to Clayworth via a covering map).
# Actually: Clayworth = Sab(G_2(2), C_3), and the cubic graph on 4032 vertices.
# We cannot construct it directly from the Langer formula without the group action.
# However, we can compute its eigenvalues using a theoretical approach:
#
# For a Sabidussi graph Sab(G, H) where G acts on cosets G/H,
# the spectrum is determined by the irreducible representations of G
# that appear in the permutation representation C[G/H].
#
# For Clayworth = Sab(G_2(2), C_3): G = G_2(2), |G| = 12096
# The permutation character on G/C_3 decomposes over the irreps of G_2(2).
# The Hecke algebra is the double-coset algebra C[C_3 \ G_2(2) / C_3].
# This is more complex than the 6-regular cases.
#
# Instead, we attempt a direct numpy computation.
# The Clayworth graph is also the incidence graph of a combinatorial structure.
# Since we cannot construct it purely from the Langer formula, we note:
#
# The Meinhold-72 graph (72 vertices, known from hecke catalog) is a quotient
# of Clayworth. The Clayworth graph fibers over the Langer graph with fiber C_2.
#
# For the paper, we compute what we can and state what remains open.

# Theoretical spectral data for Clayworth (from Hecke algebra theory):
# Clayworth is 3-regular on 4032 vertices.
# For a 3-regular Ramanujan graph, |lambda_2| <= 2*sqrt(2) ≈ 2.828
# The Biggs-Smith graph (102v, cubic) has spectrum {3^1, sqrt(2)^34, 0^17, -1^51, -sqrt(2)^?}
# Clayworth is expected to have similar algebraic structure from G_2(2) representation theory.

# We can estimate the spectrum from trace formulas if we knew the adjacency operator.
# Since we don't have the explicit adjacency list, we state theoretical results:

CLAYWORTH_THEORETICAL = {
    'n': 4032,
    'k': 3,
    'group': 'G_2(2)',
    'order_G': 12096,
    'order_H': 3,
    'genus': 505,
    'description': 'Sab(G_2(2), C_3): cubic graph on 4032 vertices, {12,3} tiling of genus-505 surface',
    'computation_status': 'adjacency_list_not_available',
    'ramanujan_bound': float(2 * np.sqrt(2)),
    'note': ('Clayworth cannot be drawn (too many vertices) and its explicit '
             'adjacency list requires GAP computation of G_2(2) coset table. '
             'The Fiedler value and full spectrum require future GAP work.'),
}
print("  n = 4032, k = 3, group = G_2(2), |G| = 12096, |H| = 3")
print("  Genus: 505 (from chi = V - E + F = 4032 - 6048 + 1008 = -1008 = 2 - 2*505)")
print("  Ramanujan bound: 2*sqrt(2) ≈ {:.6f}".format(2*np.sqrt(2)))
print("  Note: explicit adjacency list requires GAP coset enumeration of G_2(2)")
print("  The Meinhold-72 graph is a quotient of Clayworth.")

results['clayworth'] = CLAYWORTH_THEORETICAL

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Langer eigenvalues:      {groups_L}")
print(f"Dual Langer eigenvalues: {groups_DL}")
print(f"COSPECTRAL: {cospectral}")
print(f"Zhou eigenvalues:        {groups_Z}")
print(f"Langer Ramanujan: {ram_L}")
print(f"Dual Langer Ramanujan: {ram_DL}")
print(f"Zhou Ramanujan: {ram_Z}")
print(f"Clayworth: {CLAYWORTH_THEORETICAL['n']} vertices, cubic, Ramanujan bound = {CLAYWORTH_THEORETICAL['ramanujan_bound']:.4f}")
