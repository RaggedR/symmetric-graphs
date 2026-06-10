#!/usr/bin/env python3
"""Compute Z_7 orbit maps for the Langer and Klein graphs
using the SAME vertex ordering as the Lean definitions."""

from itertools import product as cart_product
import json

F2 = [0, 1]

# ============================================================
# LANGER GRAPH (63 vertices) — using Zorn product criterion
# Same formula as Q62.lean's langerAdjBool
# ============================================================

print("=" * 60)
print("  Langer graph: computing Z_7 orbits")
print("=" * 60)

# Generate Q(6,2) points in the same order as q62Indices
# q62Indices = (List.finRange 128).tail.filter (Q62form ∘ toBitVec = 0)
# toBitVec n i = (n / 2^i) % 2
def to_bit_vec(n):
    """Same as Lean toBitVec: bit i = (n / 2^i) % 2."""
    return tuple((n >> i) & 1 for i in range(7))

def q62form(x):
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

# Enumerate Q(6,2) points (indices 1..127 with Q=0)
q62_indices = [n for n in range(1, 128) if q62form(to_bit_vec(n)) == 0]
assert len(q62_indices) == 63, f"Expected 63 points, got {len(q62_indices)}"

# Build adjacency using Zorn product criterion
def dot(a, b):
    return sum(a[i]*b[i] for i in range(3)) % 2

def cross(a, b):
    return ((a[1]*b[2]+a[2]*b[1])%2, (a[2]*b[0]+a[0]*b[2])%2, (a[0]*b[1]+a[1]*b[0])%2)

def langer_adj(i, j):
    """Same formula as langerAdjBool in Q62.lean."""
    if i == j:
        return False
    x = to_bit_vec(q62_indices[i])
    y = to_bit_vec(q62_indices[j])
    alpha, gamma = x[3], y[3]
    # S1: αγ + a·d = 0
    if (alpha*gamma + x[0]*y[4] + x[1]*y[5] + x[2]*y[6]) % 2 != 0:
        return False
    # S2: αγ + b·c = 0
    if (alpha*gamma + x[4]*y[0] + x[5]*y[1] + x[6]*y[2]) % 2 != 0:
        return False
    # V1: γa + αc + b×d = 0
    a, b = (x[0],x[1],x[2]), (x[4],x[5],x[6])
    c, d = (y[0],y[1],y[2]), (y[4],y[5],y[6])
    bxd = cross(b, d)
    if (gamma*a[0] + alpha*c[0] + bxd[0]) % 2 != 0: return False
    if (gamma*a[1] + alpha*c[1] + bxd[1]) % 2 != 0: return False
    if (gamma*a[2] + alpha*c[2] + bxd[2]) % 2 != 0: return False
    # V2: γb + αd + a×c = 0
    axc = cross(a, c)
    if (gamma*b[0] + alpha*d[0] + axc[0]) % 2 != 0: return False
    if (gamma*b[1] + alpha*d[1] + axc[1]) % 2 != 0: return False
    if (gamma*b[2] + alpha*d[2] + axc[2]) % 2 != 0: return False
    return True

# Build adjacency lists
langer_adj_list = [[j for j in range(63) if langer_adj(i, j)] for i in range(63)]

# Find an automorphism of order 7
# Strategy: try permutations that map the graph to itself
# For GH(2,2), the Z_7 action comes from a cyclic subgroup of G_2(2)
# acting on Q(6,2) points. A natural Z_7 element: multiplication by
# a 7th root of unity in the split octonion algebra.
#
# Over F_2, we can try: the map (a,b) -> (a', b') where a' and b'
# are obtained by a cyclic permutation of coordinates.
# But this might not preserve Q(6,2).
#
# Simpler: brute-force search. Find a vertex v with neighbors N(v),
# then try mapping v -> each w and extending.

# Actually, let's use a cleaner approach: find a linear map on F_2^7
# that preserves Q and has order 7, then compute the induced permutation
# on the 63 points.

import numpy as np

def mat_mul_f2(A, v):
    """Matrix-vector multiply over F_2."""
    return tuple(sum(A[i][j] * v[j] for j in range(7)) % 2 for i in range(7))

def preserves_Q(A):
    """Check if linear map A preserves Q(6,2)."""
    for idx in q62_indices:
        v = to_bit_vec(idx)
        Av = mat_mul_f2(A, v)
        if q62form(Av) != 0:
            return False
    return True

def mat_order(A, max_order=100):
    """Compute order of matrix A over F_2."""
    B = [[1 if i==j else 0 for j in range(7)] for i in range(7)]
    for k in range(1, max_order+1):
        # B = B * A mod 2
        C = [[0]*7 for _ in range(7)]
        for i in range(7):
            for j in range(7):
                C[i][j] = sum(B[i][l] * A[l][j] for l in range(7)) % 2
        B = C
        if all(B[i][j] == (1 if i==j else 0) for i in range(7) for j in range(7)):
            return k
    return -1

def induced_perm(A):
    """Compute the permutation on q62_indices induced by linear map A."""
    idx_to_pos = {idx: pos for pos, idx in enumerate(q62_indices)}
    perm = []
    for pos, idx in enumerate(q62_indices):
        v = to_bit_vec(idx)
        Av = mat_mul_f2(A, v)
        # Find Av in q62_indices
        av_idx = sum(Av[i] * (2**i) for i in range(7))
        if av_idx not in idx_to_pos:
            return None  # Not a valid automorphism
        perm.append(idx_to_pos[av_idx])
    return perm

def is_automorphism(perm):
    """Check if perm is a graph automorphism of the Langer graph."""
    for i in range(63):
        for j in langer_adj_list[i]:
            if not langer_adj(perm[i], perm[j]):
                return False
    return True

# Try some matrices of order 7 over F_2
# A companion matrix for x^7 - 1 = (x-1)(x^3+x+1)(x^3+x^2+1) over F_2
# An element of order 7 in GL(7,F_2) can be a cyclic permutation matrix
# shifted by one step, but we need it to preserve Q.

# Let's try: the map that cyclically permutes coordinates but respects
# the Q pairing structure. Q pairs (0,4), (1,5), (2,6) with radical 3.
# A cyclic permutation of {0,1,2} and {4,5,6} simultaneously preserves Q.
# But that has order 3, not 7.

# For order 7: we need a 7x7 matrix over F_2 of order 7 that preserves Q.
# The group G_2(2) has elements of order 7. Let me try to find one.

# The primitive polynomial x^3 + x + 1 over F_2 has roots of order 7.
# We can embed F_8* (which has elements of order 7) into GL(7,F_2).

# Actually, let me try a different approach: find the automorphism by
# checking which permutations of the 63 vertices preserve adjacency
# and have order 7.

# Brute force is too slow for 63! permutations. Let me use the graph
# structure: pick a vertex, map it somewhere, extend greedily.

# A smarter approach: the Z_7 action partitions vertices into orbits of 7.
# If we find one orbit and the generator within that orbit, we can extend.

# Let me try: for each vertex v, compute its 7-neighborhood (vertices within
# distance 1, 2, 3). If two vertices have the same neighborhood structure,
# they're in the same orbit.

# Actually, for distance-regular graphs, the distance distribution from any
# vertex is: 1, 6, 24, 32. So all vertices look the same from a local
# perspective. We need global structure.

# Let me try a different approach: compute a matrix in GL(7,F_2) that
# preserves Q and has order 7, then compute its action on Q(6,2) points.

# The companion matrix of x^3 + x + 1 (irreducible over F_2, order 7):
# [0, 0, 1]
# [1, 0, 1]
# [0, 1, 0]
# This acts on F_2^3 with orbits of size 7 (plus the zero vector).

# To get a 7x7 matrix preserving Q, I can embed this in both the "a" and "b"
# blocks, with the radical coordinate x_3 determined.

# Q(x) = x_0*x_4 + x_1*x_5 + x_2*x_6 + x_3
# If I apply a 3x3 matrix M to (x_0,x_1,x_2) and (M^-T) to (x_4,x_5,x_6),
# then the sum x_0*x_4 + x_1*x_5 + x_2*x_6 is preserved (since a^T (M^-T b) = (M^-1 a)^T b).
# And x_3 = a·b is also preserved since a·b = (Ma)·(M^-T b).
# Wait, that's not right. We need a·b = (Ma)·(M^-T b)? Let me check:
# (Ma)^T (M^{-T}b) = a^T M^T M^{-T} b = a^T b. ✓

# So the 7x7 matrix is:
# [[M, 0, 0], [0, 1, 0], [0, 0, M^{-T}]]
# where M is a 3x3 matrix of order 7, arranged in the (0,1,2), (3), (4,5,6) blocks.

# Wait, the coordinate ordering is (x_0, x_1, x_2, x_3, x_4, x_5, x_6).
# The "a" block is (x_0, x_1, x_2) = coordinates 0,1,2.
# The "b" block is (x_4, x_5, x_6) = coordinates 4,5,6.
# The radical is x_3 = coordinate 3.

# So the 7x7 matrix should be:
# rows/cols indexed 0,1,2,3,4,5,6
# Block (0,1,2)×(0,1,2) = M
# Block (3)×(3) = 1
# Block (4,5,6)×(4,5,6) = M^{-T}
# All other blocks = 0

# Let me compute M^{-T} for M = companion of x^3 + x + 1.

M = [[0, 0, 1],
     [1, 0, 1],
     [0, 1, 0]]

# Compute M^{-1} over F_2
# det(M) = 0*0*0 + 0*1*0 + 1*1*1 - (1*0*0 + 0*1*0 + 0*1*1) = 1 - 0 = 1 (mod 2)
# So M is invertible.

# M^{-1} over F_2: using cofactors
# M^{-1} = (1/det) * adj(M)^T = adj(M)^T (since det = 1)

# Cofactor C_{ij} = (-1)^{i+j} * minor_{ij}
# Over F_2, (-1)^{i+j} = 1, so C_{ij} = minor_{ij}

# Minor(0,0) = det([[0,1],[1,0]]) = 0 - 1 = 1 (mod 2)
# Minor(0,1) = det([[1,1],[0,0]]) = 0 - 0 = 0
# Minor(0,2) = det([[1,0],[0,1]]) = 1 - 0 = 1
# Minor(1,0) = det([[0,1],[1,0]]) = 0 - 1 = 1
# Minor(1,1) = det([[0,1],[0,0]]) = 0 - 0 = 0
# Minor(1,2) = det([[0,0],[0,1]]) = 0 - 0 = 0
# Minor(2,0) = det([[0,1],[0,1]]) = 0 - 0 = 0
# Minor(2,1) = det([[0,1],[1,1]]) = 0 - 1 = 1
# Minor(2,2) = det([[0,0],[1,0]]) = 0 - 0 = 0

# adj(M) = [[1,1,0],[0,0,1],[1,0,0]]
# adj(M)^T = M^{-1} = [[1,0,1],[1,0,0],[0,1,0]]

M_inv = [[1, 1, 0],
         [0, 0, 1],
         [1, 0, 0]]

# Check M * M_inv = I mod 2
for i in range(3):
    for j in range(3):
        val = sum(M[i][k] * M_inv[k][j] for k in range(3)) % 2
        assert val == (1 if i == j else 0), f"M*M_inv[{i}][{j}] = {val}"

# M^{-T} = (M^{-1})^T
M_inv_T = [[M_inv[j][i] for j in range(3)] for i in range(3)]

# Build the 7x7 matrix
A = [[0]*7 for _ in range(7)]
for i in range(3):
    for j in range(3):
        A[i][j] = M[i][j]              # a-block
        A[i+4][j+4] = M_inv_T[i][j]    # b-block
A[3][3] = 1  # radical fixed

print("7x7 matrix A:")
for row in A:
    print(f"  {row}")

# Check: preserves Q?
assert preserves_Q(A), "Matrix does not preserve Q!"
print("✓ Preserves Q(6,2)")

# Check order
order = mat_order(A)
print(f"Order: {order}")
assert order == 7, f"Expected order 7, got {order}"
print("✓ Order 7")

# Compute induced permutation on q62_indices
perm = induced_perm(A)
assert perm is not None, "Not a valid permutation on Q(6,2) points"
assert len(set(perm)) == 63, "Permutation is not a bijection"

# Check it's a graph automorphism
assert is_automorphism(perm), "Not a graph automorphism!"
print("✓ Graph automorphism of Langer graph")

# Compute orbits
visited = [False] * 63
orbits = []
orbit_map = [0] * 63
for start in range(63):
    if visited[start]:
        continue
    orbit = []
    v = start
    for _ in range(7):
        orbit.append(v)
        visited[v] = True
        orbit_map[v] = len(orbits)
        v = perm[v]
    assert v == start, f"Orbit from {start} has wrong size"
    orbits.append(orbit)

print(f"\n{len(orbits)} orbits of size 7")
for i, orb in enumerate(orbits):
    print(f"  Orbit {i}: {orb}")

print(f"\nLean array (langerZ7Data):")
for start in range(0, 63, 10):
    chunk = orbit_map[start:start+10]
    print(f"  {', '.join(str(x) for x in chunk)},")

# Count quotient edges
q_edges = set()
for i in range(63):
    for j in langer_adj_list[i]:
        oi, oj = orbit_map[i], orbit_map[j]
        if oi != oj:
            q_edges.add((min(oi, oj), max(oi, oj)))
print(f"\nQuotient edges: {len(q_edges)}")
for e in sorted(q_edges):
    print(f"  {e}")

# ============================================================
# KLEIN GRAPH (56 vertices) — using CellularSurface edge data
# ============================================================

print("\n" + "=" * 60)
print("  Klein graph: computing Z_7 orbits")
print("=" * 60)

# Read edge data matching KleinSurface.lean
klein_src = [
    0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7, 7, 8, 8, 8,
    9, 9, 9, 10, 10, 11, 12, 13, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 18, 18,
    19, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29,
    30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
    45, 46, 47, 48
]

klein_tgt = [
    1, 6, 27, 2, 21, 3, 22, 4, 23, 5, 24, 6, 25, 26, 10, 11, 52, 11, 12, 53,
    12, 13, 54, 13, 55, 49, 50, 51, 16, 19, 40, 17, 20, 38, 18, 36, 19, 41, 20, 39,
    37, 35, 28, 43, 29, 44, 30, 45, 31, 46, 32, 47, 33, 48, 34, 42, 35, 50, 36, 52,
    37, 54, 38, 49, 39, 51, 40, 53, 41, 55, 42, 43, 44, 45, 46, 47, 48, 51, 53, 55,
    50, 52, 54, 49
]

assert len(klein_src) == 84 and len(klein_tgt) == 84

# Build adjacency
klein_adj = [set() for _ in range(56)]
for e in range(84):
    s, t = klein_src[e], klein_tgt[e]
    klein_adj[s].add(t)
    klein_adj[t].add(s)

# Verify 3-regular
for v in range(56):
    assert len(klein_adj[v]) == 3, f"Vertex {v} has degree {len(klein_adj[v])}"
print("✓ 3-regular, 84 edges")

# For Klein, the automorphism group is PGL(2,7), order 336.
# We need to find a Z_7 automorphism. Since the graph comes from PSL(2,7),
# a natural Z_7 element is the cyclic permutation of the 7 points of PG(1,7).
#
# But we don't know the labeling. Let me try a brute-force search:
# pick vertex 0, try mapping it to each vertex, extend using adjacency.
#
# A smarter approach: the Z_7 action on 56 vertices with 8 orbits of 7
# means the generator permutes vertices in cycles of length 7.
# If we find vertex 0's orbit, we can extend to all orbits.

# For a cubic graph, an automorphism must map edges to edges.
# Given v -> w, the 3 neighbors of v must map to the 3 neighbors of w.

def find_z7_automorphism():
    """Find an automorphism of order 7 of the Klein graph by backtracking."""
    n = 56

    def extend(perm, inv_perm):
        """Try to extend the partial permutation."""
        # Find an unmapped vertex that has a mapped neighbor
        for v in range(n):
            if v in perm:
                continue
            mapped_nbrs = [(u, perm[u]) for u in klein_adj[v] if u in perm]
            if not mapped_nbrs:
                continue
            # v must map to a vertex adjacent to perm[u] for all mapped neighbors u
            candidates = None
            for u, pu in mapped_nbrs:
                nbrs_pu = klein_adj[pu]
                if candidates is None:
                    candidates = set(w for w in nbrs_pu if w not in inv_perm)
                else:
                    candidates &= set(w for w in nbrs_pu if w not in inv_perm)
            if not candidates:
                return None
            # Try each candidate
            for w in candidates:
                perm2 = dict(perm)
                inv2 = dict(inv_perm)
                perm2[v] = w
                inv2[w] = v
                # Check consistency: all mapped neighbors of v must map to neighbors of w
                ok = True
                for u in klein_adj[v]:
                    if u in perm2:
                        if perm2[u] not in klein_adj[w]:
                            ok = False
                            break
                if not ok:
                    continue
                result = extend(perm2, inv2)
                if result is not None:
                    return result
            return None
        # All vertices mapped
        return perm

    # Try mapping vertex 0 to each vertex w
    for w in range(n):
        # Map 0 -> w, then try to extend
        if len(klein_adj[0]) != len(klein_adj[w]):
            continue
        # Try all 3! = 6 ways to map neighbors of 0 to neighbors of w
        import itertools
        nbrs_0 = sorted(klein_adj[0])
        nbrs_w = sorted(klein_adj[w])
        for perm_nbrs in itertools.permutations(nbrs_w):
            perm = {0: w}
            inv_perm = {w: 0}
            ok = True
            for i, u in enumerate(nbrs_0):
                v = perm_nbrs[i]
                if v in inv_perm:
                    ok = False
                    break
                perm[u] = v
                inv_perm[v] = u
            if not ok:
                continue
            # Check edge consistency among mapped vertices
            for u in perm:
                for v in klein_adj[u]:
                    if v in perm and perm[v] not in klein_adj[perm[u]]:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
            result = extend(perm, inv_perm)
            if result is not None:
                # Check it's really an automorphism
                full_perm = [result[v] for v in range(n)]
                is_auto = True
                for v in range(n):
                    for u in klein_adj[v]:
                        if full_perm[u] not in klein_adj[full_perm[v]]:
                            is_auto = False
                            break
                    if not is_auto:
                        break
                if not is_auto:
                    continue
                # Check order
                p = list(full_perm)
                pp = list(full_perm)
                for k in range(1, 8):
                    if all(pp[i] == i for i in range(n)):
                        if k == 7:
                            return full_perm
                        break
                    pp = [full_perm[pp[i]] for i in range(n)]
    return None

print("Searching for Z_7 automorphism...")
klein_perm = find_z7_automorphism()
if klein_perm is None:
    print("✗ No Z_7 automorphism found!")
else:
    print("✓ Found Z_7 automorphism")

    # Compute orbits
    visited = [False] * 56
    orbits = []
    orbit_map = [0] * 56
    for start in range(56):
        if visited[start]:
            continue
        orbit = []
        v = start
        for _ in range(7):
            orbit.append(v)
            visited[v] = True
            orbit_map[v] = len(orbits)
            v = klein_perm[v]
        assert v == start, f"Orbit from {start} has wrong size"
        orbits.append(orbit)

    print(f"\n{len(orbits)} orbits of size 7")
    for i, orb in enumerate(orbits):
        print(f"  Orbit {i}: {orb}")

    print(f"\nLean array (kleinZ7Data):")
    for start in range(0, 56, 10):
        chunk = orbit_map[start:start+10]
        print(f"  {', '.join(str(x) for x in chunk)},")

    # Count quotient edges
    q_edges = set()
    for v in range(56):
        for u in klein_adj[v]:
            oi, oj = orbit_map[v], orbit_map[u]
            if oi != oj:
                q_edges.add((min(oi, oj), max(oi, oj)))
    print(f"\nQuotient edges: {len(q_edges)}")
