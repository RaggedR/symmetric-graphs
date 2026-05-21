#!/usr/bin/env python3
"""Compute G₂(2) generators for the Langer graph on the 63 points of Q(6,2).

Uses the same vertex ordering as Lean's q62Indices (CSSCode/Q62.lean).
Outputs Lean 4 code for CSSCode/G2Action.lean:
  - Two generators as Equiv.Perm (Fin 63) lookup tables
  - BFS witness words mapping vertex 0 to each vertex
"""

from itertools import product as cart_product
from collections import deque

# ============================================================
# Q(6,2) point enumeration (same as Q62.lean)
# ============================================================

def to_bit_vec(n):
    """Same as Lean toBitVec: bit i = (n / 2^i) % 2."""
    return tuple((n >> i) & 1 for i in range(7))

def q62form(x):
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

q62_indices = [n for n in range(1, 128) if q62form(to_bit_vec(n)) == 0]
assert len(q62_indices) == 63

# ============================================================
# Langer adjacency (Zorn product criterion, same as Q62.lean)
# ============================================================

def langer_adj(i, j):
    if i == j:
        return False
    x = to_bit_vec(q62_indices[i])
    y = to_bit_vec(q62_indices[j])
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

adj_list = [[j for j in range(63) if langer_adj(i, j)] for i in range(63)]

# Verify 6-regular
for i in range(63):
    assert len(adj_list[i]) == 6, f"Vertex {i} has degree {len(adj_list[i])}"

# ============================================================
# Matrix utilities over F_2
# ============================================================

def mat_mul_f2(A, v):
    return tuple(sum(A[i][j] * v[j] for j in range(7)) % 2 for i in range(7))

def mat_mat_f2(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) % 2 for j in range(n)] for i in range(n)]

def mat_order(A, max_order=200):
    n = len(A)
    I = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    B = [row[:] for row in I]
    for k in range(1, max_order+1):
        B = mat_mat_f2(B, A)
        if B == I:
            return k
    return -1

def preserves_Q(A):
    for idx in q62_indices:
        v = to_bit_vec(idx)
        Av = mat_mul_f2(A, v)
        if q62form(Av) != 0:
            return False
    return True

def induced_perm(A):
    """Compute the permutation on Fin 63 indices induced by linear map A."""
    idx_to_pos = {idx: pos for pos, idx in enumerate(q62_indices)}
    perm = []
    for pos, idx in enumerate(q62_indices):
        v = to_bit_vec(idx)
        Av = mat_mul_f2(A, v)
        av_idx = sum(Av[i] * (2**i) for i in range(7))
        if av_idx not in idx_to_pos:
            return None
        perm.append(idx_to_pos[av_idx])
    return perm

def is_automorphism(perm):
    for i in range(63):
        for j in adj_list[i]:
            if perm[j] not in adj_list[perm[i]]:
                return False
    return True

def compose_perm(p, q):
    """Compose permutations: (p ∘ q)(i) = p(q(i))."""
    return [p[q[i]] for i in range(len(p))]

def invert_perm(p):
    inv = [0] * len(p)
    for i, v in enumerate(p):
        inv[v] = i
    return inv

def perm_order(p):
    q = list(p)
    for k in range(1, 200):
        if all(q[i] == i for i in range(len(q))):
            return k
        q = compose_perm(p, q)
    return -1

# ============================================================
# Generator 1: order-7 (from compute_orbits.py)
# ============================================================

# Companion matrix of x^3 + x + 1 in the a-block, M^{-T} in b-block
M = [[0, 0, 1],
     [1, 0, 1],
     [0, 1, 0]]

M_inv = [[1, 1, 0],
         [0, 0, 1],
         [1, 0, 0]]

M_inv_T = [[M_inv[j][i] for j in range(3)] for i in range(3)]

A1 = [[0]*7 for _ in range(7)]
for i in range(3):
    for j in range(3):
        A1[i][j] = M[i][j]
        A1[i+4][j+4] = M_inv_T[i][j]
A1[3][3] = 1

gen1 = induced_perm(A1)
assert gen1 is not None
assert is_automorphism(gen1)
assert perm_order(gen1) == 7
print(f"Generator 1: order {perm_order(gen1)}")

# ============================================================
# Generator 2: search for a graph automorphism that with gen1
# generates G₂(2) (order 12096) on the 63 points
# ============================================================

print("Searching for generator 2...")

# Transvections of Q(6,2) generate O(7,2), NOT G₂(2).
# G₂(2) = Aut(GH(2,2)) = Aut(Langer graph).
# So we search for graph automorphisms directly by backtracking.

def orbit_under_gens(gens, start=0):
    """Compute the orbit of `start` under a list of permutations."""
    seen = {start}
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for g in gens:
            w = g[v]
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return seen

def group_order_by_enumeration(gens, n=63, max_size=15000):
    """Compute |<gens>| by enumerating all group elements."""
    identity = list(range(n))
    elements = {tuple(identity)}
    queue = deque([identity])
    while queue:
        g = queue.popleft()
        for gen in gens:
            for h in [compose_perm(gen, g), compose_perm(g, gen)]:
                th = tuple(h)
                if th not in elements:
                    elements.add(th)
                    queue.append(h)
                    if len(elements) > max_size:
                        return len(elements)
    return len(elements)

def find_graph_automorphism(target_vertex):
    """Find a graph automorphism of the Langer graph mapping 0 -> target_vertex.
    Uses backtracking with the distance-regular structure."""
    import itertools

    nbrs_0 = sorted(adj_list[0])
    nbrs_t = sorted(adj_list[target_vertex])

    # For each way to map neighbors of 0 to neighbors of target:
    for perm_nbrs in itertools.permutations(nbrs_t):
        perm = {0: target_vertex}
        inv_perm = {target_vertex: 0}
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

        # Check edge consistency among already-mapped vertices
        for u in list(perm.keys()):
            for w in adj_list[u]:
                if w in perm and perm[w] not in adj_list[perm[u]]:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue

        # Try to extend greedily
        result = _extend_automorphism(perm, inv_perm)
        if result is not None:
            return result
    return None

def _extend_automorphism(perm, inv_perm):
    """Extend partial automorphism by constraint propagation + backtracking."""
    n = 63
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if v in perm:
                continue
            # Find mapped neighbors of v
            mapped_nbrs = [(u, perm[u]) for u in adj_list[v] if u in perm]
            if not mapped_nbrs:
                continue

            # Candidates: must be adjacent to images of all mapped neighbors
            # and not already used
            candidates = None
            for u, pu in mapped_nbrs:
                nbrs_pu = set(adj_list[pu])
                avail = set(w for w in nbrs_pu if w not in inv_perm)
                if candidates is None:
                    candidates = avail
                else:
                    candidates &= avail

            if candidates is None:
                continue
            if len(candidates) == 0:
                return None  # dead end
            if len(candidates) == 1:
                w = candidates.pop()
                # Check full consistency
                ok = True
                for u in adj_list[v]:
                    if u in perm and perm[u] not in adj_list[w]:
                        ok = False
                        break
                for u in range(n):
                    if u not in perm and u != v:
                        continue
                    if u == v:
                        continue
                    if perm.get(u) == w:
                        ok = False
                        break
                if not ok:
                    return None
                perm[v] = w
                inv_perm[w] = v
                changed = True

    # Check if complete
    if len(perm) == n:
        # Verify it's a valid automorphism
        full_perm = [perm[i] for i in range(n)]
        if is_automorphism(full_perm):
            return full_perm
        return None

    # Backtrack on the first unmapped vertex with constraints
    for v in range(n):
        if v in perm:
            continue
        mapped_nbrs = [(u, perm[u]) for u in adj_list[v] if u in perm]
        if not mapped_nbrs:
            continue

        candidates = None
        for u, pu in mapped_nbrs:
            avail = set(w for w in adj_list[pu] if w not in inv_perm)
            if candidates is None:
                candidates = avail
            else:
                candidates &= avail

        if not candidates:
            return None

        for w in sorted(candidates):
            ok = True
            for u in adj_list[v]:
                if u in perm and perm[u] not in adj_list[w]:
                    ok = False
                    break
            if not ok:
                continue

            perm2 = dict(perm)
            inv2 = dict(inv_perm)
            perm2[v] = w
            inv2[w] = v
            result = _extend_automorphism(perm2, inv2)
            if result is not None:
                return result
        return None

    return None

# Strategy: map 0 to a vertex in a different gen1-orbit.
# gen1 orbits: 0 is in orbit 0 = {0,1,3,2,5,6,4} (from compute_orbits.py).
# We want gen2 to map 0 to a vertex NOT in orbit 0.

# First compute gen1 orbits
visited = [False] * 63
gen1_orbits = []
for start in range(63):
    if visited[start]:
        continue
    orbit = []
    v = start
    for _ in range(7):
        orbit.append(v)
        visited[v] = True
        v = gen1[v]
    gen1_orbits.append(orbit)

print(f"Gen1 orbits: {[orb[0] for orb in gen1_orbits]}")

# Try mapping 0 to each vertex not in orbit 0
found = False
for target in [7, 8, 9, 10, 11, 12, 13, 14]:  # representatives of other orbits
    print(f"  Trying 0 -> {target}...")
    auto = find_graph_automorphism(target)
    if auto is None:
        continue

    print(f"  Found automorphism! Order: {perm_order(auto)}")

    # Check if gen1 + auto generate full G₂(2)
    orbit = orbit_under_gens([gen1, auto, invert_perm(gen1), invert_perm(auto)])
    if len(orbit) < 63:
        print(f"  Orbit size {len(orbit)}, not transitive")
        continue

    print(f"  Transitive! Computing group order...")
    order = group_order_by_enumeration([gen1, auto, invert_perm(gen1), invert_perm(auto)])
    print(f"  Group order: {order}")

    if order == 12096:
        print("  *** Full G₂(2) generated! ***")
        best_gen2 = auto
        found = True
        break
    else:
        print(f"  (subgroup of order {order}, trying next...)")

if not found:
    print("ERROR: Could not find gen2!")
    exit(1)

gen2 = best_gen2
gen1_inv = invert_perm(gen1)
gen2_inv = invert_perm(gen2)

print(f"\nGenerator 1: order {perm_order(gen1)}")
print(f"  {gen1}")
print(f"Generator 2: order {perm_order(gen2)}")
print(f"  {gen2}")

# ============================================================
# BFS witness words: for each vertex v, find shortest word in
# {gen1, gen1^-1, gen2, gen2^-1} mapping 0 -> v
# ============================================================

print("\nComputing BFS witness words...")

# Generators indexed 0..3: gen1, gen1^-1, gen2, gen2^-1
gens = [gen1, gen1_inv, gen2, gen2_inv]

# BFS from vertex 0
witness_word = {0: []}  # vertex -> list of generator indices
queue = deque([0])

while len(witness_word) < 63:
    v = queue.popleft()
    for gi, g in enumerate(gens):
        w = g[v]
        if w not in witness_word:
            witness_word[w] = witness_word[v] + [gi]
            queue.append(w)

# Verify all vertices reached
assert len(witness_word) == 63, f"Only reached {len(witness_word)} vertices"

# Verify each word is correct
def apply_word(word, start=0):
    v = start
    for gi in word:
        v = gens[gi][v]
    return v

for v in range(63):
    assert apply_word(witness_word[v]) == v, f"Word for vertex {v} is incorrect"

max_len = max(len(w) for w in witness_word.values())
print(f"All 63 vertices reached. Max word length: {max_len}")

# ============================================================
# Output Lean 4 code
# ============================================================

print("\n" + "=" * 60)
print("  Lean 4 output for CSSCode/G2Action.lean")
print("=" * 60)

def fmt_perm(p, name):
    """Format a permutation as a Lean Array (Fin 63)."""
    lines = [f"private def {name} : Array (Fin 63) := #["]
    for start in range(0, 63, 9):
        chunk = p[start:start+9]
        line = "  " + ", ".join(f"⟨{v}, by omega⟩" for v in chunk)
        if start + 9 < 63:
            line += ","
        lines.append(line)
    lines.append("]")
    return "\n".join(lines)

print()
print(fmt_perm(gen1, "g2gen1Fwd"))
print()
print(fmt_perm(gen1_inv, "g2gen1Inv"))
print()
print(fmt_perm(gen2, "g2gen2Fwd"))
print()
print(fmt_perm(gen2_inv, "g2gen2Inv"))

# Output witness words
print()
print("/-- BFS witness words: witnessWord v is a list of generator indices")
print("    (0=σ₁, 1=σ₁⁻¹, 2=σ₂, 3=σ₂⁻¹) such that applying them maps 0 → v. -/")
print("def witnessWord : Fin 63 → List (Fin 4)")

for v in range(63):
    word = witness_word[v]
    if not word:
        word_str = "[]"
    else:
        word_str = "[" + ", ".join(str(gi) for gi in word) + "]"
    sep = "\n" if v == 0 else ""
    prefix = f"  | ⟨{v}, _⟩ => {word_str}"
    print(prefix)
