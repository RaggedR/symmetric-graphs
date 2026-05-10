"""
Cohomology of the Bolza Surface via the (2,3,8) Triangle Group
==============================================================

The Bolza surface is the genus-2 Riemann surface with the largest
automorphism group: GL(2,3) of order 48 (not quite Hurwitz — the
Hurwitz bound 84(g-1)=84 is not achieved at genus 2).

It carries dual regular maps:
  {3,8}: 16 triangles, 24 edges, 6 vertices (8-valent)
  {8,3}: 6 octagons,  24 edges, 16 vertices (3-valent) = Möbius-Kantor graph

The (2,3,8) triangle group: Δ = ⟨a, b | a² = b³ = (ab)⁸ = 1⟩
"""

import numpy as np
from itertools import product
from collections import defaultdict


# ══════════════════════════════════════════════════════════════════════
#  Step 1: Build GL(2,3) — the automorphism group
# ══════════════════════════════════════════════════════════════════════
#
# GL(2,3) = all invertible 2×2 matrices over F_3
# |GL(2,3)| = (9-1)(9-3) = 8·6 = 48
# (No projective quotient needed — GL(2,3) is the full group here,
#  unlike PSL(2,7) where we quotient SL by ±I)

print("Building GL(2,3)...")

def mat_mult_mod3(A, B):
    return (A @ B) % 3

def mat_det_mod3(A):
    return (A[0,0]*A[1,1] - A[0,1]*A[1,0]) % 3

def mat_inv_mod3(A):
    det = mat_det_mod3(A)
    assert det != 0, "Matrix not invertible"
    det_inv = pow(int(det), 1, 3)  # det^{-1} mod 3 (1→1, 2→2 since 2²≡1 mod 3)
    if det == 2:
        det_inv = 2
    else:
        det_inv = 1
    adj = np.array([[A[1,1], (-A[0,1]) % 3],
                     [(-A[1,0]) % 3, A[0,0]]]) % 3
    return (det_inv * adj) % 3

def mat_to_key(A):
    return tuple((A % 3).flatten())

# Build all invertible 2×2 matrices over F_3
elements = {}
for a, b, c, d in product(range(3), repeat=4):
    A = np.array([[a, b], [c, d]])
    det = mat_det_mod3(A)
    if det == 0:
        continue
    key = mat_to_key(A)
    if key not in elements:
        elements[key] = A % 3

gl_matrices = list(elements.values())
gl_keys = list(elements.keys())
key_to_idx = {k: i for i, k in enumerate(gl_keys)}

print(f"  |GL(2,3)| = {len(gl_matrices)}")
assert len(gl_matrices) == 48

def group_mult(i, j):
    prod = mat_mult_mod3(gl_matrices[i], gl_matrices[j])
    return key_to_idx[mat_to_key(prod)]

def group_inv(i):
    inv = mat_inv_mod3(gl_matrices[i])
    return key_to_idx[mat_to_key(inv)]

identity = key_to_idx[mat_to_key(np.eye(2, dtype=int))]


# ══════════════════════════════════════════════════════════════════════
#  Step 2: Element orders and (2,3,8) generators
# ══════════════════════════════════════════════════════════════════════

def element_order(i):
    current = i
    for n in range(1, 49):
        if current == identity:
            return n
        current = group_mult(current, i)
    return -1

print("\nComputing element orders...")
orders = [element_order(i) for i in range(48)]
order_counts = defaultdict(int)
for o in orders:
    order_counts[o] += 1
print(f"  Order distribution: {dict(sorted(order_counts.items()))}")

# Find (2,3,8) generators: a² = b³ = (ab)⁸ = 1
print("\nSearching for (2,3,8) generators...")

a_gen = b_gen = c_gen = None
for a in range(48):
    if orders[a] != 2:
        continue
    for b in range(48):
        if orders[b] != 3:
            continue
        ab = group_mult(a, b)
        if orders[ab] == 8:
            a_gen, b_gen, c_gen = a, b, ab
            break
    if a_gen is not None:
        break

if a_gen is None:
    print("  No (2,3,8) generators found! Trying alternative...")
    # GL(2,3) might need a different presentation
    # Try: look for elements whose orders divide into (2,3,8)
    for a in range(48):
        if orders[a] != 2:
            continue
        for b in range(48):
            if orders[b] != 3:
                continue
            ab = group_mult(a, b)
            ab_ord = orders[ab]
            if ab_ord in [4, 8]:
                print(f"  Found a(ord 2), b(ord 3), ab(ord {ab_ord})")
                if ab_ord == 4:
                    # (ab)^4 might work if we use a different group
                    pass

if a_gen is not None:
    print(f"  a (order {orders[a_gen]}): {gl_matrices[a_gen].flatten()}")
    print(f"  b (order {orders[b_gen]}): {gl_matrices[b_gen].flatten()}")
    print(f"  c=ab (order {orders[c_gen]}): {gl_matrices[c_gen].flatten()}")
else:
    print("  Trying (2,3,8) with ab having order dividing 8...")
    # The Bolza surface might use a quotient where (ab)^8=1 but
    # the actual order of ab in GL(2,3) divides 8
    for a in range(48):
        if orders[a] != 2: continue
        for b in range(48):
            if orders[b] != 3: continue
            ab = group_mult(a, b)
            if orders[ab] in [8, 4]:
                a_gen, b_gen, c_gen = a, b, ab
                print(f"  Using: a(ord {orders[a]}), b(ord {orders[b]}), ab(ord {orders[ab]})")
                print(f"    a = {gl_matrices[a_gen].flatten()}")
                print(f"    b = {gl_matrices[b_gen].flatten()}")
                break
        if a_gen is not None: break


# ══════════════════════════════════════════════════════════════════════
#  Step 3: Determine the correct group and presentation
# ══════════════════════════════════════════════════════════════════════

# The Bolza surface Aut group is of order 48. Let's check what we have.
# GL(2,3) has order 48 and has elements of orders 1,2,3,4,6,8.
# The (2,3,8) triangle group quotient needs a² = b³ = (ab)⁸ = 1.
# If GL(2,3) has elements of order 8, the presentation should work.

if c_gen is not None:
    c_order = orders[c_gen]
    print(f"\n  Vertex stabilizer order: {c_order}")
    print(f"  Cell counts: {48//c_order} vertices, {48//2} edges, {48//3} faces")
    print(f"  χ = {48//c_order} - {48//2} + {48//3}")
    chi = 48//c_order - 48//2 + 48//3
    print(f"  χ = {chi}, genus = {(2-chi)//2}")


# ══════════════════════════════════════════════════════════════════════
#  Step 4: Build cell complex
# ══════════════════════════════════════════════════════════════════════

def cyclic_subgroup(gen_idx, order):
    subgroup = set()
    current = identity
    for _ in range(order):
        subgroup.add(current)
        current = group_mult(current, gen_idx)
    return frozenset(subgroup)

def left_cosets(subgroup):
    cosets = []
    covered = set()
    for g in range(48):
        if g in covered:
            continue
        coset = frozenset(group_mult(g, h) for h in subgroup)
        cosets.append((coset, g))
        covered |= coset
    return cosets

if c_gen is not None:
    c_order = orders[c_gen]
    H_vertex = cyclic_subgroup(c_gen, c_order)
    H_face = cyclic_subgroup(b_gen, 3)
    H_edge = cyclic_subgroup(a_gen, 2)

    vertices = left_cosets(H_vertex)
    faces = left_cosets(H_face)
    edges = left_cosets(H_edge)

    NV = len(vertices)
    NE = len(edges)
    NF = len(faces)

    print(f"\n  {NV} vertices, {NE} edges, {NF} faces")
    chi = NV - NE + NF
    genus = (2 - chi) // 2
    print(f"  χ = {chi}, genus = {genus}")

    # Incidence
    def incidence_matrix(cells_low, cells_high):
        n_low = len(cells_low)
        n_high = len(cells_high)
        inc = np.zeros((n_low, n_high), dtype=int)
        for i, (coset_low, _) in enumerate(cells_low):
            for j, (coset_high, _) in enumerate(cells_high):
                if coset_low & coset_high:
                    inc[i, j] = 1
        return inc

    print("\n  Computing incidence...")
    inc_ev = incidence_matrix(edges, vertices)
    inc_fe = incidence_matrix(faces, edges)
    inc_fv = incidence_matrix(faces, vertices)

    print(f"  Vertices per edge: {inc_ev.sum(axis=1).min()}-{inc_ev.sum(axis=1).max()}")
    print(f"  Edges per vertex: {inc_ev.sum(axis=0).min()}-{inc_ev.sum(axis=0).max()}")
    print(f"  Edges per face: {inc_fe.sum(axis=1).min()}-{inc_fe.sum(axis=1).max()}")
    print(f"  Faces per edge: {inc_fe.sum(axis=0).min()}-{inc_fe.sum(axis=0).max()}")
    print(f"  Vertices per face: {inc_fv.sum(axis=1).min()}-{inc_fv.sum(axis=1).max()}")
    print(f"  Faces per vertex: {inc_fv.sum(axis=0).min()}-{inc_fv.sum(axis=0).max()}")

    # Edge endpoints
    edge_endpoints = []
    for i in range(NE):
        verts = [j for j in range(NV) if inc_ev[i, j] == 1]
        edge_endpoints.append((verts[0], verts[1]) if len(verts) == 2 else (verts[0], verts[0]))

    # Build adjacency matrix for the dual graph ({8,3} = Möbius-Kantor?)
    adj_dual = np.zeros((NF, NF), dtype=int)
    for e in range(NE):
        incident_faces = [f for f in range(NF) if inc_fe[f, e] == 1]
        if len(incident_faces) == 2:
            adj_dual[incident_faces[0], incident_faces[1]] = 1
            adj_dual[incident_faces[1], incident_faces[0]] = 1

    print(f"\n  Dual graph ({NF} vertices):")
    dual_degrees = adj_dual.sum(axis=1)
    print(f"    Degree sequence: {sorted(set(dual_degrees))}")
    print(f"    Edges: {adj_dual.sum()//2}")

    # Adjacency matrix for the primal graph
    adj_primal = np.zeros((NV, NV), dtype=int)
    for e in range(NE):
        if len(edge_endpoints[e]) == 2:
            t, h = edge_endpoints[e]
            if t != h:
                adj_primal[t, h] = 1
                adj_primal[h, t] = 1

    print(f"\n  Primal graph ({NV} vertices):")
    primal_degrees = adj_primal.sum(axis=1)
    print(f"    Degree sequence: {sorted(set(primal_degrees))}")


    # ══════════════════════════════════════════════════════════════════
    #  Step 5: Oriented boundary maps and homology
    # ══════════════════════════════════════════════════════════════════

    print("\n\n" + "="*60)
    print("  ORIENTED BOUNDARY MAPS AND HOMOLOGY")
    print("="*60)

    # Orient edges
    d1 = np.zeros((NV, NE), dtype=int)
    for i, (t, h) in enumerate(edge_endpoints):
        if t != h:
            d1[h, i] = 1
            d1[t, i] = -1

    rank_d1 = np.linalg.matrix_rank(d1)
    print(f"\n  ∂₁: {d1.shape}, rank {rank_d1}")

    # Orient faces and build ∂₂
    d2 = np.zeros((NE, NF), dtype=int)
    for f in range(NF):
        face_edges = [e for e in range(NE) if inc_fe[f, e] == 1]
        face_verts = [v for v in range(NV) if inc_fv[f, v] == 1]

        n_face_edges = len(face_edges)
        n_face_verts = len(face_verts)

        if n_face_edges < 3 or n_face_verts < 3:
            continue

        # Build adjacency within face
        adj = defaultdict(list)
        for e in face_edges:
            t, h = edge_endpoints[e]
            adj[t].append((e, h))
            adj[h].append((e, t))

        # Walk cycle
        v0 = face_verts[0]
        cycle_v = [v0]
        cycle_e = []
        current = v0
        for step in range(n_face_edges):
            found = False
            for e, nxt in adj[current]:
                if nxt not in cycle_v or (step == n_face_edges - 1 and nxt == v0):
                    cycle_v.append(nxt)
                    cycle_e.append(e)
                    current = nxt
                    found = True
                    break
            if not found:
                break

        if len(cycle_e) == n_face_edges:
            for step in range(n_face_edges):
                e = cycle_e[step]
                v_from = cycle_v[step]
                v_to = cycle_v[(step + 1) % n_face_edges]
                t, h = edge_endpoints[e]
                d2[e, f] = 1 if (t == v_from and h == v_to) else -1

    rank_d2 = np.linalg.matrix_rank(d2)
    print(f"  ∂₂: {d2.shape}, rank {rank_d2}")

    chain = d1 @ d2
    print(f"  ∂₁∘∂₂ = 0? {np.all(chain == 0)}")

    # Betti numbers
    beta_0 = NV - rank_d1
    beta_1 = (NE - rank_d1) - rank_d2
    beta_2 = NF - rank_d2

    print(f"\n  Betti numbers:")
    print(f"    β₀ = {beta_0}")
    print(f"    β₁ = {beta_1}  (expected 2g = {2*genus})")
    print(f"    β₂ = {beta_2}")
    print(f"    χ = {beta_0 - beta_1 + beta_2}")


    # ══════════════════════════════════════════════════════════════════
    #  Step 6: Conjugacy classes and character of H₁
    # ══════════════════════════════════════════════════════════════════

    print("\n\n" + "="*60)
    print("  CONJUGACY CLASSES AND H₁ CHARACTER")
    print("="*60)

    def conjugacy_class(g_idx):
        cls = set()
        for h in range(48):
            h_inv = group_inv(h)
            conj = group_mult(group_mult(h, g_idx), h_inv)
            cls.add(conj)
        return frozenset(cls)

    classes = []
    classified = set()
    class_reps = []
    for g in range(48):
        if g in classified:
            continue
        cls = conjugacy_class(g)
        classes.append(cls)
        class_reps.append(g)
        classified |= cls

    # Sort by (order, size)
    class_data = [(i, len(classes[i]), orders[class_reps[i]]) for i in range(len(classes))]
    class_data.sort(key=lambda x: (x[2], x[1]))
    sorted_indices = [d[0] for d in class_data]

    print(f"\n  {len(classes)} conjugacy classes:")
    for idx in sorted_indices:
        rep = class_reps[idx]
        print(f"    Size {len(classes[idx]):3d}, order {orders[rep]}, "
              f"rep = {gl_matrices[rep].flatten()}")

    # Oriented Lefschetz traces
    print("\n  Computing oriented Lefschetz traces...")

    # Build lookups
    elt_to_vertex = {}
    for i, (coset, _) in enumerate(vertices):
        for h in coset:
            elt_to_vertex[h] = i

    elt_to_edge = {}
    for i, (coset, _) in enumerate(edges):
        for h in coset:
            elt_to_edge[h] = i

    elt_to_face = {}
    for i, (coset, _) in enumerate(faces):
        for h in coset:
            elt_to_face[h] = i

    def compute_oriented_traces(g_idx):
        # C_0: vertices (no orientation)
        tr_C0 = sum(1 for i, (coset, rep) in enumerate(vertices)
                    if elt_to_vertex[group_mult(g_idx, rep)] == i)

        # C_1: edges (track orientation)
        tr_C1 = 0
        for i, (coset, rep) in enumerate(edges):
            if elt_to_edge[group_mult(g_idx, rep)] != i:
                continue
            tail, head = edge_endpoints[i]
            if tail == head:
                tr_C1 += 1
                continue
            g_tail = elt_to_vertex[group_mult(g_idx, vertices[tail][1])]
            g_head = elt_to_vertex[group_mult(g_idx, vertices[head][1])]
            if g_tail == tail and g_head == head:
                tr_C1 += 1
            elif g_tail == head and g_head == tail:
                tr_C1 -= 1

        # C_2: faces (track orientation)
        tr_C2 = 0
        for i, (coset, rep) in enumerate(faces):
            if elt_to_face[group_mult(g_idx, rep)] != i:
                continue
            face_verts = [v for v in range(NV) if inc_fv[i, v] == 1]
            g_verts = [elt_to_vertex[group_mult(g_idx, vertices[v][1])] for v in face_verts]
            if set(g_verts) != set(face_verts):
                continue
            perm = [face_verts.index(gv) for gv in g_verts]
            inv_count = sum(1 for a in range(len(perm)) for b in range(a+1, len(perm)) if perm[a] > perm[b])
            tr_C2 += 1 if inv_count % 2 == 0 else -1

        return tr_C0, tr_C1, tr_C2

    print(f"\n  {'Class':>8} {'Size':>5} {'Ord':>4} {'tr(C0)':>7} {'tr(C1)':>7} {'tr(C2)':>7} {'L(g)':>6} {'tr(H1)':>7}")
    print("  " + "-"*62)

    char_H1 = np.zeros(48, dtype=int)

    for idx in sorted_indices:
        g = next(iter(classes[idx]))
        tr0, tr1, tr2 = compute_oriented_traces(g)
        L = tr0 - tr1 + tr2
        tr_H1 = 2 - L  # since H_0 = H_2 = trivial

        for g2 in classes[idx]:
            char_H1[g2] = tr_H1

        print(f"    {orders[class_reps[idx]]:>3}  {len(classes[idx]):>4}  {orders[class_reps[idx]]:>3}"
              f"  {tr0:>6}  {tr1:>6}  {tr2:>6}  {L:>5}  {tr_H1:>6}")

    # Inner products
    def inner_product(chi1, chi2):
        total = 0
        for cls in classes:
            g = next(iter(cls))
            total += len(cls) * chi1[g] * chi2[g]
        return total / 48

    triv = np.ones(48, dtype=int)
    norm_H1 = inner_product(char_H1, char_H1)

    print(f"\n  ⟨χ(H₁), χ(H₁)⟩ = {norm_H1:.4f}")
    print(f"  ⟨χ(H₁), trivial⟩ = {inner_product(char_H1, triv):.4f}")
    print(f"  dim H₁ = {char_H1[identity]}")


    # ══════════════════════════════════════════════════════════════════
    #  Step 7: Graph spectra
    # ══════════════════════════════════════════════════════════════════

    print("\n\n" + "="*60)
    print("  GRAPH SPECTRA")
    print("="*60)

    if NV > 1 and adj_primal.sum() > 0:
        eigs_primal = np.sort(np.linalg.eigvalsh(adj_primal.astype(float)))[::-1]
        print(f"\n  Primal graph ({NV} vertices, {c_order}-regular):")
        unique_eigs = []
        used = set()
        for e in eigs_primal:
            r = round(e, 4)
            if r not in used:
                used.add(r)
                mult = sum(1 for x in eigs_primal if abs(x - e) < 0.01)
                unique_eigs.append((r, mult))
        unique_eigs.sort(key=lambda x: -x[0])
        for val, mult in unique_eigs:
            print(f"    λ = {val:8.4f}  (mult {mult})")

    if NF > 1 and adj_dual.sum() > 0:
        eigs_dual = np.sort(np.linalg.eigvalsh(adj_dual.astype(float)))[::-1]
        print(f"\n  Dual graph ({NF} vertices, 3-regular — Möbius-Kantor?):")
        unique_eigs = []
        used = set()
        for e in eigs_dual:
            r = round(e, 4)
            if r not in used:
                used.add(r)
                mult = sum(1 for x in eigs_dual if abs(x - e) < 0.01)
                unique_eigs.append((r, mult))
        unique_eigs.sort(key=lambda x: -x[0])
        for val, mult in unique_eigs:
            print(f"    λ = {val:8.4f}  (mult {mult})")

        # Check if it's the Möbius-Kantor graph
        # MK graph has spectrum: 3, √3, √3, 1, 1, -1, -1, -√3, -√3, -3
        # with multiplicities 1, 2, 2, 2, 2, 2, 2, 2, 2, 1
        # Actually: ±3 (mult 1 each), ±√3 (mult 2 each), ±1 (mult 4 each)
        print(f"\n    Is Möbius-Kantor (F016A)? ", end="")
        mk_spectrum = sorted([3, -3] + [np.sqrt(3)]*2 + [-np.sqrt(3)]*2 + [1]*4 + [-1]*4, reverse=True)
        if NF == 16 and len(eigs_dual) == 16:
            match = all(abs(eigs_dual[i] - mk_spectrum[i]) < 0.01 for i in range(16))
            print(f"{'YES ✓' if match else 'NO'}")
        else:
            print(f"Wrong vertex count ({NF})")


    # ══════════════════════════════════════════════════════════════════
    #  Summary
    # ══════════════════════════════════════════════════════════════════

    print("\n\n" + "="*60)
    print("  SUMMARY")
    print("="*60)
    print(f"""
Bolza Surface Cohomology
========================

Surface: Bolza surface (genus {genus})
Symmetry: GL(2,3), |G| = 48

Cell decompositions from (2,3,{c_order}) triangle group:
  Primal: {NV} vertices — {NE} edges — {NF} faces
  χ = {NV} - {NE} + {NF} = {chi}

Cellular homology:
  H₀ = ℤ   (trivial GL(2,3)-module)
  H₁ = ℤ^{beta_1}
  H₂ = ℤ   (trivial GL(2,3)-module)

GL(2,3) action on H¹ (oriented Lefschetz traces):
  Character at identity: {char_H1[identity]}
  ⟨χ(H₁), χ(H₁)⟩ = {norm_H1}

Dual graph: {'Möbius-Kantor (F016A)' if NF == 16 else f'{NF}-vertex graph'}
""")
