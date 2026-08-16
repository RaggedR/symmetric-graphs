#!/usr/bin/env python3
"""
Spectral analysis of voltage graphs: cyclic covers and cubic residue graphs.

The voltage graph (derived graph) construction takes a base graph whose edges
carry labels (voltages) in a cyclic group Z_m and lifts it to a graph on
m × (base vertices) vertices.

We study:
  - CR(p): the cubic residue graph — a Z_p-cover of K₂ with voltages = {1, g, g²}
    where g is a primitive cube root of unity mod p (p ≡ 1 mod 3).
  - GP(n, k): generalised Petersen graphs as Z_n-covers of a theta-with-loops base.
  - PG(2,3) projection (Foster F26A): same as CR(13).

Spectral theory: for a Z_m-cover of a graph G, the adjacency eigenvalues are
obtained by evaluating the "voltage Laplacian" for each character χ_j(g) = ω^(jg)
where ω = e^(2πi/m). For the K₂ base with voltage set V, eigenvalue for character j is:
    λ_j = Σ_{v ∈ V} 2·cos(2πjv/m)
(the factor of 2 comes from counting both v and m-v as voltages in the undirected lift).
"""

import numpy as np
from math import gcd
from collections import Counter


# ---------------------------------------------------------------------------
# Primitive cube root of unity mod p
# ---------------------------------------------------------------------------

def primitive_cube_root_mod(p):
    """Return g with g³ ≡ 1 (mod p) and g ≠ 1, for prime p ≡ 1 (mod 3).

    Such g exists iff 3 | (p-1). We search g in 2..p-1.
    """
    if (p - 1) % 3 != 0:
        raise ValueError(f"p = {p} is not ≡ 1 (mod 3), no cube root of unity exists")
    for g in range(2, p):
        if pow(g, 3, p) == 1:
            return g
    raise ValueError(f"No primitive cube root of unity found mod {p}")


def cubic_residues_mod(p):
    """Return {1, g, g²} for the primitive cube root g mod p."""
    g = primitive_cube_root_mod(p)
    return sorted({1, g, g * g % p})


# ---------------------------------------------------------------------------
# Voltage graph lifting
# ---------------------------------------------------------------------------

def lift_voltage_graph(n_base, edges, m):
    """Build the adjacency matrix of the derived graph from a voltage graph.

    Parameters
    ----------
    n_base : int
        Number of vertices in the base graph (vertices indexed 1..n_base).
    edges : list of (i, j, v)
        Each entry is a directed edge from vertex i to vertex j with voltage
        v ∈ Z_m.  For an undirected edge between i and j with voltage v,
        include BOTH (i, j, v) and (j, i, -v mod m).
        For a self-loop at i with voltage v, include (i, i, v).
    m : int
        Order of the voltage group Z_m.

    Returns
    -------
    A : ndarray of shape (n_base * m, n_base * m)
        Adjacency matrix. Vertex (i, g) for i=1..n_base, g=0..m-1 has
        index (i-1)*m + g.
    """
    N = n_base * m
    A = np.zeros((N, N), dtype=float)

    def idx(i, g):
        # i is 1-indexed base vertex, g is 0..m-1
        return (i - 1) * m + g % m

    for (i, j, v) in edges:
        for g in range(m):
            u = idx(i, g)
            w = idx(j, (g + v) % m)
            A[u, w] = 1.0

    return A


def undirected_edge_list(i, j, v, m):
    """Return both directed copies of an undirected edge with voltage v."""
    return [(i, j, v % m), (j, i, (-v) % m)]


def self_loop_edge(i, v, m):
    """Return the edge list for a self-loop at vertex i with voltage v.

    A self-loop (i, i, v) in the voltage graph gives edges
    (i, g) — (i, (g+v) mod m) for all g. We also need the reverse
    (i, g) — (i, (g-v) mod m), since the graph is undirected.
    These are distinct unless 2v ≡ 0 (mod m).
    """
    edges = []
    for g in range(0, 1):  # template; lift_voltage_graph iterates over all g
        pass
    # Include both directions
    return [(i, i, v % m), (i, i, (-v) % m)]


# ---------------------------------------------------------------------------
# Specific graph constructions
# ---------------------------------------------------------------------------

def build_CR(p):
    """Build the cubic residue graph CR(p) for prime p ≡ 1 (mod 3).

    CR(p) is the Cayley graph on Z_p with connection set = cubic residues mod p.
    Equivalently it is a Z_p-cover of K₂ (2 vertices, 3 edges between them)
    where the three edges carry voltages {1, g, g²}.

    We construct it directly as a Cayley graph on Z_p.
    """
    S = cubic_residues_mod(p)
    # Cayley graph on Z_p with generators S ∪ {p-s : s ∈ S}
    # (since the graph is undirected, -s ≡ p-s is also a generator)
    connection_set = set()
    for s in S:
        connection_set.add(s)
        connection_set.add((-s) % p)

    N = p
    A = np.zeros((N, N), dtype=float)
    for v in range(N):
        for s in connection_set:
            w = (v + s) % p
            A[v, w] = 1.0
    return A, sorted(S)


def build_CR_via_voltage(p):
    """Build CR(p) as a Z_p-cover of K₂ via voltage graph construction.

    Base graph: 2 vertices (1 and 2), three undirected edges between them
    with voltages {1, g, g²}.  This gives a 3-regular bipartite graph
    with parts {(1,g)} and {(2,g)} for g ∈ Z_p.

    Note: the Cayley-graph construction (build_CR) gives the same spectrum
    but is not bipartite. This voltage lift IS bipartite.
    """
    S = cubic_residues_mod(p)
    edges = []
    for v in S:
        edges += undirected_edge_list(1, 2, v, p)
    A = lift_voltage_graph(2, edges, p)
    return A, sorted(S)


def build_generalised_petersen(n, k):
    """Build the generalised Petersen graph GP(n, k).

    GP(n, k) has 2n vertices: outer cycle O = {0..n-1} and inner star
    I = {n..2n-1}.  Edges:
      - outer cycle: i — (i+1) mod n
      - spokes: i — (n+i)
      - inner star: (n+i) — (n+(i+k) mod n)
    """
    N = 2 * n
    A = np.zeros((N, N), dtype=float)
    for i in range(n):
        # outer cycle
        j = (i + 1) % n
        A[i, j] = A[j, i] = 1.0
        # spoke
        A[i, n + i] = A[n + i, i] = 1.0
        # inner star
        j2 = (i + k) % n
        A[n + i, n + j2] = A[n + j2, n + i] = 1.0
    return A


def build_petersen_via_voltage():
    """Build the Petersen graph as a Z₅-cover of a theta-with-loops base.

    Base graph: 2 vertices.
      - self-loop on vertex 1 with voltage 1  (outer 5-cycle)
      - self-loop on vertex 2 with voltage 2  (inner pentagram)
      - edge 1-2 with voltage 0               (spokes)
    The self-loop (i,i,v) contributes edges (i,g)-(i,(g+v) mod m).
    """
    m = 5
    edges = []
    # Self-loop at 1 with voltage 1 (and its reverse, voltage -1 ≡ 4)
    edges += self_loop_edge(1, 1, m)
    # Self-loop at 2 with voltage 2 (and its reverse, voltage -2 ≡ 3)
    edges += self_loop_edge(2, 2, m)
    # Undirected edge 1-2 with voltage 0
    edges += undirected_edge_list(1, 2, 0, m)
    A = lift_voltage_graph(2, edges, m)
    return A


def build_cube_via_voltage():
    """Build the 3-cube Q₃ as a Z₄-cover of a theta-with-loops base.

    GP(4,1) base: self-loop at 1 (voltage 1), self-loop at 2 (voltage 1),
    edge 1-2 (voltage 0).
    """
    m = 4
    edges = []
    edges += self_loop_edge(1, 1, m)
    edges += self_loop_edge(2, 1, m)
    edges += undirected_edge_list(1, 2, 0, m)
    A = lift_voltage_graph(2, edges, m)
    return A


def build_mobius_kantor_via_voltage():
    """Build the Möbius-Kantor graph GP(8,3) as a Z₈-cover.

    Base: self-loop at 1 (voltage 1), self-loop at 2 (voltage 3),
    edge 1-2 (voltage 0).
    """
    m = 8
    edges = []
    edges += self_loop_edge(1, 1, m)
    edges += self_loop_edge(2, 3, m)
    edges += undirected_edge_list(1, 2, 0, m)
    A = lift_voltage_graph(2, edges, m)
    return A


# ---------------------------------------------------------------------------
# Spectral computations
# ---------------------------------------------------------------------------

def compute_spectrum(A):
    """Return sorted adjacency eigenvalues (descending) for symmetric matrix A."""
    eigs = np.linalg.eigvalsh(A)
    return np.sort(eigs)[::-1]


def spectrum_with_multiplicities(eigs, tol=1e-8):
    """Round eigenvalues and collect multiplicities.

    Returns list of (value, multiplicity) sorted descending by value.
    Values are rounded to 6 decimal places for display.
    """
    rounded = np.round(eigs, 6)
    counts = Counter(rounded)
    return sorted(counts.items(), key=lambda x: -x[0])


def spectral_gap(eigs):
    """Return λ₁ - λ₂ (spectral gap)."""
    if len(eigs) < 2:
        return 0.0
    return eigs[0] - eigs[1]


def fiedler_value(A):
    """Algebraic connectivity: smallest nonzero Laplacian eigenvalue."""
    n = A.shape[0]
    d = A.sum(axis=1)
    L = np.diag(d) - A
    lap_eigs = np.sort(np.linalg.eigvalsh(L))
    # Find first eigenvalue > tol
    for e in lap_eigs:
        if e > 1e-8:
            return e
    return 0.0


def is_ramanujan(eigs, degree):
    """Check Ramanujan: all non-trivial eigenvalues satisfy |λ| ≤ 2√(d-1)."""
    bound = 2.0 * np.sqrt(degree - 1)
    trivial = eigs[0]  # largest eigenvalue = degree
    for e in eigs:
        if abs(e) < trivial - 1e-8:  # non-trivial
            if abs(e) > bound + 1e-8:
                return False, bound
    return True, bound


def compute_girth(A):
    """Compute girth via BFS from each vertex. Returns ∞ if no cycle."""
    n = A.shape[0]
    adj = [list(np.where(A[v] > 0.5)[0]) for v in range(n)]
    girth = float('inf')
    for start in range(n):
        # BFS; track (current, parent, distance)
        dist = [-1] * n
        dist[start] = 0
        parent = [-1] * n
        queue = [start]
        head = 0
        found = False
        while head < len(queue) and not found:
            v = queue[head]; head += 1
            for w in adj[v]:
                if dist[w] == -1:
                    dist[w] = dist[v] + 1
                    parent[w] = v
                    queue.append(w)
                elif parent[v] != w:
                    # Found a cycle of length dist[v] + dist[w] + 1
                    cycle_len = dist[v] + dist[w] + 1
                    girth = min(girth, cycle_len)
                    found = True
                    break
        if girth <= 3:
            break  # can't get smaller
    return girth


def theoretical_cr_eigenvalues(p, num_eigs=20):
    """Compute CR(p) eigenvalues analytically from character sums.

    For the Cayley graph on Z_p with connection set S = {s, -s : s ∈ cubic_residues},
    eigenvalue for character j is: Σ_{s ∈ S_full} ω^(j*s) where ω = e^(2πi/p).
    """
    S = cubic_residues_mod(p)
    # Full symmetric connection set
    S_full = set()
    for s in S:
        S_full.add(s)
        S_full.add((-s) % p)

    omega = 2 * np.pi / p
    eigs = []
    for j in range(p):
        lam = sum(np.cos(omega * j * s) for s in S_full)
        eigs.append(lam)
    return np.sort(eigs)[::-1]


# ---------------------------------------------------------------------------
# Full graph analysis
# ---------------------------------------------------------------------------

def analyse_graph(name, A, degree, description=""):
    """Compute and print full spectral summary for a graph."""
    n = A.shape[0]
    eigs = compute_spectrum(A)
    spec = spectrum_with_multiplicities(eigs)
    gap = spectral_gap(eigs)
    fiedler = fiedler_value(A)
    ram, ram_bound = is_ramanujan(eigs, degree)
    girth = compute_girth(A)
    lambda2 = eigs[1] if len(eigs) > 1 else None

    print(f"\n{'='*65}")
    print(f"  {name}")
    if description:
        print(f"  {description}")
    print(f"  {n} vertices, {degree}-regular")
    print(f"{'='*65}")

    print(f"\n  Adjacency spectrum:")
    print(f"  {'Eigenvalue':>12}  {'Multiplicity':>12}")
    print(f"  {'─'*12}  {'─'*12}")
    for val, mult in spec:
        print(f"  {val:>12.6f}  {mult:>12}")

    print(f"\n  Key invariants:")
    print(f"    λ₁ (largest)     = {eigs[0]:.6f}")
    print(f"    λ₂               = {lambda2:.6f}")
    print(f"    λ_min            = {eigs[-1]:.6f}")
    print(f"    Spectral gap     = {gap:.6f}")
    print(f"    Fiedler value    = {fiedler:.6f}")
    print(f"    Girth            = {girth}")
    print(f"    Ramanujan bound  = 2√(d-1) = {ram_bound:.6f}")
    print(f"    Is Ramanujan?    {'YES' if ram else 'NO'}")

    return {
        'name': name,
        'n': n,
        'degree': degree,
        'girth': girth,
        'lambda1': float(eigs[0]),
        'lambda2': float(lambda2),
        'lambda_min': float(eigs[-1]),
        'spectral_gap': float(gap),
        'fiedler': float(fiedler),
        'ramanujan_bound': float(ram_bound),
        'is_ramanujan': ram,
        'spectrum': [(float(v), int(m)) for v, m in spec],
    }


# ---------------------------------------------------------------------------
# Analytical eigenvalue formula for K₂-cover with voltage set V
# ---------------------------------------------------------------------------

def k2_cover_eigenvalues_analytical(p, V):
    """Eigenvalues of Z_p-cover of K₂ with voltages V ⊆ Z_p.

    For character χ_j(g) = e^(2πi·j·g/p), the 2×2 voltage Laplacian is:
        M_j = [[0, f(j)], [conj(f(j)), 0]]
    where f(j) = Σ_{v ∈ V} ω^(j·v), ω = e^(2πi/p).

    Eigenvalues of M_j are ±|f(j)|. As j ranges over Z_p, we get 2p eigenvalues.
    """
    omega = 2 * np.pi / p
    all_eigs = []
    for j in range(p):
        fj = sum(np.exp(1j * omega * j * v) for v in V)
        magnitude = abs(fj)
        all_eigs.append(magnitude)
        all_eigs.append(-magnitude)
    return np.sort(all_eigs)[::-1]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    results = []

    # -----------------------------------------------------------------------
    # 1. Generalised Petersen graphs as voltage covers
    # -----------------------------------------------------------------------

    print("\n" + "#"*65)
    print("  SECTION 1: Generalised Petersen Graphs as Voltage Covers")
    print("#"*65)

    # Petersen graph GP(5,2)
    A_pet = build_generalised_petersen(5, 2)
    results.append(analyse_graph(
        "Petersen graph GP(5,2)",
        A_pet, 3,
        "Z₅-cover of theta-with-loops (voltages 1, 2, 0)"
    ))

    # Cube GP(4,1)
    A_cube = build_generalised_petersen(4, 1)
    results.append(analyse_graph(
        "3-Cube Q₃ = GP(4,1)",
        A_cube, 3,
        "Z₄-cover of theta-with-loops (voltages 1, 1, 0)"
    ))

    # Möbius-Kantor GP(8,3)
    A_mk = build_generalised_petersen(8, 3)
    results.append(analyse_graph(
        "Möbius-Kantor graph GP(8,3)",
        A_mk, 3,
        "Z₈-cover of theta-with-loops (voltages 1, 3, 0)"
    ))

    # -----------------------------------------------------------------------
    # 2. Cubic residue graphs CR(p)
    # -----------------------------------------------------------------------

    print("\n" + "#"*65)
    print("  SECTION 2: Cubic Residue Graphs CR(p)")
    print("#"*65)

    cr_primes = [7, 13, 19, 31, 37, 43]
    cr_names = {
        7: "CR(7) = Heawood graph",
        13: "CR(13) = Foster F26A / PG(2,3) projection",
        19: "CR(19)",
        31: "CR(31)",
        37: "CR(37)",
        43: "CR(43)",
    }

    for p in cr_primes:
        A_cr, voltages = build_CR(p)
        name = cr_names.get(p, f"CR({p})")
        desc = f"Cayley graph on Z_{p}, connection set = cubic residues {voltages} ∪ negatives"
        results.append(analyse_graph(name, A_cr, 6 if p > 3 else 3, desc))
        # Note: CR(p) as Cayley on Z_p has degree |S ∪ -S|.
        # S = {1, g, g²} and -S = {p-1, p-g, p-g²}.
        # If none of s ≡ -s (mod p), |connection set| = 6, so 6-regular.
        # (For cubic arc-transitive graphs the VOLTAGE LIFT of K₂ gives 3-regular.)

    # -----------------------------------------------------------------------
    # 3. K₂-cover (voltage lift) versions — these are 3-regular and bipartite
    # -----------------------------------------------------------------------

    print("\n" + "#"*65)
    print("  SECTION 3: CR(p) as Z_p-covers of K₂ (3-regular, bipartite)")
    print("#"*65)

    for p in cr_primes:
        A_lift, voltages = build_CR_via_voltage(p)
        name = f"CR({p}) voltage lift of K₂"
        desc = f"2p={2*p} vertices, voltages {voltages}"
        results.append(analyse_graph(name, A_lift, 3, desc))

    # -----------------------------------------------------------------------
    # 4. Analytical eigenvalue check for K₂-covers
    # -----------------------------------------------------------------------

    print("\n" + "#"*65)
    print("  SECTION 4: Analytical eigenvalues for Z_p-covers of K₂")
    print("#"*65)

    print("\n  Analytical formula: eigenvalues = ±|Σ_{v ∈ V} ω^(j·v)|")
    print("  for j = 0, 1, ..., p-1, where ω = e^(2πi/p).\n")

    for p in cr_primes:
        S = cubic_residues_mod(p)
        ana_eigs = k2_cover_eigenvalues_analytical(p, S)
        A_lift, _ = build_CR_via_voltage(p)
        num_eigs = compute_spectrum(A_lift)
        match = np.allclose(np.sort(ana_eigs)[::-1], num_eigs, atol=1e-8)
        print(f"  CR({p}) lift: analytical vs numerical match = {'PASS ✓' if match else 'FAIL ✗'}")
        # Show the three distinct |f(j)| values
        ana_pos = np.sort(ana_eigs[ana_eigs > 1e-8])[::-1]
        print(f"    Largest analytical eigenvalues: {ana_pos[:5].round(6)}")

    # -----------------------------------------------------------------------
    # 5. Summary table
    # -----------------------------------------------------------------------

    print("\n\n" + "="*80)
    print("  SUMMARY TABLE — All Computed Graphs")
    print("="*80)
    header = f"  {'Graph':<40} {'|V|':>5} {'d':>3} {'λ₂':>9} {'gap':>7} {'Girth':>5} {'Ram?':>5}"
    print(header)
    print("  " + "─"*78)
    for r in results:
        ram_str = "YES" if r['is_ramanujan'] else "NO"
        print(f"  {r['name']:<40} {r['n']:>5} {r['degree']:>3} "
              f"{r['lambda2']:>9.4f} {r['spectral_gap']:>7.4f} "
              f"{r['girth']:>5} {ram_str:>5}")

    # -----------------------------------------------------------------------
    # 6. CR(p) spectral evolution table
    # -----------------------------------------------------------------------

    print("\n\n" + "="*80)
    print("  CR(p) VOLTAGE LIFT SPECTRAL EVOLUTION  (3-regular, bipartite)")
    print("="*80)
    print(f"  {'p':>4}  {'|V|':>5}  {'λ₂':>9}  {'λ_min':>9}  "
          f"{'gap':>7}  {'Fiedler':>8}  {'Girth':>5}  {'Ram?':>5}")
    print("  " + "─"*65)

    for r in results:
        if "voltage lift" in r['name']:
            p = int(r['name'].split("(")[1].split(")")[0])
            ram_str = "YES" if r['is_ramanujan'] else "NO"
            print(f"  {p:>4}  {r['n']:>5}  {r['lambda2']:>9.4f}  {r['lambda_min']:>9.4f}  "
                  f"{r['spectral_gap']:>7.4f}  {r['fiedler']:>8.4f}  "
                  f"{r['girth']:>5}  {ram_str:>5}")

    # -----------------------------------------------------------------------
    # 7. Weil bound check
    # -----------------------------------------------------------------------

    print("\n\n" + "="*80)
    print("  WEIL BOUND CHECK for CR(p) voltage lifts")
    print("  |λ| ≤ 2√p  (Weil's theorem on character sums over F_p)")
    print("="*80)
    print(f"  {'p':>4}  {'2√p':>8}  {'max|λ_nontrivial|':>18}  {'Weil satisfied?':>16}")
    print("  " + "─"*50)

    for p in cr_primes:
        A_lift, S = build_CR_via_voltage(p)
        eigs = compute_spectrum(A_lift)
        weil_bound = 2.0 * np.sqrt(p)
        # non-trivial = all except ±3 (the largest eigenvalues)
        nontrivial = [e for e in eigs if abs(abs(e) - 3.0) > 1e-6]
        max_nontrivial = max(abs(e) for e in nontrivial) if nontrivial else 0.0
        satisfied = max_nontrivial <= weil_bound + 1e-8
        print(f"  {p:>4}  {weil_bound:>8.4f}  {max_nontrivial:>18.6f}  "
              f"{'YES' if satisfied else 'NO':>16}")
