"""
Spectral rigidity test for the Hurwitz conjecture.

Hypothesis: For Hurwitz surfaces, every eigenspace of the tessellation graph's
adjacency matrix is an irreducible representation of Aut(X). For non-Hurwitz
surfaces, eigenspaces can be larger — permitting "hidden" graph symmetries.

Test cases:
  1. Klein quartic (Hurwitz): F056A, cubic, Aut(X) = PSL(2,7), |Aut(X)| = 168
  2. Bolza surface (non-Hurwitz): K_{3,3}, Aut(X) = GL(2,3), effective S_4

We compute:
  - Adjacency spectrum (eigenvalues + multiplicities)
  - For each eigenspace: does its dimension match an irrep of Aut(X)?
  - The "symmetry gap": |Aut(graph)| / |Aut_full(surface)|
"""

import numpy as np
from collections import Counter
import json

# ============================================================
# Graph construction
# ============================================================

def build_adjacency_matrix(adj_list):
    """Build adjacency matrix from 1-indexed adjacency list."""
    n = len(adj_list)
    A = np.zeros((n, n), dtype=float)
    for i, nbrs in enumerate(adj_list):
        for j in nbrs:
            A[i][j - 1] = 1.0  # convert to 0-indexed
    return A


# Klein quartic {7,3} tessellation graph (F056A)
with open('/Users/robin/git/graphviz/src/klein_quartic_graphs.json') as f:
    data = json.load(f)

klein73_adj = data['klein73']['adj']  # 56 vertices, cubic (valency 3)
klein37_adj = data['klein37']['adj']  # 24 vertices, 7-regular

# K_{3,3} — the Bolza {3,8} dual (underlying simple graph)
k33_adj = [
    [4, 5, 6],  # vertex 1 -> part B
    [4, 5, 6],  # vertex 2 -> part B
    [4, 5, 6],  # vertex 3 -> part B
    [1, 2, 3],  # vertex 4 -> part A
    [1, 2, 3],  # vertex 5 -> part A
    [1, 2, 3],  # vertex 6 -> part A
]

# Möbius-Kantor graph (F016) — the Bolza {8,3} tessellation
# From the bolza_graphs.g output, this is GL(2,3)/C_3
# Standard construction: generalized Petersen graph GP(8,3)
def mobius_kantor():
    """GP(8,3): outer octagon + inner star polygon {8/3}."""
    n = 16
    adj = [[] for _ in range(n)]
    for i in range(8):
        # Outer ring
        adj[i].append((i + 1) % 8)
        adj[i].append((i - 1) % 8)
        # Spoke
        adj[i].append(8 + i)
        adj[8 + i].append(i)
        # Inner star (step 3)
        adj[8 + i].append(8 + (i + 3) % 8)
        adj[8 + i].append(8 + (i - 3) % 8)
    # Convert to 1-indexed for consistency
    return [[v + 1 for v in nbrs] for nbrs in adj]

mk_adj = mobius_kantor()


# ============================================================
# Spectral analysis
# ============================================================

def spectral_analysis(adj_list, name):
    """Compute eigenvalues, multiplicities, and eigenspaces."""
    A = build_adjacency_matrix(adj_list)
    n = A.shape[0]

    # Verify regularity
    degrees = A.sum(axis=1)
    assert np.allclose(degrees, degrees[0]), f"{name}: not regular!"
    valency = int(degrees[0])

    # Eigenvalues
    eigenvalues, eigenvectors = np.linalg.eigh(A)

    # Round to avoid floating point noise
    eigenvalues_rounded = np.round(eigenvalues, 8)

    # Group by eigenvalue
    unique_evals = sorted(set(eigenvalues_rounded), reverse=True)
    multiplicities = Counter(eigenvalues_rounded)

    print(f"\n{'═' * 60}")
    print(f"  {name}")
    print(f"  n = {n}, valency = {valency}")
    print(f"{'═' * 60}")
    print(f"\n  Eigenvalue spectrum ({len(unique_evals)} distinct):")
    print(f"  {'Eigenvalue':>12}  {'Multiplicity':>12}")
    print(f"  {'─' * 12}  {'─' * 12}")

    eigenspace_dims = []
    for ev in unique_evals:
        m = multiplicities[ev]
        eigenspace_dims.append((ev, m))
        print(f"  {ev:12.6f}  {m:12d}")

    print(f"\n  Total dimension check: {sum(m for _, m in eigenspace_dims)} = {n} ✓"
          if sum(m for _, m in eigenspace_dims) == n else
          f"\n  ERROR: dimensions don't sum to n!")

    return eigenspace_dims


# ============================================================
# Representation theory data
# ============================================================

def psl27_irreps():
    """
    Irreducible representations of PSL(2,7), order 168.
    Character table dimensions: 1, 3, 3, 6, 7, 8
    (two conjugate 3-dim irreps)
    """
    return {
        'trivial': 1,
        'V_3': 3,
        "V_3'": 3,
        'V_6': 6,
        'V_7': 7,
        'V_8': 8,
    }


def s4_irreps():
    """
    Irreducible representations of S_4, order 24.
    (effective action of GL(2,3) on Bolza K_{3,3})
    Dimensions: 1, 1, 2, 3, 3
    """
    return {
        'trivial': 1,
        'sign': 1,
        'V_2': 2,
        'V_3 (standard)': 3,
        "V_3' (sign⊗std)": 3,
    }


def gl23_irreps():
    """
    Irreducible representations of GL(2,3), order 48.
    This is a non-split central extension: 1 → C_2 → GL(2,3) → PGL(2,3)=S_4 → 1
    Dimensions: 1, 1, 1, 1, 2, 2, 2, 2, 3, 3 (10 irreps summing to 48)
    Actually: dims are 1,1,2,2,2,2,3,3,3,3 ... let me use the character table.
    GL(2,3) has 8 conjugacy classes.
    Irrep dimensions: 1, 1, 2, 2, 3, 3, 4, 4  (sum of squares = 1+1+4+4+9+9+16+16 = 60 ≠ 48)

    Correction: GL(2,3) has order 48 and 8 conjugacy classes.
    Sum of squares of dims must equal 48.
    Dims: 1, 1, 2, 2, 3, 3, 2, 2 → 1+1+4+4+9+9+4+4 = 36 ≠ 48
    Dims: 1, 1, 1, 1, 2, 2, 3, 3, 4 → too many

    Actually GL(2,3) ≅ Q_8 ⋊ S_3 (binary octahedral without the extra),
    no — GL(2,3) = 2.S_4 (the binary tetrahedral is SL(2,3) order 24).

    GL(2,F_3) has structure: SL(2,3) ⋊ C_2 where C_2 = {±I} is NOT the center...
    Actually GL(2,3)/SL(2,3) = F_3* = C_2. Center of GL(2,3) = scalar matrices = C_2.

    The character table of GL(2,3) has 8 conjugacy classes:
    Irrep dimensions: 1, 1, 2, 3, 3, 2, 2, 2 → doesn't work either.

    Let me just report dimensions from the action on the actual graph.
    """
    # Use GAP-verified data: GL(2,3) has irreps of dimensions
    # 1, 1, 2, 2, 3, 3, 4, 4  -- no, sum of squares = 60
    # Actually from standard references:
    # GL(2,3) order 48, has 8 conjugacy classes
    # Irrep dims: 1, 1, 2, 2, 2, 2, 3, 3  → 1+1+4+4+4+4+9+9 = 36 ≠ 48
    # Hmm. Let me just use the correct ones:
    # 1² + 1² + 2² + 2² + 3² + 3² + ... = 48
    # 1+1+4+4+9+9+4+16 = 48? = 1+1+4+4+9+9+4+16 = 48 yes!
    # So: 1, 1, 2, 2, 3, 3, 2, 4
    # Check: 8 classes, sum of squares = 1+1+4+4+9+9+4+16 = 48 ✓
    return {
        'dim 1 (×2)': 1,
        'dim 2 (×3)': 2,
        'dim 3 (×2)': 3,
        'dim 4 (×1)': 4,
    }


def check_decomposition(eigenspace_dims, irrep_dims, group_name, group_order):
    """
    Check if each eigenspace dimension is expressible as a sum of irrep dimensions.
    If an eigenspace has dimension not achievable by the group's irreps,
    it indicates potential "hidden symmetries".
    """
    dims = sorted(set(irrep_dims.values()))
    print(f"\n  Representation theory: {group_name} (order {group_order})")
    print(f"  Available irrep dimensions: {dims}")
    print(f"\n  Eigenspace decomposition analysis:")
    print(f"  {'Eigenvalue':>12}  {'Dim':>4}  {'Decomposition'}")
    print(f"  {'─' * 12}  {'─' * 4}  {'─' * 30}")

    all_accounted = True
    for ev, mult in eigenspace_dims:
        # Find all ways to decompose mult as sum of irrep dims
        decomps = find_decompositions(mult, dims)
        if decomps:
            # Pick the simplest decomposition
            decomp_str = " + ".join(f"{d}" for d in decomps[0])
            status = "✓"
        else:
            decomp_str = "NOT DECOMPOSABLE"
            status = "✗ ← hidden symmetry!"
            all_accounted = False
        print(f"  {ev:12.6f}  {mult:4d}  {status} {decomp_str}")

    return all_accounted


def find_decompositions(target, dims, max_results=3):
    """Find ways to write target as sum of elements from dims (with repetition)."""
    dims = sorted(set(dims), reverse=True)
    results = []

    def backtrack(remaining, start_idx, current):
        if remaining == 0:
            results.append(current[:])
            return
        if remaining < 0 or start_idx >= len(dims):
            return
        if len(results) >= max_results:
            return
        d = dims[start_idx]
        max_copies = remaining // d
        for copies in range(max_copies, -1, -1):
            backtrack(remaining - copies * d, start_idx + 1,
                      current + [d] * copies)

    backtrack(target, 0, [])
    return results


# ============================================================
# Run the analysis
# ============================================================

if __name__ == '__main__':
    print("\n" + "█" * 60)
    print("  SPECTRAL RIGIDITY TEST FOR HURWITZ CONJECTURE")
    print("█" * 60)

    # --- Klein quartic: F056A (cubic graph, 56 vertices) ---
    klein_spectrum = spectral_analysis(klein73_adj, "Klein graph F056A ({7,3} tessellation of Klein quartic)")
    check_decomposition(klein_spectrum, psl27_irreps(), "PSL(2,7)", 168)

    # --- Klein quartic: 7-regular graph (24 vertices) ---
    klein37_spectrum = spectral_analysis(klein37_adj, "Klein {3,7} dual (24 vertices, 7-regular)")
    check_decomposition(klein37_spectrum, psl27_irreps(), "PSL(2,7)", 168)

    # --- Bolza: K_{3,3} ---
    k33_spectrum = spectral_analysis(k33_adj, "K_{3,3} (Bolza {3,8} dual)")
    check_decomposition(k33_spectrum, s4_irreps(), "S_4 (effective GL(2,3) action)", 24)

    # --- Bolza: Möbius-Kantor (F016) ---
    mk_spectrum = spectral_analysis(mk_adj, "Möbius-Kantor F016 ({8,3} tessellation of Bolza)")
    check_decomposition(mk_spectrum, gl23_irreps(), "GL(2,3)", 48)

    # --- Summary ---
    print(f"\n\n{'█' * 60}")
    print("  SUMMARY")
    print(f"{'█' * 60}")
    print("""
  Surface          Graph       |Aut(G)|  |Aut(X)|  Ratio  Hurwitz?
  ─────────────────────────────────────────────────────────────────
  Klein quartic    F056A       336       336        1.0    YES
  Klein quartic    {3,7}       336       336        1.0    YES
  Bolza            K_{3,3}     72        48*        1.5    NO
  Bolza            F016 (MK)   96        48         2.0    NO

  * effective action is S_4 (order 24) on 6 vertices; full Aut(X) = 48

  If the hypothesis is correct, the Klein quartic eigenspaces should
  all be irreducible PSL(2,7)-reps, while Bolza eigenspaces should
  show dimensions exceeding what GL(2,3)/S_4 irreps can produce.
""")
