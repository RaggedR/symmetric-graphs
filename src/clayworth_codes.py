#!/usr/bin/env python3
"""
Compute error-correcting codes from the Clayworth graphs.

Clayworth-3: {12,3} tiling of genus-505 surface, Sab(G_2(2), C_3)
  4032 vertices, 6048 edges, cubic, girth 12

Clayworth-12: {3,12} dual tiling, Sab(G_2(2), C_12)
  1008 vertices, 6048 edges, 12-regular, girth 4

Together they define a CSS quantum code via the surface tiling/dual construction.
"""

import json
import sys
import os

# Import F_2 linear algebra from graph_codes.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graph_codes import (
    incidence_matrix_f2, rref_f2, null_space_f2, row_space_f2,
    is_bipartite, biadjacency_matrix_f2, compute_girth, edges_from_adj
)


def rank_f2(matrix):
    """Compute rank of a matrix over F_2."""
    if not matrix:
        return 0
    _, _, r = rref_f2(matrix)
    return r


def boundary_matrices(adj_primal, n_primal, adj_dual, n_dual, n_edges):
    """
    Construct the boundary operators for the CSS code from a surface tiling.

    A tiling of a genus-g surface gives a chain complex:
        C_2 --d2--> C_1 --d1--> C_0

    where:
      C_0 = F_2^V (vertices of primal = faces of dual)
      C_1 = F_2^E (edges, shared between primal and dual)
      C_2 = F_2^F (faces of primal = vertices of dual)

    d1 = vertex-edge incidence matrix of primal (V x E)
    d2 = face-edge incidence matrix = vertex-edge incidence of dual, transposed (E x F)

    The CSS code has:
      X-checks = rows of d1 (vertex stabilizers)
      Z-checks = columns of d2 = rows of incidence matrix of dual

    Code parameters:
      n = E (number of qubits = edges)
      k = dim H_1(surface; F_2) = 2g (first Betti number for orientable surface)
      d_X >= girth of primal
      d_Z >= girth of dual
      d = min(d_X, d_Z)
    """
    # d1: vertex-edge incidence of primal (n_primal x n_edges)
    B_primal, edges_primal = incidence_matrix_f2(adj_primal, n_primal)

    # d2^T: vertex-edge incidence of dual (n_dual x n_edges)
    # But the dual's edges are the SAME edges, just with different endpoints
    # For the CSS code we need: the face-edge incidence matrix
    # Each face of the primal corresponds to a vertex of the dual
    # A face is bounded by the edges incident to the corresponding dual vertex
    B_dual, edges_dual = incidence_matrix_f2(adj_dual, n_dual)

    return B_primal, edges_primal, B_dual, edges_dual


def css_code_params(adj3, n3, adj12, n12):
    """
    Compute CSS quantum code parameters from the Clayworth tiling pair.

    The genus-505 surface has first Betti number beta_1 = 2 * 505 = 1010.
    This should equal n - rank(H_X) - rank(H_Z) where:
      H_X = vertex-edge incidence of Clayworth-3 (primal)
      H_Z = vertex-edge incidence of Clayworth-12 (dual)
    """
    print("=" * 60)
    print("  CSS Quantum Code from Clayworth Tiling")
    print("=" * 60)
    print()

    # Primal: Clayworth-3 ({12,3} tiling)
    # Vertices = V, Faces = F (dodecagons)
    V = n3      # 4032
    F = n12     # 1008 (faces = vertices of dual)
    E_expected = 6048

    print(f"  Primal tiling: {{{12},{3}}}")
    print(f"    V = {V} (vertices, degree 3)")
    print(f"    F = {F} (dodecagonal faces)")
    print(f"    E = {E_expected} (edges)")
    print(f"    chi = V - E + F = {V - E_expected + F}")
    print(f"    genus = {(2 - (V - E_expected + F)) // 2}")
    print()

    # Compute H_X: vertex-edge incidence of primal (V x E matrix)
    print("  Computing H_X (primal vertex-edge incidence)...")
    B_primal, edges_p = incidence_matrix_f2(adj3, n3)
    E = len(edges_p)
    print(f"    Matrix size: {V} x {E}")

    rank_X = rank_f2(B_primal)
    print(f"    rank(H_X) = {rank_X}")
    print(f"    Expected: V - 1 = {V - 1} (connected graph)")
    print()

    # Compute H_Z: face-edge incidence = vertex-edge incidence of dual (F x E matrix)
    # BUT: the dual graph's edges must be indexed the same way as the primal's edges.
    #
    # The dual's vertices are the primal's faces. Two dual vertices are adjacent
    # iff the corresponding primal faces share an edge. Each primal edge separates
    # exactly two faces (on a closed surface).
    #
    # Since both graphs tile the SAME surface with the SAME 6048 edges,
    # we can compute rank(H_Z) from the dual's vertex-edge incidence.
    # But we need the edges indexed consistently.
    #
    # For the CSS code, the key identity is:
    #   k = E - rank(H_X) - rank(H_Z) = 2g (first Betti number)
    #
    # We can compute rank(H_Z) = E - rank(H_X) - 2g = E - (V-1) - 2g
    # Or compute it directly from the dual.

    print("  Computing H_Z (dual vertex-edge incidence)...")
    B_dual, edges_d = incidence_matrix_f2(adj12, n12)
    E_dual = len(edges_d)
    print(f"    Matrix size: {F} x {E_dual}")

    rank_Z = rank_f2(B_dual)
    print(f"    rank(H_Z) = {rank_Z}")
    print(f"    Expected: F - 1 = {F - 1} (connected graph)")
    print()

    # CSS code parameters
    # n (qubits) = E
    # k (logical qubits) = E - rank(H_X) - rank(H_Z)
    # For a closed orientable surface of genus g: k = 2g
    n_qubits = E
    k_logical = E - rank_X - rank_Z
    genus = (2 - (V - E + F)) // 2

    print("  --- CSS Code Parameters ---")
    print(f"    n (physical qubits) = {n_qubits}")
    print(f"    k (logical qubits)  = {k_logical}")
    print(f"    Expected k = 2g = {2 * genus}")
    print()

    # Distance bounds
    # d_X >= girth of primal (Clayworth-3, girth 12)
    # d_Z >= girth of dual (Clayworth-12, girth 4 from GAP output)
    girth_primal = 12  # computed by GAP
    girth_dual = 4     # computed by GAP
    d_lower = min(girth_primal, girth_dual)

    print(f"    d_X >= girth(primal) = {girth_primal}")
    print(f"    d_Z >= girth(dual)   = {girth_dual}")
    print(f"    d >= min(d_X, d_Z)   >= {d_lower}")
    print()
    print(f"  CSS Code: [[{n_qubits}, {k_logical}, >= {d_lower}]]")
    print()

    # Classical codes from each graph individually
    print("  --- Classical Cycle Codes ---")
    print()

    # Cycle code of Clayworth-3
    k_cycle3 = E - rank_X
    print(f"  Clayworth-3 cycle code:")
    print(f"    [{E}, {k_cycle3}, {girth_primal}]")
    print(f"    Rate: {k_cycle3/E:.4f}")
    print()

    # Cycle code of Clayworth-12
    k_cycle12 = E_dual - rank_Z
    print(f"  Clayworth-12 cycle code:")
    print(f"    [{E_dual}, {k_cycle12}, {girth_dual}]")
    print(f"    Rate: {k_cycle12/E_dual:.4f}")
    print()

    # Check if bipartite
    bip3, _, _ = is_bipartite(adj3, n3)
    bip12, _, _ = is_bipartite(adj12, n12)
    print(f"  Clayworth-3 bipartite: {bip3}")
    print(f"  Clayworth-12 bipartite: {bip12}")

    return {
        'n_qubits': n_qubits,
        'k_logical': k_logical,
        'd_lower': d_lower,
        'genus': genus,
        'girth_primal': girth_primal,
        'girth_dual': girth_dual,
        'rank_X': rank_X,
        'rank_Z': rank_Z,
    }


def main():
    # Load Clayworth graphs
    json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'clayworth_graphs.json')
    if not os.path.exists(json_path):
        print(f"ERROR: {json_path} not found. Run clayworth_build.g first.")
        sys.exit(1)

    with open(json_path) as f:
        data = json.load(f)

    cw3 = data['clayworth3']
    cw12 = data['clayworth12']

    print(f"Loaded Clayworth-3: {cw3['n']} vertices")
    print(f"Loaded Clayworth-12: {cw12['n']} vertices")
    print()

    # Compute CSS code
    result = css_code_params(cw3['adj'], cw3['n'], cw12['adj'], cw12['n'])

    print()
    print("=" * 60)
    print("  Summary")
    print("=" * 60)
    nq = result['n_qubits']
    rx, rz = result['rank_X'], result['rank_Z']
    gp, gd = result['girth_primal'], result['girth_dual']
    print(f"""
  The Clayworth graphs tile a genus-{result['genus']} surface with
  {nq} edges. The CSS quantum code encodes
  {result['k_logical']} logical qubits into {nq} physical qubits.

  This is an LDPC code: X-stabilizers have weight 3 (from the
  cubic primal graph) and Z-stabilizers have weight 12 (from
  the 12-regular dual). The asymmetry {gp} vs {gd} in
  the girths means X-errors are much harder to create than
  Z-errors.

  Classical codes:
    Clayworth-3 cycle code:  [{nq}, {nq - rx}, {gp}]
    Clayworth-12 cycle code: [{nq}, {nq - rz}, {gd}]
""")


if __name__ == "__main__":
    main()
