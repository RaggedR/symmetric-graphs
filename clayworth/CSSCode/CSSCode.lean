/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# CSS Quantum Codes from Surface Tilings

A tiling of a closed orientable surface of genus g by polygons gives a
CSS (Calderbank-Shor-Steane) quantum error-correcting code. The X-checks
come from the vertex-edge incidence matrix of the primal tiling, and the
Z-checks come from the face-edge incidence matrix (equivalently, the
vertex-edge incidence of the dual tiling).

The main theorem: the number of logical qubits equals 2g, the first
Betti number of the surface. This follows purely from:
  - rank(H_X) = V - 1  (connected primal graph)
  - rank(H_Z) = F - 1  (connected dual graph)
  - Euler's formula: V - E + F = 2 - 2g
-/

import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.Matrix.Rank

open Matrix

/-- A tiling of a closed orientable surface of genus `g`.
    The Euler characteristic is V - E + F = 2 - 2g.
    We use ℤ throughout to avoid natural number subtraction issues. -/
structure SurfaceTiling where
  /-- Number of vertices in the primal tiling -/
  V : ℕ
  /-- Number of edges (shared between primal and dual) -/
  E : ℕ
  /-- Number of faces in the primal tiling (= vertices of dual) -/
  F : ℕ
  /-- Genus of the surface -/
  g : ℕ
  /-- The primal graph is non-degenerate -/
  hV : 0 < V
  /-- The dual graph is non-degenerate -/
  hF : 0 < F
  /-- Euler's formula for a closed orientable surface -/
  euler : (V : ℤ) - E + F = 2 - 2 * g

/-- A CSS quantum code constructed from a surface tiling.
    - `H_X` is the vertex-edge incidence matrix (X-stabilizers)
    - `H_Z` is the face-edge incidence matrix (Z-stabilizers)
    - The CSS condition requires H_X · H_Zᵀ = 0 (every face boundary is a cycle) -/
structure CSSFromTiling (T : SurfaceTiling) where
  /-- X-check matrix: vertex-edge incidence of primal (V × E over F₂) -/
  H_X : Matrix (Fin T.V) (Fin T.E) (ZMod 2)
  /-- Z-check matrix: face-edge incidence / vertex-edge incidence of dual (F × E over F₂) -/
  H_Z : Matrix (Fin T.F) (Fin T.E) (ZMod 2)
  /-- CSS orthogonality: the image of H_Zᵀ lies in the kernel of H_X.
      Equivalently, every face boundary (row of H_Z) is a cycle (in ker H_X).
      This encodes the topological fact that ∂² = 0 in the chain complex. -/
  orthogonal : H_X * H_Zᵀ = 0

namespace CSSFromTiling

variable {T : SurfaceTiling}

/-- The number of physical qubits = number of edges -/
def nQubits (_ : CSSFromTiling T) : ℕ := T.E

/-- The number of logical qubits as an integer: k = E - rank(H_X) - rank(H_Z).
    We use ℤ to make the arithmetic with Euler's formula clean.
    This equals dim(H₁(surface; F₂)) = 2g ≥ 0. -/
noncomputable def k (css : CSSFromTiling T) : ℤ :=
  (T.E : ℤ) - css.H_X.rank - css.H_Z.rank

end CSSFromTiling

/-! ## The Main Theorem

For a surface tiling where both the primal and dual graphs are connected,
the parity check matrices have rank V-1 and F-1 respectively. Then:

  k = E - (V - 1) - (F - 1) = E - V - F + 2 = 2g

This is purely arithmetic given the rank hypotheses and Euler's formula.
-/

/-- **The CSS surface code theorem**: the number of logical qubits encoded by a
    CSS code from a genus-g surface tiling equals 2g (the first Betti number).

    The hypotheses `hX` and `hZ` state that both the primal and dual graphs are
    connected (rank of incidence matrix = number of vertices - 1). These would
    follow from `SimpleGraph.Connected` and a theorem about incidence matrix rank
    over F₂, which is not yet in Mathlib. -/
theorem css_k_eq_2g (T : SurfaceTiling) (css : CSSFromTiling T)
    (hX : (css.H_X.rank : ℤ) = T.V - 1)
    (hZ : (css.H_Z.rank : ℤ) = T.F - 1) :
    css.k = 2 * T.g := by
  unfold CSSFromTiling.k
  have := T.euler
  linarith
