/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# From Cellular Surfaces to CSS Quantum Codes

Assembly layer: given a `CellularSurface` with the appropriate hypotheses
(connectivity, Euler's formula, two-sides condition), construct a `CSSFromTiling`
and prove that it encodes exactly 2g logical qubits.

This connects the chain complex ∂² = 0 (Phase 3) and the rank theorems (Phase 4)
with the existing `css_k_eq_2g` theorem.
-/

import CSSCode.CSSCode
import CSSCode.IncMatrixRank

open Matrix CellularSurface

/-! ### The surface tiling from a CellularSurface -/

/-- Extract a `SurfaceTiling` from a `CellularSurface` with given genus. -/
def CellularSurface.toSurfaceTiling (S : CellularSurface) (g : ℕ)
    (hV : 0 < S.nV) (hF : 0 < S.nF)
    (euler : (S.nV : ℤ) - S.nE + S.nF = 2 - 2 * g) : SurfaceTiling where
  V := S.nV
  E := S.nE
  F := S.nF
  g := g
  hV := hV
  hF := hF
  euler := euler

/-! ### The CSS code from a CellularSurface -/

/-- Construct a `CSSFromTiling` from a `CellularSurface`.

    - `H_X = d1` (vertex-edge incidence, the X-stabilizers)
    - `H_Z = d2ᵀ` (face-edge incidence transposed, the Z-stabilizers)
    - Orthogonality `H_X * H_Zᵀ = 0` follows from `d1 * d2 = 0` (Phase 3) -/
def CellularSurface.toCSSFromTiling (S : CellularSurface) (g : ℕ)
    (hV : 0 < S.nV) (hF : 0 < S.nF)
    (euler : (S.nV : ℤ) - S.nE + S.nF = 2 - 2 * g) :
    CSSFromTiling (S.toSurfaceTiling g hV hF euler) where
  H_X := S.d1
  H_Z := S.d2ᵀ
  orthogonal := by
    -- H_X * H_Zᵀ = d1 * (d2ᵀ)ᵀ = d1 * d2 = 0
    rw [Matrix.transpose_transpose]
    exact S.d1_mul_d2_eq_zero

/-! ### The main theorem: k = 2g from a CellularSurface -/

/-- **The CSS surface code theorem from combinatorial data.**

    Given a `CellularSurface` satisfying:
    - primal graph is connected
    - dual graph is connected
    - each edge borders exactly 2 faces (two-sides condition)
    - Euler's formula holds with genus g

    the resulting CSS quantum code encodes exactly `2g` logical qubits. -/
theorem CellularSurface.css_k_eq_2g_from_surface (S : CellularSurface) (g : ℕ)
    (hV : 1 < S.nV) (hF : 1 < S.nF)
    (euler : (S.nV : ℤ) - S.nE + S.nF = 2 - 2 * g)
    (hconn : S.toSimpleGraph.Connected)
    (hconn_dual : S.toDualSimpleGraph.Connected)
    (htwo : ∀ e, (Finset.univ.filter fun f : Fin S.nF => S.d2 e f ≠ 0).card = 2) :
    (S.toCSSFromTiling g (by omega) (by omega) euler).k = 2 * g := by
  apply css_k_eq_2g
  · exact S.d1_rank_eq hconn hV
  · simp only [toCSSFromTiling, rank_transpose]
    exact S.d2_rank_eq hconn_dual hF htwo
