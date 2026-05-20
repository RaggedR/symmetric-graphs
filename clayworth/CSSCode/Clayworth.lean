/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Clayworth CSS Code

The Clayworth graphs are a dual pair of tilings of a genus-505 surface,
constructed from the exceptional Lie group G₂(2) (order 12096) via the
Sabidussi coset construction:

  Clayworth-3:  Sab(G₂(2), C₃)  — {12,3} tiling — 4032 vertices, degree 3, girth 12
  Clayworth-12: Sab(G₂(2), C₁₂) — {3,12} tiling — 1008 vertices, degree 12, girth 4

Together they give a CSS quantum code [[6048, 1010, ≥ 4]].

This file defines the tiling parameters, verifies Euler's formula by computation,
and proves that the CSS code encodes exactly 1010 = 2 × 505 logical qubits
(conditional on the rank hypotheses, which were verified numerically in Python).
-/

import CSSCode.CSSCode

/-- The Clayworth surface tiling of genus 505.

  - V = 4032 vertices (degree 3, primal {12,3} tiling)
  - E = 6048 edges
  - F = 1008 dodecagonal faces (= vertices of dual {3,12} tiling)
  - g = 505

  Euler check: 4032 - 6048 + 1008 = -1008 = 2 - 2 × 505. -/
def clayworthTiling : SurfaceTiling where
  V := 4032
  E := 6048
  F := 1008
  g := 505
  hV := by norm_num
  hF := by norm_num
  euler := by norm_num

/-- The CSS code from the Clayworth tiling encodes exactly 1010 logical qubits.

  This follows from `css_k_eq_2g` applied to the Clayworth tiling parameters.
  The rank hypotheses (rank(H_X) = 4031, rank(H_Z) = 1007) were verified
  by Gaussian elimination over F₂ in Python (see `clayworth_codes.py`). -/
theorem clayworth_k_eq_1010 (css : CSSFromTiling clayworthTiling)
    (hX : (css.H_X.rank : ℤ) = 4031)
    (hZ : (css.H_Z.rank : ℤ) = 1007) :
    css.k = 1010 :=
  css_k_eq_2g clayworthTiling css hX hZ

/-! ## Additional Clayworth properties

The following are stated as definitions / lemmas for documentation.
The numerical values were computed in GAP and Python. -/

/-- The Clayworth code has 6048 physical qubits (one per edge). -/
theorem clayworth_nQubits (css : CSSFromTiling clayworthTiling) :
    css.nQubits = 6048 := rfl

/-- Girth of the primal graph (lower bound on X-distance).
    Computed by GAP: the shortest cycle in Clayworth-3 has length 12. -/
def clayworth_girth_primal : ℕ := 12

/-- Girth of the dual graph (lower bound on Z-distance).
    Computed by GAP: the shortest cycle in Clayworth-12 has length 4. -/
def clayworth_girth_dual : ℕ := 4

/-! ## The Klein quartic: genus 3, PSL(2,7)

The Klein quartic is the Hurwitz surface of smallest genus: PSL(2,7) achieves
the 84(g-1) bound on automorphisms. The triangle group Δ(2,3,7) gives a {7,3}
tiling by heptagons with dual {3,7} tiling by triangles. -/

/-- The Klein quartic tiling of genus 3.
  V = 56, E = 84, F = 24. χ = 56 - 84 + 24 = -4 = 2 - 2×3. -/
def kleinQuarticTiling : SurfaceTiling where
  V := 56
  E := 84
  F := 24
  g := 3
  hV := by norm_num
  hF := by norm_num
  euler := by norm_num

/-- The Klein quartic CSS code encodes 6 logical qubits. -/
theorem klein_k_eq_6 (css : CSSFromTiling kleinQuarticTiling)
    (hX : (css.H_X.rank : ℤ) = 55)
    (hZ : (css.H_Z.rank : ℤ) = 23) :
    css.k = 6 :=
  css_k_eq_2g kleinQuarticTiling css hX hZ

-- The Heawood / Fano plane tiling for comparison: genus 1 (torus)
/-- The Heawood graph tiles a torus (genus 1).
    V = 14, E = 21, F = 7 heptagonal faces.
    χ = 14 - 21 + 7 = 0 = 2 - 2×1. -/
def heawoodTiling : SurfaceTiling where
  V := 14
  E := 21
  F := 7
  g := 1
  hV := by norm_num
  hF := by norm_num
  euler := by norm_num

/-- The Heawood CSS code encodes 2 logical qubits.
    The Heawood graph is the incidence graph of PG(2,2) = Fano plane,
    and its CSS code on the torus is the simplest non-trivial surface code. -/
theorem heawood_k_eq_2 (css : CSSFromTiling heawoodTiling)
    (hX : (css.H_X.rank : ℤ) = 13)
    (hZ : (css.H_Z.rank : ℤ) = 6) :
    css.k = 2 :=
  css_k_eq_2g heawoodTiling css hX hZ
