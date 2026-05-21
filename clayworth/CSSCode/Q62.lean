/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Parabolic Quadric Q(6,2) and the G₂(2) CSS Code

The 63 points of the parabolic quadric Q(6,2) in PG(6,2) are the nonzero
vectors x ∈ 𝔽₂⁷ satisfying Q(x) = x₀x₄ + x₁x₅ + x₂x₆ + x₃² = 0.

These 63 points are the points of the split Cayley hexagon GH(2,2),
whose automorphism group is the exceptional Chevalley group G₂(2) of
order 12096. The incidence graph of GH(2,2) is the Tutte 12-cage.

Using these points as columns of a 7 × 63 parity-check matrix H,
we verify:
  1. Q(6,2) has exactly 63 points
  2. H · Hᵀ = 0 over 𝔽₂ (the self-orthogonality condition)

This gives a self-orthogonal [63, 56, 3] classical code, and hence
a CSS quantum code [[63, 49, 3]] encoding 49 logical qubits.
-/

import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Fin.VecNotation

open Matrix

/-! ## The quadratic form Q on 𝔽₂⁷ -/

/-- The quadratic form Q(x) = x₀x₄ + x₁x₅ + x₂x₆ + x₃² on 𝔽₂⁷.
    Points of Q(6,2) are the nonzero vectors with Q(x) = 0. -/
def Q62form (x : Fin 7 → ZMod 2) : ZMod 2 :=
  x 0 * x 4 + x 1 * x 5 + x 2 * x 6 + x 3 * x 3

/-- A vector in 𝔽₂⁷ is on Q(6,2) iff Q(x) = 0. -/
def onQ62 (x : Fin 7 → ZMod 2) : Prop := Q62form x = 0

instance (x : Fin 7 → ZMod 2) : Decidable (onQ62 x) :=
  inferInstanceAs (Decidable (Q62form x = 0))

/-! ## Enumerate the 63 points

We represent each nonzero vector of 𝔽₂⁷ by its index 1..127,
converting to/from bit vectors. -/

/-- Convert a natural number (1..127) to its bit vector in 𝔽₂⁷. -/
def toBitVec (n : Fin 128) : Fin 7 → ZMod 2 :=
  fun i => if n.val / (2 ^ i.val) % 2 = 0 then 0 else 1

/-- The list of indices (1..127) whose bit vectors lie on Q(6,2). -/
def q62Indices : List (Fin 128) :=
  (List.finRange 128).tail.filter fun n => Q62form (toBitVec n) = 0

/-- Q(6,2) has exactly 63 points. -/
theorem q62_card : q62Indices.length = 63 := by native_decide

/-! ## The parity-check matrix H

H is the 7 × 63 matrix whose j-th column is the j-th point of Q(6,2).
We define it by looking up the j-th element of q62Indices. -/

/-- Get the n-th point of Q(6,2) (as an index into Fin 128). -/
def q62Point (j : Fin 63) : Fin 128 :=
  q62Indices.get (j.cast (by native_decide))

/-- The parity-check matrix: H[i][j] = (j-th Q(6,2) point)[i-th coordinate]. -/
def Q62_H : Matrix (Fin 7) (Fin 63) (ZMod 2) :=
  fun i j => toBitVec (q62Point j) i

/-! ## The main theorem: self-orthogonality -/

/-- **Self-orthogonality of Q(6,2)**: H · Hᵀ = 0 over 𝔽₂.

    This is the CSS condition: it says that for any two coordinates i, k ∈ {0,...,6},
    the number of Q(6,2) points with both coordinate i and coordinate k equal to 1
    is even. Equivalently, the code C = ker(H) contains its dual C⊥ = rowspan(H).

    Consequence: C gives a CSS quantum code [[63, 49, 3]] with automorphism
    group G₂(2), the exceptional Chevalley group of order 12096. -/
theorem Q62_self_orthogonal : Q62_H * Q62_Hᵀ = (0 : Matrix (Fin 7) (Fin 7) (ZMod 2)) := by
  native_decide
