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
import Mathlib.Combinatorics.SimpleGraph.Basic

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

/-! ## The Langer graph: collinearity graph of GH(2,2)

The **Langer graph** is the collinearity graph of the split Cayley hexagon GH(2,2)
on the 63 points of Q(6,2). Two points x, y are adjacent iff they are joined by a
line of the hexagon.

Algebraically, each point x = (x₀,...,x₆) on Q(6,2) corresponds to the Zorn vector
matrix M_x = [[x₃, (x₀,x₁,x₂)], [(x₄,x₅,x₆), x₃]] in the split octonion algebra
over 𝔽₂. Two points are collinear in GH(2,2) iff M_x · M_y = 0.

The Zorn multiplication (Springer–Veldkamp, *Octonions, Jordan Algebras and
Exceptional Groups*, §1.8) for imaginary elements with α = a·b, γ = c·d gives:

  M_x · M_y = [[αγ + a·d,  γa + αc + b×d],
               [γb + αd + a×c,  αγ + b·c]]

where a = (x₀,x₁,x₂), b = (x₄,x₅,x₆), c = (y₀,y₁,y₂), d = (y₄,y₅,y₆),
and × is the cross product in 𝔽₂³. Setting all eight components to zero gives
eight conditions (two scalar, six vector) characterising collinearity.

The resulting graph is 6-regular with 189 edges — the unique distance-regular graph
with intersection array {6, 4, 4; 1, 1, 3}. -/

/-- Adjacency test for the Langer graph via the Zorn product criterion.
    Returns `true` iff `i ≠ j` and all eight entries of M_x · M_y vanish. -/
def langerAdjBool (i j : Fin 63) : Bool :=
  if i == j then false else
  let x := toBitVec (q62Point i)
  let y := toBitVec (q62Point j)
  let α := x 3; let γ := y 3
  -- Scalar entries: αγ + a·d = 0 and αγ + b·c = 0
  (α * γ + x 0 * y 4 + x 1 * y 5 + x 2 * y 6 == 0) &&
  (α * γ + x 4 * y 0 + x 5 * y 1 + x 6 * y 2 == 0) &&
  -- (1,2) vector: γa + αc + b×d = 0 (three components)
  (γ * x 0 + α * y 0 + x 5 * y 6 + x 6 * y 5 == 0) &&
  (γ * x 1 + α * y 1 + x 6 * y 4 + x 4 * y 6 == 0) &&
  (γ * x 2 + α * y 2 + x 4 * y 5 + x 5 * y 4 == 0) &&
  -- (2,1) vector: γb + αd + a×c = 0 (three components)
  (γ * x 4 + α * y 4 + x 1 * y 2 + x 2 * y 1 == 0) &&
  (γ * x 5 + α * y 5 + x 2 * y 0 + x 0 * y 2 == 0) &&
  (γ * x 6 + α * y 6 + x 0 * y 1 + x 1 * y 0 == 0)

/-- The **Langer graph**: collinearity graph of the split Cayley hexagon GH(2,2)
    on the 63 points of Q(6,2).

    63 vertices, 6-regular, distance-regular with intersection array {6, 4, 4; 1, 1, 3}.
    Aut ≅ G₂(2), the exceptional Chevalley group of order 12096. -/
def langerSimpleGraph : SimpleGraph (Fin 63) where
  Adj i j := langerAdjBool i j
  symm i j := by simp only [langerAdjBool]; revert i j; native_decide
  loopless := ⟨fun i => by simp only [langerAdjBool]; revert i; native_decide⟩

instance langerDecAdj : DecidableRel langerSimpleGraph.Adj :=
  fun i j => inferInstanceAs (Decidable (langerAdjBool i j))

/-- The Langer graph is 6-regular. -/
theorem langer_regular :
    ∀ v : Fin 63, (Finset.univ.filter fun w => langerSimpleGraph.Adj v w).card = 6 := by
  native_decide

/-- The Langer graph has 189 undirected edges. -/
theorem langer_edges :
    (Finset.univ.filter fun p : Fin 63 × Fin 63 =>
      p.1 < p.2 ∧ langerSimpleGraph.Adj p.1 p.2).card = 189 := by
  native_decide

/-! ## Primitivity

The Langer graph is **primitive**: Aut(Langer) ≅ G₂(2) (order 12096) acts
transitively on the 63 vertices with no non-trivial block system. This means
there is no symmetric quotient — the Langer graph cannot be expressed as a
non-trivial covering of a smaller graph.

This is consistent with the Langer graph being the collinearity graph of
GH(2,2): the 63 points of the generalized hexagon carry a primitive action
of G₂(2). -/
