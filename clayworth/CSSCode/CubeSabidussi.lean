/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Cube Graph is a Sabidussi Coset Graph

The **cube graph** Q₃ (8 vertices, cubic) is a Sabidussi coset graph:
`Q₃ ≃g Sab(S₄, C₃, D)`. We prove this by exhibiting two generators of
S₄ ≅ Aut(Q₃)⁺ (the rotation group) acting vertex-transitively on Q₃ and
applying the Sabidussi representation theorem.

**Generators** (computed by GAP):
- σ₁: 90° rotation around z-axis, order 4: (0,2,6,4)(1,3,7,5)
- σ₂: 120° rotation around (1,1,1) diagonal, order 3: (1,2,4)(3,6,5)
Together ⟨σ₁, σ₂⟩ ≅ S₄ (order 24), stabilizer of vertex 0 ≅ C₃ (order 3).

## References

- Sabidussi, *Vertex-transitive graphs*, Monatsh. Math. 68 (1964)
-/

import CSSCode.SabidussiWitness
import CSSCode.CubeQuotient

/-! ## S₄ generators on the 8 cube vertices -/

private def cubeGen1Fwd : Array (Fin 8) := #[2, 3, 6, 7, 0, 1, 4, 5]
private def cubeGen1Inv : Array (Fin 8) := #[4, 5, 0, 1, 6, 7, 2, 3]
private def cubeGen2Fwd : Array (Fin 8) := #[0, 2, 4, 6, 1, 3, 5, 7]
private def cubeGen2Inv : Array (Fin 8) := #[0, 4, 1, 5, 2, 6, 3, 7]

private theorem cubeGen1Fwd_size : cubeGen1Fwd.size = 8 := by native_decide
private theorem cubeGen1Inv_size : cubeGen1Inv.size = 8 := by native_decide
private theorem cubeGen2Fwd_size : cubeGen2Fwd.size = 8 := by native_decide
private theorem cubeGen2Inv_size : cubeGen2Inv.size = 8 := by native_decide

/-- σ₁: 90° rotation around z-axis, order 4.
    In cycle notation: (0,2,6,4)(1,3,7,5). -/
def cubeGen1 : Equiv.Perm (Fin 8) where
  toFun i := cubeGen1Fwd[i.val]'(by have := cubeGen1Fwd_size; omega)
  invFun i := cubeGen1Inv[i.val]'(by have := cubeGen1Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

/-- σ₂: 120° rotation around the (1,1,1)↔(0,0,0) diagonal, order 3.
    In cycle notation: (1,2,4)(3,6,5). Fixes vertices 0 and 7. -/
def cubeGen2 : Equiv.Perm (Fin 8) where
  toFun i := cubeGen2Fwd[i.val]'(by have := cubeGen2Fwd_size; omega)
  invFun i := cubeGen2Inv[i.val]'(by have := cubeGen2Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

/-! ## Adjacency preservation -/

/-- σ₁ preserves cube adjacency. -/
theorem cubeGen1_adj : ∀ u v : Fin 8,
    cubeGraph.Adj u v → cubeGraph.Adj (cubeGen1 u) (cubeGen1 v) := by
  native_decide

/-- σ₁⁻¹ preserves cube adjacency. -/
theorem cubeGen1_inv_adj : ∀ u v : Fin 8,
    cubeGraph.Adj u v → cubeGraph.Adj (cubeGen1.symm u) (cubeGen1.symm v) := by
  native_decide

/-- σ₂ preserves cube adjacency. -/
theorem cubeGen2_adj : ∀ u v : Fin 8,
    cubeGraph.Adj u v → cubeGraph.Adj (cubeGen2 u) (cubeGen2 v) := by
  native_decide

/-- σ₂⁻¹ preserves cube adjacency. -/
theorem cubeGen2_inv_adj : ∀ u v : Fin 8,
    cubeGraph.Adj u v → cubeGraph.Adj (cubeGen2.symm u) (cubeGen2.symm v) := by
  native_decide

/-! ## BFS witness words -/

/-- BFS witness words from vertex 0 to each vertex.
    Word encoding: 0=σ₁, 1=σ₁⁻¹, 2=σ₂, 3=σ₂⁻¹. -/
-- genOrInv encoding: 0=g1, 1=g2, 2=g1⁻¹, 3=g2⁻¹
private def cubeWitnessWords : Fin 8 → List (Fin 4) := fun v =>
  match v with
  | 0 => []           -- identity
  | 1 => [0, 3]       -- g1, g2⁻¹
  | 2 => [0]          -- g1
  | 3 => [0, 0, 3]    -- g1, g1, g2⁻¹
  | 4 => [2]          -- g1⁻¹
  | 5 => [0, 0, 1]    -- g1, g1, g2
  | 6 => [0, 0]       -- g1, g1
  | 7 => [0, 0, 1, 2] -- g1, g1, g2, g1⁻¹

/-- The two generators as a `Fin 2`-indexed family. -/
def cubeGens : Fin 2 → Equiv.Perm (Fin 8)
  | 0 => cubeGen1
  | 1 => cubeGen2

/-- Witness words are correct: each maps vertex 0 to the target. -/
theorem cubeWitness_correct :
    ∀ v : Fin 8, applyWordGen cubeGens (cubeWitnessWords v) 0 = v := by
  native_decide

/-! ## The Sabidussi witness and the main theorem -/

/-- Sabidussi witness for the cube graph: S₄ acting via rotations. -/
def cubeSabidussiWitness : SabidussiWitness cubeGraph where
  numGens := 2
  gens := cubeGens
  gens_preserve_adj := by
    intro i
    match i with
    | ⟨0, _⟩ => exact cubeGen1_adj
    | ⟨1, _⟩ => exact cubeGen2_adj
  gens_inv_preserve_adj := by
    intro i
    match i with
    | ⟨0, _⟩ => exact cubeGen1_inv_adj
    | ⟨1, _⟩ => exact cubeGen2_inv_adj
  base := 0
  witnessWords := cubeWitnessWords
  witness_correct := cubeWitness_correct

/-- **The cube graph is a Sabidussi coset graph.**

The rotation group S₄ = ⟨σ₁, σ₂⟩ acts vertex-transitively on Q₃.
The stabilizer of vertex 0 is C₃ = ⟨σ₂⟩ (order 3), giving
  Q₃ ≃g Sab(S₄, C₃, D)
with 8 = |S₄|/|C₃| = 24/3 vertices. -/
noncomputable def cubeSabidussiIso := cubeSabidussiWitness.toSabidussiIso
