/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Dodec: Sabidussi coset graph

The **dodecahedron** (20 vertices, cubic) is Sab(A₅, C₃).
|A₅| = 60, |C₃| = 3, 60/3 = 20 vertices.
-/

import CSSCode.SabidussiWitness

/-! ## Graph definition -/

private def dodecEdges : Array (Fin 20 × Fin 20) := #[
  (0, 8),
  (0, 12),
  (0, 16),
  (1, 4),
  (1, 11),
  (1, 14),
  (2, 7),
  (2, 10),
  (2, 18),
  (3, 5),
  (3, 15),
  (3, 17),
  (4, 13),
  (4, 16),
  (5, 8),
  (5, 19),
  (6, 11),
  (6, 15),
  (6, 18),
  (7, 9),
  (7, 12),
  (8, 13),
  (9, 14),
  (9, 16),
  (10, 15),
  (10, 19),
  (11, 17),
  (12, 19),
  (13, 17),
  (14, 18)
]

private def dodecAdjBool (u v : Fin 20) : Bool :=
  dodecEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

def dodecGraph : SimpleGraph (Fin 20) where
  Adj u v := dodecAdjBool u v
  symm u v := by simp only [dodecAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [dodecAdjBool]; revert u; native_decide⟩

instance dodecDecAdj : DecidableRel dodecGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (dodecAdjBool u v))

/-! ## Generators -/

private def dodecGen1Fwd : Array (Fin 20) := #[
  4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 1, 0, 3, 2, 13, 12, 15, 14
]
private def dodecGen1Inv : Array (Fin 20) := #[
  13, 12, 15, 14, 0, 1, 2, 3, 4, 5, 6, 7, 17, 16, 19, 18, 8, 9, 10, 11
]

private def dodecGen2Fwd : Array (Fin 20) := #[
  5, 7, 6, 4, 12, 13, 14, 15, 8, 10, 11, 9, 3, 0, 2, 1, 19, 16, 18, 17
]
private def dodecGen2Inv : Array (Fin 20) := #[
  13, 15, 14, 12, 3, 0, 2, 1, 8, 11, 9, 10, 4, 5, 6, 7, 17, 19, 18, 16
]

private theorem dodecGen1Fwd_size : dodecGen1Fwd.size = 20 := by native_decide
private theorem dodecGen1Inv_size : dodecGen1Inv.size = 20 := by native_decide
private theorem dodecGen2Fwd_size : dodecGen2Fwd.size = 20 := by native_decide
private theorem dodecGen2Inv_size : dodecGen2Inv.size = 20 := by native_decide

def dodecGen1 : Equiv.Perm (Fin 20) where
  toFun i := dodecGen1Fwd[i.val]'(by have := dodecGen1Fwd_size; omega)
  invFun i := dodecGen1Inv[i.val]'(by have := dodecGen1Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def dodecGen2 : Equiv.Perm (Fin 20) where
  toFun i := dodecGen2Fwd[i.val]'(by have := dodecGen2Fwd_size; omega)
  invFun i := dodecGen2Inv[i.val]'(by have := dodecGen2Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def dodecGens : Fin 2 → Equiv.Perm (Fin 20)
  | 0 => dodecGen1
  | 1 => dodecGen2

theorem dodecGen1_adj : ∀ u v : Fin 20,
    dodecGraph.Adj u v → dodecGraph.Adj (dodecGen1 u) (dodecGen1 v) := by
  native_decide

theorem dodecGen1_inv_adj : ∀ u v : Fin 20,
    dodecGraph.Adj u v → dodecGraph.Adj (dodecGen1.symm u) (dodecGen1.symm v) := by
  native_decide

theorem dodecGen2_adj : ∀ u v : Fin 20,
    dodecGraph.Adj u v → dodecGraph.Adj (dodecGen2 u) (dodecGen2 v) := by
  native_decide

theorem dodecGen2_inv_adj : ∀ u v : Fin 20,
    dodecGraph.Adj u v → dodecGraph.Adj (dodecGen2.symm u) (dodecGen2.symm v) := by
  native_decide

/-! ## BFS witness words -/

private def dodecWitnessWords : Fin 20 → List (Fin 4) := fun v =>
  match v with
  | 0 => []
  | 1 => [1, 2]
  | 2 => [0, 3, 2, 1]
  | 3 => [0, 3]
  | 4 => [0]
  | 5 => [1]
  | 6 => [0, 3, 2, 3]
  | 7 => [0, 3, 0]
  | 8 => [0, 0]
  | 9 => [1, 0]
  | 10 => [1, 0, 1]
  | 11 => [1, 0, 3]
  | 12 => [0, 1]
  | 13 => [2]
  | 14 => [0, 3, 2]
  | 15 => [1, 2, 3]
  | 16 => [2, 2]
  | 17 => [0, 1, 2]
  | 18 => [1, 0, 1, 0]
  | 19 => [2, 2, 1]
  | _ => []  -- unreachable, needed for exhaustiveness

theorem dodecWitness_correct :
    ∀ v : Fin 20, applyWordGen dodecGens (dodecWitnessWords v) 0 = v := by
  native_decide

/-! ## Sabidussi witness -/

def dodecSabidussiWitness : SabidussiWitness dodecGraph where
  numGens := 2
  gens := dodecGens
  gens_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact dodecGen1_adj
    | ⟨1, _⟩ => exact dodecGen2_adj
  gens_inv_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact dodecGen1_inv_adj
    | ⟨1, _⟩ => exact dodecGen2_inv_adj
  base := 0
  witnessWords := dodecWitnessWords
  witness_correct := dodecWitness_correct

noncomputable def dodecSabidussiIso := dodecSabidussiWitness.toSabidussiIso
