/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Heawood: Sabidussi coset graph

The **Heawood graph** (14 vertices, cubic) is Sab(PGL(2,7), H₂₄).
|PGL(2,7)| = 336, |H₂₄| = 24, 336/24 = 14 vertices.
-/

import CSSCode.SabidussiWitness

/-! ## Graph definition -/

private def heawoodEdges : Array (Fin 14 × Fin 14) := #[
  (0, 3),
  (0, 8),
  (0, 10),
  (1, 4),
  (1, 12),
  (1, 13),
  (2, 5),
  (2, 6),
  (2, 8),
  (3, 7),
  (3, 12),
  (4, 8),
  (4, 11),
  (5, 9),
  (5, 12),
  (6, 7),
  (6, 13),
  (7, 11),
  (9, 10),
  (9, 11),
  (10, 13)
]

private def heawoodAdjBool (u v : Fin 14) : Bool :=
  heawoodEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

def heawGraph : SimpleGraph (Fin 14) where
  Adj u v := heawoodAdjBool u v
  symm u v := by simp only [heawoodAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [heawoodAdjBool]; revert u; native_decide⟩

instance heawoodDecAdj : DecidableRel heawGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (heawoodAdjBool u v))

/-! ## Generators -/

private def heawoodGen1Fwd : Array (Fin 14) := #[
  5, 0, 1, 2, 3, 4, 13, 6, 12, 11, 9, 7, 8, 10
]
private def heawoodGen1Inv : Array (Fin 14) := #[
  1, 2, 3, 4, 5, 0, 7, 11, 12, 10, 13, 9, 8, 6
]

private def heawoodGen2Fwd : Array (Fin 14) := #[
  7, 1, 2, 11, 13, 8, 5, 9, 6, 0, 3, 10, 4, 12
]
private def heawoodGen2Inv : Array (Fin 14) := #[
  9, 1, 2, 10, 12, 6, 8, 0, 5, 7, 11, 3, 13, 4
]

private theorem heawoodGen1Fwd_size : heawoodGen1Fwd.size = 14 := by native_decide
private theorem heawoodGen1Inv_size : heawoodGen1Inv.size = 14 := by native_decide
private theorem heawoodGen2Fwd_size : heawoodGen2Fwd.size = 14 := by native_decide
private theorem heawoodGen2Inv_size : heawoodGen2Inv.size = 14 := by native_decide

def heawoodGen1 : Equiv.Perm (Fin 14) where
  toFun i := heawoodGen1Fwd[i.val]'(by have := heawoodGen1Fwd_size; omega)
  invFun i := heawoodGen1Inv[i.val]'(by have := heawoodGen1Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def heawoodGen2 : Equiv.Perm (Fin 14) where
  toFun i := heawoodGen2Fwd[i.val]'(by have := heawoodGen2Fwd_size; omega)
  invFun i := heawoodGen2Inv[i.val]'(by have := heawoodGen2Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def heawoodGens : Fin 2 → Equiv.Perm (Fin 14)
  | 0 => heawoodGen1
  | 1 => heawoodGen2

theorem heawoodGen1_adj : ∀ u v : Fin 14,
    heawGraph.Adj u v → heawGraph.Adj (heawoodGen1 u) (heawoodGen1 v) := by
  native_decide

theorem heawoodGen1_inv_adj : ∀ u v : Fin 14,
    heawGraph.Adj u v → heawGraph.Adj (heawoodGen1.symm u) (heawoodGen1.symm v) := by
  native_decide

theorem heawoodGen2_adj : ∀ u v : Fin 14,
    heawGraph.Adj u v → heawGraph.Adj (heawoodGen2 u) (heawoodGen2 v) := by
  native_decide

theorem heawoodGen2_inv_adj : ∀ u v : Fin 14,
    heawGraph.Adj u v → heawGraph.Adj (heawoodGen2.symm u) (heawoodGen2.symm v) := by
  native_decide

/-! ## BFS witness words -/

private def heawoodWitnessWords : Fin 14 → List (Fin 4) := fun v =>
  match v with
  | 0 => []
  | 1 => [2]
  | 2 => [2, 2]
  | 3 => [0, 0, 0]
  | 4 => [0, 0]
  | 5 => [0]
  | 6 => [0, 3]
  | 7 => [1]
  | 8 => [0, 1]
  | 9 => [3]
  | 10 => [3, 2]
  | 11 => [1, 2]
  | 12 => [0, 0, 3]
  | 13 => [0, 0, 1]
  | _ => []  -- unreachable, needed for exhaustiveness

theorem heawoodWitness_correct :
    ∀ v : Fin 14, applyWordGen heawoodGens (heawoodWitnessWords v) 0 = v := by
  native_decide

/-! ## Sabidussi witness -/

def heawoodSabidussiWitness : SabidussiWitness heawGraph where
  numGens := 2
  gens := heawoodGens
  gens_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact heawoodGen1_adj
    | ⟨1, _⟩ => exact heawoodGen2_adj
  gens_inv_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact heawoodGen1_inv_adj
    | ⟨1, _⟩ => exact heawoodGen2_inv_adj
  base := 0
  witnessWords := heawoodWitnessWords
  witness_correct := heawoodWitness_correct

noncomputable def heawoodSabidussiIso := heawoodSabidussiWitness.toSabidussiIso
