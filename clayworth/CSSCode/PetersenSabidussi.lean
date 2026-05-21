/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Petersen: Sabidussi coset graph

The **Petersen graph** (10 vertices, cubic) is Sab(A₅, S₃).
|A₅| = 60, |S₃| = 6, 60/6 = 10 vertices.
-/

import CSSCode.SabidussiWitness

/-! ## Graph definition -/

private def petersenEdges : Array (Fin 10 × Fin 10) := #[
  (0, 3),
  (0, 4),
  (0, 7),
  (1, 6),
  (1, 7),
  (1, 8),
  (2, 3),
  (2, 5),
  (2, 8),
  (3, 6),
  (4, 8),
  (4, 9),
  (5, 7),
  (5, 9),
  (6, 9)
]

private def petersenAdjBool (u v : Fin 10) : Bool :=
  petersenEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

def petGraph : SimpleGraph (Fin 10) where
  Adj u v := petersenAdjBool u v
  symm u v := by simp only [petersenAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [petersenAdjBool]; revert u; native_decide⟩

instance petersenDecAdj : DecidableRel petGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (petersenAdjBool u v))

/-! ## Generators -/

private def petersenGen1Fwd : Array (Fin 10) := #[
  1, 3, 5, 7, 8, 9, 0, 6, 2, 4
]
private def petersenGen1Inv : Array (Fin 10) := #[
  6, 0, 8, 1, 9, 2, 7, 3, 4, 5
]

private def petersenGen2Fwd : Array (Fin 10) := #[
  2, 4, 6, 3, 5, 1, 0, 8, 9, 7
]
private def petersenGen2Inv : Array (Fin 10) := #[
  6, 5, 0, 3, 1, 4, 2, 9, 7, 8
]

private theorem petersenGen1Fwd_size : petersenGen1Fwd.size = 10 := by native_decide
private theorem petersenGen1Inv_size : petersenGen1Inv.size = 10 := by native_decide
private theorem petersenGen2Fwd_size : petersenGen2Fwd.size = 10 := by native_decide
private theorem petersenGen2Inv_size : petersenGen2Inv.size = 10 := by native_decide

def petersenGen1 : Equiv.Perm (Fin 10) where
  toFun i := petersenGen1Fwd[i.val]'(by have := petersenGen1Fwd_size; omega)
  invFun i := petersenGen1Inv[i.val]'(by have := petersenGen1Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def petersenGen2 : Equiv.Perm (Fin 10) where
  toFun i := petersenGen2Fwd[i.val]'(by have := petersenGen2Fwd_size; omega)
  invFun i := petersenGen2Inv[i.val]'(by have := petersenGen2Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def petersenGens : Fin 2 → Equiv.Perm (Fin 10)
  | 0 => petersenGen1
  | 1 => petersenGen2

theorem petersenGen1_adj : ∀ u v : Fin 10,
    petGraph.Adj u v → petGraph.Adj (petersenGen1 u) (petersenGen1 v) := by
  native_decide

theorem petersenGen1_inv_adj : ∀ u v : Fin 10,
    petGraph.Adj u v → petGraph.Adj (petersenGen1.symm u) (petersenGen1.symm v) := by
  native_decide

theorem petersenGen2_adj : ∀ u v : Fin 10,
    petGraph.Adj u v → petGraph.Adj (petersenGen2 u) (petersenGen2 v) := by
  native_decide

theorem petersenGen2_inv_adj : ∀ u v : Fin 10,
    petGraph.Adj u v → petGraph.Adj (petersenGen2.symm u) (petersenGen2.symm v) := by
  native_decide

/-! ## BFS witness words -/

private def petersenWitnessWords : Fin 10 → List (Fin 4) := fun v =>
  match v with
  | 0 => []
  | 1 => [0]
  | 2 => [1]
  | 3 => [0, 0]
  | 4 => [0, 1]
  | 5 => [0, 3]
  | 6 => [2]
  | 7 => [2, 2]
  | 8 => [1, 2]
  | 9 => [0, 1, 2]

theorem petersenWitness_correct :
    ∀ v : Fin 10, applyWordGen petersenGens (petersenWitnessWords v) 0 = v := by
  native_decide

/-! ## Sabidussi witness -/

def petersenSabidussiWitness : SabidussiWitness petGraph where
  numGens := 2
  gens := petersenGens
  gens_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact petersenGen1_adj
    | ⟨1, _⟩ => exact petersenGen2_adj
  gens_inv_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact petersenGen1_inv_adj
    | ⟨1, _⟩ => exact petersenGen2_inv_adj
  base := 0
  witnessWords := petersenWitnessWords
  witness_correct := petersenWitness_correct

noncomputable def petersenSabidussiIso := petersenSabidussiWitness.toSabidussiIso
