/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Mk: Sabidussi coset graph

The **Möbius-Kantor graph** (16 vertices, cubic) is Sab(GL(2,3), C₃).
|GL(2,3)| = 48, |C₃| = 3, 48/3 = 16 vertices.
-/

import CSSCode.SabidussiWitness

/-! ## Graph definition -/

private def mkEdges : Array (Fin 16 × Fin 16) := #[
  (0, 5),
  (0, 11),
  (0, 13),
  (1, 6),
  (1, 8),
  (1, 14),
  (2, 7),
  (2, 9),
  (2, 15),
  (3, 4),
  (3, 10),
  (3, 12),
  (4, 9),
  (4, 13),
  (5, 10),
  (5, 14),
  (6, 11),
  (6, 15),
  (7, 8),
  (7, 12),
  (8, 13),
  (9, 14),
  (10, 15),
  (11, 12)
]

private def mkAdjBool (u v : Fin 16) : Bool :=
  mkEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

def mkSabGraph : SimpleGraph (Fin 16) where
  Adj u v := mkAdjBool u v
  symm u v := by simp only [mkAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [mkAdjBool]; revert u; native_decide⟩

instance mkDecAdj : DecidableRel mkSabGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (mkAdjBool u v))

/-! ## Generators -/

private def mkGen1Fwd : Array (Fin 16) := #[
  1, 0, 3, 2, 7, 6, 5, 4, 13, 12, 15, 14, 9, 8, 11, 10
]
private def mkGen1Inv : Array (Fin 16) := #[
  1, 0, 3, 2, 7, 6, 5, 4, 13, 12, 15, 14, 9, 8, 11, 10
]

private def mkGen2Fwd : Array (Fin 16) := #[
  6, 7, 4, 5, 14, 15, 12, 13, 8, 9, 10, 11, 0, 1, 2, 3
]
private def mkGen2Inv : Array (Fin 16) := #[
  12, 13, 14, 15, 2, 3, 0, 1, 8, 9, 10, 11, 6, 7, 4, 5
]

private theorem mkGen1Fwd_size : mkGen1Fwd.size = 16 := by native_decide
private theorem mkGen1Inv_size : mkGen1Inv.size = 16 := by native_decide
private theorem mkGen2Fwd_size : mkGen2Fwd.size = 16 := by native_decide
private theorem mkGen2Inv_size : mkGen2Inv.size = 16 := by native_decide

def mkGen1 : Equiv.Perm (Fin 16) where
  toFun i := mkGen1Fwd[i.val]'(by have := mkGen1Fwd_size; omega)
  invFun i := mkGen1Inv[i.val]'(by have := mkGen1Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def mkGen2 : Equiv.Perm (Fin 16) where
  toFun i := mkGen2Fwd[i.val]'(by have := mkGen2Fwd_size; omega)
  invFun i := mkGen2Inv[i.val]'(by have := mkGen2Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def mkGens : Fin 2 → Equiv.Perm (Fin 16)
  | 0 => mkGen1
  | 1 => mkGen2

theorem mkGen1_adj : ∀ u v : Fin 16,
    mkSabGraph.Adj u v → mkSabGraph.Adj (mkGen1 u) (mkGen1 v) := by
  native_decide

theorem mkGen1_inv_adj : ∀ u v : Fin 16,
    mkSabGraph.Adj u v → mkSabGraph.Adj (mkGen1.symm u) (mkGen1.symm v) := by
  native_decide

theorem mkGen2_adj : ∀ u v : Fin 16,
    mkSabGraph.Adj u v → mkSabGraph.Adj (mkGen2 u) (mkGen2 v) := by
  native_decide

theorem mkGen2_inv_adj : ∀ u v : Fin 16,
    mkSabGraph.Adj u v → mkSabGraph.Adj (mkGen2.symm u) (mkGen2.symm v) := by
  native_decide

/-! ## BFS witness words -/

private def mkWitnessWords : Fin 16 → List (Fin 4) := fun v =>
  match v with
  | 0 => []
  | 1 => [0]
  | 2 => [0, 1, 0, 3]
  | 3 => [1, 0, 3]
  | 4 => [0, 1, 0]
  | 5 => [1, 0]
  | 6 => [1]
  | 7 => [0, 1]
  | 8 => [0, 3, 0]
  | 9 => [3, 0]
  | 10 => [1, 0, 1, 0]
  | 11 => [0, 1, 0, 1, 0]
  | 12 => [3]
  | 13 => [0, 3]
  | 14 => [0, 1, 0, 1]
  | 15 => [1, 0, 1]
  | _ => []  -- unreachable, needed for exhaustiveness

theorem mkWitness_correct :
    ∀ v : Fin 16, applyWordGen mkGens (mkWitnessWords v) 0 = v := by
  native_decide

/-! ## Sabidussi witness -/

def mkSabidussiWitness : SabidussiWitness mkSabGraph where
  numGens := 2
  gens := mkGens
  gens_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact mkGen1_adj
    | ⟨1, _⟩ => exact mkGen2_adj
  gens_inv_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact mkGen1_inv_adj
    | ⟨1, _⟩ => exact mkGen2_inv_adj
  base := 0
  witnessWords := mkWitnessWords
  witness_correct := mkWitness_correct

noncomputable def mkSabidussiIso := mkSabidussiWitness.toSabidussiIso
