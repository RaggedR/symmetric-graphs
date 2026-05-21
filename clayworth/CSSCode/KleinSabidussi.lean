/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Klein: Sabidussi coset graph

The **Klein quartic graph** (56 vertices, cubic) is Sab(PSL(2,7), C₃).
|PSL(2,7)| = 168, |C₃| = 3, 168/3 = 56 vertices.
-/

import CSSCode.SabidussiWitness

/-! ## Graph definition -/

private def kleinEdges : Array (Fin 56 × Fin 56) := #[
  (0, 11),
  (0, 22),
  (0, 49),
  (1, 7),
  (1, 14),
  (1, 28),
  (2, 19),
  (2, 43),
  (2, 52),
  (3, 12),
  (3, 36),
  (3, 45),
  (4, 21),
  (4, 32),
  (4, 42),
  (5, 18),
  (5, 24),
  (5, 35),
  (6, 33),
  (6, 38),
  (6, 50),
  (7, 25),
  (7, 51),
  (8, 21),
  (8, 33),
  (8, 35),
  (9, 31),
  (9, 47),
  (9, 52),
  (10, 16),
  (10, 26),
  (10, 43),
  (11, 39),
  (11, 44),
  (12, 15),
  (12, 34),
  (13, 19),
  (13, 40),
  (13, 53),
  (14, 27),
  (14, 37),
  (15, 22),
  (15, 42),
  (16, 38),
  (16, 54),
  (17, 25),
  (17, 30),
  (17, 50),
  (18, 46),
  (18, 51),
  (19, 29),
  (20, 33),
  (20, 39),
  (20, 47),
  (21, 54),
  (22, 40),
  (23, 38),
  (23, 45),
  (23, 52),
  (24, 29),
  (24, 47),
  (25, 48),
  (26, 28),
  (26, 41),
  (27, 31),
  (27, 55),
  (28, 44),
  (29, 49),
  (30, 40),
  (30, 45),
  (31, 36),
  (32, 37),
  (32, 53),
  (34, 46),
  (34, 54),
  (35, 48),
  (36, 51),
  (37, 43),
  (39, 55),
  (41, 46),
  (41, 49),
  (42, 55),
  (44, 50),
  (48, 53)
]

private def kleinAdjBool (u v : Fin 56) : Bool :=
  kleinEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

def kleinSabGraph : SimpleGraph (Fin 56) where
  Adj u v := kleinAdjBool u v
  symm u v := by simp only [kleinAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [kleinAdjBool]; revert u; native_decide⟩

instance kleinDecAdj : DecidableRel kleinSabGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (kleinAdjBool u v))

/-! ## Generators -/

private def kleinGen1Fwd : Array (Fin 56) := #[
  5, 1, 6, 2, 0, 4, 3, 14, 15, 16, 17, 18, 19, 20, 28, 29, 30, 31, 32, 33, 34, 22, 24, 23, 21, 27, 25, 26, 7, 8,
  9, 10, 11, 12, 13, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 35, 36, 37, 38, 39, 40, 41
]
private def kleinGen1Inv : Array (Fin 56) := #[
  4, 1, 3, 6, 5, 0, 2, 28, 29, 30, 31, 32, 33, 34, 7, 8, 9, 10, 11, 12, 13, 24, 21, 23, 22, 26, 27, 25, 14, 15,
  16, 17, 18, 19, 20, 49, 50, 51, 52, 53, 54, 55, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48
]

private def kleinGen2Fwd : Array (Fin 56) := #[
  7, 8, 9, 10, 11, 12, 13, 21, 22, 23, 24, 25, 26, 27, 33, 28, 29, 32, 34, 31, 30, 0, 1, 2, 3, 4, 5, 6, 35, 36,
  37, 38, 39, 40, 41, 15, 16, 20, 19, 17, 14, 18, 44, 47, 48, 43, 46, 45, 42, 51, 53, 54, 52, 55, 49, 50
]
private def kleinGen2Inv : Array (Fin 56) := #[
  21, 22, 23, 24, 25, 26, 27, 0, 1, 2, 3, 4, 5, 6, 40, 35, 36, 39, 41, 38, 37, 7, 8, 9, 10, 11, 12, 13, 15, 16,
  20, 19, 17, 14, 18, 28, 29, 30, 31, 32, 33, 34, 48, 45, 42, 47, 46, 43, 44, 54, 55, 49, 52, 50, 51, 53
]

private theorem kleinGen1Fwd_size : kleinGen1Fwd.size = 56 := by native_decide
private theorem kleinGen1Inv_size : kleinGen1Inv.size = 56 := by native_decide
private theorem kleinGen2Fwd_size : kleinGen2Fwd.size = 56 := by native_decide
private theorem kleinGen2Inv_size : kleinGen2Inv.size = 56 := by native_decide

def kleinGen1 : Equiv.Perm (Fin 56) where
  toFun i := kleinGen1Fwd[i.val]'(by have := kleinGen1Fwd_size; omega)
  invFun i := kleinGen1Inv[i.val]'(by have := kleinGen1Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def kleinGen2 : Equiv.Perm (Fin 56) where
  toFun i := kleinGen2Fwd[i.val]'(by have := kleinGen2Fwd_size; omega)
  invFun i := kleinGen2Inv[i.val]'(by have := kleinGen2Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

def kleinGens : Fin 2 → Equiv.Perm (Fin 56)
  | 0 => kleinGen1
  | 1 => kleinGen2

theorem kleinGen1_adj : ∀ u v : Fin 56,
    kleinSabGraph.Adj u v → kleinSabGraph.Adj (kleinGen1 u) (kleinGen1 v) := by
  native_decide

theorem kleinGen1_inv_adj : ∀ u v : Fin 56,
    kleinSabGraph.Adj u v → kleinSabGraph.Adj (kleinGen1.symm u) (kleinGen1.symm v) := by
  native_decide

theorem kleinGen2_adj : ∀ u v : Fin 56,
    kleinSabGraph.Adj u v → kleinSabGraph.Adj (kleinGen2 u) (kleinGen2 v) := by
  native_decide

theorem kleinGen2_inv_adj : ∀ u v : Fin 56,
    kleinSabGraph.Adj u v → kleinSabGraph.Adj (kleinGen2.symm u) (kleinGen2.symm v) := by
  native_decide

/-! ## BFS witness words -/

private def kleinWitnessWords : Fin 56 → List (Fin 4) := fun v =>
  match v with
  | 0 => []
  | 1 => [3, 0, 1]
  | 2 => [3, 2, 1, 0]
  | 3 => [3, 2, 1]
  | 4 => [2]
  | 5 => [0]
  | 6 => [0, 3, 2, 1]
  | 7 => [1]
  | 8 => [3, 0, 3]
  | 9 => [3, 2, 1, 0, 1]
  | 10 => [3, 2, 3]
  | 11 => [2, 1]
  | 12 => [0, 1]
  | 13 => [0, 3, 2, 3]
  | 14 => [1, 0]
  | 15 => [1, 2, 3]
  | 16 => [1, 2, 3, 0, 3]
  | 17 => [2, 1, 2, 3]
  | 18 => [2, 1, 0]
  | 19 => [0, 1, 0]
  | 20 => [0, 3, 2, 3, 0]
  | 21 => [3]
  | 22 => [3, 0]
  | 23 => [3, 2, 1, 0, 3]
  | 24 => [3, 2]
  | 25 => [2, 3]
  | 26 => [0, 3]
  | 27 => [0, 3, 2]
  | 28 => [1, 2]
  | 29 => [1, 2, 3, 0]
  | 30 => [0, 3, 2, 3, 0, 1]
  | 31 => [0, 1, 0, 1]
  | 32 => [2, 1, 2]
  | 33 => [0, 1, 2]
  | 34 => [2, 1, 0, 1]
  | 35 => [1, 2, 1]
  | 36 => [1, 2, 3, 0, 1]
  | 37 => [0, 3, 2, 3, 0, 3]
  | 38 => [0, 1, 0, 3]
  | 39 => [2, 1, 2, 1]
  | 40 => [1, 0, 3]
  | 41 => [2, 1, 0, 3]
  | 42 => [1, 2, 1, 0]
  | 43 => [1, 0, 3, 0, 3]
  | 44 => [1, 2, 1, 0, 1]
  | 45 => [0, 1, 0, 3, 0]
  | 46 => [2, 1, 2, 1, 0]
  | 47 => [1, 0, 3, 0]
  | 48 => [1, 2, 1, 0, 3]
  | 49 => [1, 2, 1, 2]
  | 50 => [1, 0, 3, 0, 3, 0]
  | 51 => [1, 0, 3, 2, 3]
  | 52 => [0, 1, 0, 3, 2]
  | 53 => [2, 1, 2, 1, 2]
  | 54 => [1, 0, 3, 2]
  | 55 => [2, 1, 0, 3, 2]
  | _ => []  -- unreachable, needed for exhaustiveness

theorem kleinWitness_correct :
    ∀ v : Fin 56, applyWordGen kleinGens (kleinWitnessWords v) 0 = v := by
  native_decide

/-! ## Sabidussi witness -/

def kleinSabidussiWitness : SabidussiWitness kleinSabGraph where
  numGens := 2
  gens := kleinGens
  gens_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact kleinGen1_adj
    | ⟨1, _⟩ => exact kleinGen2_adj
  gens_inv_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact kleinGen1_inv_adj
    | ⟨1, _⟩ => exact kleinGen2_inv_adj
  base := 0
  witnessWords := kleinWitnessWords
  witness_correct := kleinWitness_correct

noncomputable def kleinSabidussiIso := kleinSabidussiWitness.toSabidussiIso
