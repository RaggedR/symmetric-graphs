/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Petersen Graph as a Quotient of the Dodecahedron

The dodecahedron (20 vertices, cubic) has an antipodal involution that
identifies each vertex with its diametrically opposite vertex (distance 5).
The quotient by this Z₂ action is the Petersen graph (10 vertices, cubic).

## References

- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import CSSCode.QuotientGraph

/-! ## The dodecahedron graph -/

/-- Edge list for the dodecahedron (30 edges, 0-indexed). -/
private def dodecEdges : Array (Fin 20 × Fin 20) := #[
  (0,1), (0,4), (0,5), (1,2), (1,7), (2,3), (2,9), (3,4), (3,11), (4,13),
  (5,6), (5,14), (6,7), (6,15), (7,8), (8,9), (8,16), (9,10), (10,11), (10,17),
  (11,12), (12,13), (12,18), (13,14), (14,19), (15,16), (15,19), (16,17), (17,18), (18,19)
]

/-- Decidable adjacency for the dodecahedron. -/
private def dodecAdjBool (u v : Fin 20) : Bool :=
  dodecEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

/-- The dodecahedron as a SimpleGraph on 20 vertices. -/
def dodecahedronGraph : SimpleGraph (Fin 20) where
  Adj u v := dodecAdjBool u v
  symm u v := by simp only [dodecAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [dodecAdjBool]; revert u; native_decide⟩

instance dodecahedronDecAdj : DecidableRel dodecahedronGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (dodecAdjBool u v))

/-! ## The antipodal map -/

/-- Antipodal pairs on the dodecahedron (computed by BFS, distance 5).
    Maps each of 20 vertices to one of 10 equivalence classes. -/
private def antipodalData : Array (Fin 10) :=
  #[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 3, 4, 0, 1, 2]

private theorem antipodalData_size : antipodalData.size = 20 := by native_decide

/-- The antipodal quotient map: Fin 20 → Fin 10. -/
def antipodalMap (v : Fin 20) : Fin 10 :=
  have := antipodalData_size
  antipodalData[v.val]'(by omega)

/-! ## The Petersen graph -/

/-- The **Petersen graph**, defined as the quotient of the dodecahedron
    by the antipodal involution. -/
def petersenGraph : SimpleGraph (Fin 10) :=
  dodecahedronGraph.quotientGraph antipodalMap

instance petersenDecAdj : DecidableRel petersenGraph.Adj := by
  intro i j
  unfold petersenGraph SimpleGraph.quotientGraph
  simp only
  exact instDecidableAnd

/-! ## Verification -/

/-- The Petersen graph is 3-regular. -/
theorem petersen_regular :
    ∀ v : Fin 10, (Finset.univ.filter fun w => petersenGraph.Adj v w).card = 3 := by
  native_decide

/-- The dodecahedron is 3-regular. -/
theorem dodecahedron_regular :
    ∀ v : Fin 20, (Finset.univ.filter fun w => dodecahedronGraph.Adj v w).card = 3 := by
  native_decide

/-- The Petersen graph has 15 edges. -/
theorem petersen_edges :
    (Finset.univ.filter fun p : Fin 10 × Fin 10 =>
      p.1 < p.2 ∧ petersenGraph.Adj p.1 p.2).card = 15 := by
  native_decide

/-- The dodecahedron has 30 edges. -/
theorem dodecahedron_edges :
    (Finset.univ.filter fun p : Fin 20 × Fin 20 =>
      p.1 < p.2 ∧ dodecahedronGraph.Adj p.1 p.2).card = 30 := by
  native_decide

