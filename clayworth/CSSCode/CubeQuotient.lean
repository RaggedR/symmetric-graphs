/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Tetrahedron (K₄) as a Quotient of the Cube

The 3-cube Q₃ (8 vertices, cubic) has an antipodal involution v ↦ v̄
that identifies each vertex with its bitwise complement (distance 3).
The quotient by this Z₂ action is K₄, the complete graph on 4 vertices —
equivalently, the 1-skeleton of the tetrahedron.

## References

- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import CSSCode.QuotientGraph

/-! ## The cube graph (Q₃) -/

/-- Edge list for the 3-cube (12 edges). Vertices labeled by binary strings:
    0=000, 1=001, 2=010, 3=011, 4=100, 5=101, 6=110, 7=111.
    Adjacent iff Hamming distance 1. -/
private def cubeEdges : Array (Fin 8 × Fin 8) := #[
  (0,1), (0,2), (0,4), (1,3), (1,5), (2,3), (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)
]

private def cubeAdjBool (u v : Fin 8) : Bool :=
  cubeEdges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

/-- The 3-cube as a SimpleGraph on 8 vertices. -/
def cubeGraph : SimpleGraph (Fin 8) where
  Adj u v := cubeAdjBool u v
  symm u v := by simp only [cubeAdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [cubeAdjBool]; revert u; native_decide⟩

instance cubeDecAdj : DecidableRel cubeGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (cubeAdjBool u v))

/-! ## The antipodal map -/

/-- Antipodal map on Q₃: v ↦ 7-v (bitwise complement).
    Pairs: {0,7}, {1,6}, {2,5}, {3,4}. -/
private def cubeAntipodalData : Array (Fin 4) :=
  #[0, 1, 2, 3, 3, 2, 1, 0]

private theorem cubeAntipodalData_size : cubeAntipodalData.size = 8 := by native_decide

def cubeAntipodalMap (v : Fin 8) : Fin 4 :=
  have := cubeAntipodalData_size
  cubeAntipodalData[v.val]'(by omega)

/-! ## The tetrahedron (K₄) -/

/-- The **tetrahedron** (K₄), defined as the quotient of the 3-cube
    by the antipodal involution. -/
def tetrahedronGraph : SimpleGraph (Fin 4) :=
  cubeGraph.quotientGraph cubeAntipodalMap

instance tetrahedronDecAdj : DecidableRel tetrahedronGraph.Adj := by
  intro i j
  unfold tetrahedronGraph SimpleGraph.quotientGraph
  simp only
  exact instDecidableAnd

/-! ## Verification -/

/-- The cube is 3-regular. -/
theorem cube_regular :
    ∀ v : Fin 8, (Finset.univ.filter fun w => cubeGraph.Adj v w).card = 3 := by
  native_decide

/-- The tetrahedron is 3-regular (K₄). -/
theorem tetrahedron_regular :
    ∀ v : Fin 4, (Finset.univ.filter fun w => tetrahedronGraph.Adj v w).card = 3 := by
  native_decide

/-- The cube has 12 edges. -/
theorem cube_edges :
    (Finset.univ.filter fun p : Fin 8 × Fin 8 =>
      p.1 < p.2 ∧ cubeGraph.Adj p.1 p.2).card = 12 := by
  native_decide

/-- The tetrahedron has 6 edges (complete graph K₄). -/
theorem tetrahedron_edges :
    (Finset.univ.filter fun p : Fin 4 × Fin 4 =>
      p.1 < p.2 ∧ tetrahedronGraph.Adj p.1 p.2).card = 6 := by
  native_decide
