/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Dodecahedron is a Double Cover of the Petersen Graph

The dodecahedron is a double cover of the Petersen graph. This covering
arises algebraically from the Sabidussi framework:

  Dodecahedron = Sab(A₅, C₃, D)   — |A₅|/|C₃| = 60/3 = 20 vertices
  Petersen     = Sab(A₅, S₃, D')  — |A₅|/|S₃| = 60/6 = 10 vertices

where C₃ ≤ S₃ are both subgroups of A₅ (the icosahedral rotation group
acting on 20 vertices). The generalized coset projection A₅/C₃ → A₅/S₃
(`proj_adj_of_subset` with D ⊆ D') gives the covering map.

**Subtlety**: The global antipodal involution (swapping ALL distance-5
pairs simultaneously) is NOT in A₅ — it lives in the full automorphism
group Aut(Dodec) ≅ A₅ × Z₂. However, A₅ does contain elements that
swap individual antipodal pairs (e.g., σ₃ swaps 0↔6), and the setwise
stabilizer K = Stab_{A₅}({0,6}) ≅ S₃ has C₃ ≤ K with [K:C₃] = 2.
The 10 cosets A₅/K are exactly the 10 antipodal pairs.

## References

- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import CSSCode.DodecSabidussi
import CSSCode.PetersenSabidussi
import CSSCode.CosetGraph

/-! ## The covering map

The covering map φ identifies antipodal pairs: each pair of vertices
at distance 5 in the dodecahedron maps to the same Petersen vertex.

  Fibers: {0,6}→0, {1,19}→1, {2,13}→2, {3,9}→9, {4,10}→8,
          {5,14}→6, {7,17}→5, {8,18}→3, {11,12}→7, {15,16}→4 -/

private def coveringData : Array (Fin 10) :=
  #[0, 1, 2, 9, 8, 6, 0, 5, 3, 9, 8, 7, 7, 2, 6, 4, 4, 5, 3, 1]

private theorem coveringData_size : coveringData.size = 20 := by native_decide

/-- The covering map φ : dodecahedron → Petersen, identifying antipodal pairs. -/
def dodecCoveringMap (v : Fin 20) : Fin 10 :=
  coveringData[v.val]'(by have := coveringData_size; omega)

/-! ## Graph homomorphism properties -/

/-- φ is a graph homomorphism: adjacent dodecahedron vertices map to
    adjacent Petersen vertices. Equivalently, no edge of the dodecahedron
    connects two vertices in the same antipodal pair. -/
theorem dodecCovering_adj :
    ∀ u v : Fin 20, dodecGraph.Adj u v →
      petGraph.Adj (dodecCoveringMap u) (dodecCoveringMap v) := by
  native_decide

/-- φ is surjective: every Petersen vertex has a preimage. -/
theorem dodecCovering_surjective :
    ∀ w : Fin 10, ∃ v : Fin 20, dodecCoveringMap v = w := by
  native_decide

/-- φ is 2-to-1: each Petersen vertex has exactly 2 preimages
    (an antipodal pair). -/
theorem dodecCovering_fiber_card :
    ∀ w : Fin 10,
      (Finset.univ.filter fun v => dodecCoveringMap v = w).card = 2 := by
  native_decide

/-! ## The S₃ generator: an A₅ element swapping 0↔6

The stabilizer Stab(0) ≅ C₃ ≤ A₅ fixes vertex 0. The setwise
stabilizer Stab_{A₅}({0,6}) ≅ S₃ contains C₃ and additionally
the element σ₃ which swaps 0↔6 (and acts on other vertices in a way
compatible with the A₅ structure).

Note: σ₃ is NOT the global antipodal map — it swaps 0↔6 but permutes
other vertices differently. The global antipodal map lives in
Aut(Dodec) \ A₅. -/

private def sigma3Fwd : Array (Fin 20) :=
  #[6, 5, 7, 4, 3, 1, 0, 2, 11, 10, 9, 8, 18, 17, 19, 16, 15, 13, 12, 14]
private def sigma3Inv : Array (Fin 20) := sigma3Fwd  -- σ₃ is an involution

private theorem sigma3Fwd_size : sigma3Fwd.size = 20 := by native_decide
private theorem sigma3Inv_size : sigma3Inv.size = 20 := by native_decide

/-- σ₃: an involution in A₅ that swaps vertices 0 and 6.
    Together with Stab(0) ≅ C₃, it generates K = Stab_{A₅}({0,6}) ≅ S₃. -/
def dodecSigma3 : Equiv.Perm (Fin 20) where
  toFun i := sigma3Fwd[i.val]'(by have := sigma3Fwd_size; omega)
  invFun i := sigma3Inv[i.val]'(by have := sigma3Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

/-- σ₃ is an involution (order 2). -/
theorem dodecSigma3_sq : dodecSigma3 * dodecSigma3 = 1 := by
  ext i; revert i; native_decide

/-- σ₃ swaps vertices 0 and 6. -/
theorem dodecSigma3_zero : dodecSigma3 0 = 6 := by native_decide
theorem dodecSigma3_six : dodecSigma3 6 = 0 := by native_decide

/-- σ₃ preserves dodecahedron adjacency. -/
theorem dodecSigma3_adj :
    ∀ u v : Fin 20, dodecGraph.Adj u v →
      dodecGraph.Adj (dodecSigma3 u) (dodecSigma3 v) := by
  native_decide

/-- σ₃ is in the generated group A₅: it equals the word g₁ · g₂⁻¹ · g₁⁻¹ · g₂⁻¹. -/
theorem dodecSigma3_eq_word :
    ∀ i : Fin 20, dodecSigma3 i = applyWordGen dodecGens [0, 3, 2, 3] i := by
  native_decide

/-- σ₃ equals the word permutation as an `Equiv.Perm`. -/
private theorem dodecSigma3_eq_applyWordGen :
    dodecSigma3 = applyWordGen dodecGens [0, 3, 2, 3] :=
  Equiv.ext (dodecSigma3_eq_word ·)

/-- **σ₃ is in the Sabidussi group** A₅ ≤ Perm(Fin 20).

    This uses `applyWordGen_mem_closure` from the Sabidussi framework:
    any word in the generators lies in `Subgroup.closure (Set.range gens)`.
    Combined with the word equality, this witnesses σ₃ as an algebraic
    element of A₅, not just an ad hoc permutation. -/
theorem dodecSigma3_mem_group :
    dodecSigma3 ∈ dodecSabidussiWitness.group := by
  rw [dodecSigma3_eq_applyWordGen]
  exact applyWordGen_mem_closure dodecGens [0, 3, 2, 3]

/-- **σ₃ is NOT in the point stabilizer** Stab(0) ≅ C₃.

    This is the key algebraic fact: since σ₃ ∈ A₅ but σ₃ ∉ Stab(0),
    the subgroup K = ⟨Stab(0), σ₃⟩ is strictly larger than Stab(0).
    The index [A₅ : K] = 10 gives exactly 10 cosets — the Petersen vertices.
    The inclusion Stab(0) ≤ K is the algebraic content of the covering
    Sab(A₅, C₃, D) → Sab(A₅, S₃, D') via `proj_adj_of_subset`. -/
theorem dodecSigma3_not_in_stabilizer :
    (⟨dodecSigma3, dodecSigma3_mem_group⟩ : ↥dodecSabidussiWitness.group) ∉
      MulAction.stabilizer dodecSabidussiWitness.group (0 : Fin 20) := by
  simp only [MulAction.mem_stabilizer_iff]
  show dodecSigma3 0 ≠ 0
  native_decide

/-- The dodecahedron is a Sabidussi coset graph: `dodecGraph ≃g Sab(A₅, C₃, D)`.
    The covering structure derives from C₃ ≤ S₃ in the coset lattice of A₅. -/
noncomputable def dodecCosetGraphIso := dodecSabidussiWitness.toSabidussiIso

/-! ## Equivariance: the A₅ action descends to the Petersen graph

The covering map is equivariant with respect to the A₅ action: each
generator of A₅ on the dodecahedron induces a well-defined permutation
of the Petersen vertices. -/

/-- Induced action of dodecGen1 (order 5) on Petersen vertices. -/
private def inducedGen1Data : Array (Fin 10) := #[8, 6, 0, 4, 2, 7, 9, 1, 3, 5]
private theorem inducedGen1Data_size : inducedGen1Data.size = 10 := by native_decide

def inducedPetGen1 (v : Fin 10) : Fin 10 :=
  inducedGen1Data[v.val]'(by have := inducedGen1Data_size; omega)

/-- Induced action of dodecGen2 (order 3) on Petersen vertices. -/
private def inducedGen2Data : Array (Fin 10) := #[6, 5, 0, 3, 1, 4, 2, 9, 7, 8]
private theorem inducedGen2Data_size : inducedGen2Data.size = 10 := by native_decide

def inducedPetGen2 (v : Fin 10) : Fin 10 :=
  inducedGen2Data[v.val]'(by have := inducedGen2Data_size; omega)

/-- φ is equivariant for generator 1: φ(g₁(v)) = g₁'(φ(v)). -/
theorem dodecCovering_equivariant_gen1 :
    ∀ v : Fin 20,
      dodecCoveringMap (dodecGen1 v) = inducedPetGen1 (dodecCoveringMap v) := by
  native_decide

/-- φ is equivariant for generator 2: φ(g₂(v)) = g₂'(φ(v)). -/
theorem dodecCovering_equivariant_gen2 :
    ∀ v : Fin 20,
      dodecCoveringMap (dodecGen2 v) = inducedPetGen2 (dodecCoveringMap v) := by
  native_decide

/-! ## Regularity -/

/-- The dodecahedron is 3-regular. -/
theorem dodec_regular :
    ∀ v : Fin 20,
      (Finset.univ.filter fun w => dodecGraph.Adj v w).card = 3 := by
  native_decide

/-- The Petersen graph is 3-regular. -/
theorem petersen_regular :
    ∀ v : Fin 10,
      (Finset.univ.filter fun w => petGraph.Adj v w).card = 3 := by
  native_decide

/-- The dodecahedron has 30 edges. -/
theorem dodec_edge_count :
    (Finset.univ.filter fun p : Fin 20 × Fin 20 =>
      p.1 < p.2 ∧ dodecGraph.Adj p.1 p.2).card = 30 := by
  native_decide

/-- The Petersen graph has 15 edges. -/
theorem petersen_edge_count :
    (Finset.univ.filter fun p : Fin 10 × Fin 10 =>
      p.1 < p.2 ∧ petGraph.Adj p.1 p.2).card = 15 := by
  native_decide
