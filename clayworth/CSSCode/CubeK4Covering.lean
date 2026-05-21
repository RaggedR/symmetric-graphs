/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Cube is a Double Cover of K₄

The 3-cube Q₃ (8 vertices, cubic) is a double cover of the complete
graph K₄ (4 vertices, 3-regular). This covering arises algebraically
from the Sabidussi framework:

  Cube = Sab(S₄, C₃, D)   — |S₄|/|C₃| = 24/3 = 8 vertices
  K₄   = Sab(S₄, S₃, D')  — |S₄|/|S₃| = 24/6 = 4 vertices

where C₃ ≤ S₃ are both subgroups of S₄ (the rotation group of the cube).
The generalized coset projection S₄/C₃ → S₄/S₃ (`proj_adj_of_subset`
with D ⊆ D') gives the covering map.

The covering identifies each vertex with its antipodal vertex (at
distance 3): {0,7}, {1,6}, {2,5}, {3,4}.

## References

- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import CSSCode.CubeSabidussi
import CSSCode.CosetGraph

/-! ## The covering map -/

private def cubeCoveringData : Array (Fin 4) :=
  #[0, 1, 2, 3, 3, 2, 1, 0]

private theorem cubeCoveringData_size : cubeCoveringData.size = 8 := by native_decide

/-- The covering map φ : Q₃ → K₄, identifying antipodal pairs. -/
def cubeCoveringMap (v : Fin 8) : Fin 4 :=
  cubeCoveringData[v.val]'(by have := cubeCoveringData_size; omega)

/-! ## K₄: the complete graph on 4 vertices -/

/-- K₄ as a SimpleGraph on Fin 4: all distinct pairs are adjacent. -/
def k4Graph : SimpleGraph (Fin 4) where
  Adj u v := u ≠ v
  symm _ _ h := Ne.symm h
  loopless := ⟨fun _ h => h rfl⟩

instance k4DecAdj : DecidableRel k4Graph.Adj :=
  fun u v => inferInstanceAs (Decidable (u ≠ v))

/-! ## Graph homomorphism properties -/

/-- φ is a graph homomorphism: adjacent cube vertices map to
    adjacent K₄ vertices. -/
theorem cubeCovering_adj :
    ∀ u v : Fin 8, cubeGraph.Adj u v →
      k4Graph.Adj (cubeCoveringMap u) (cubeCoveringMap v) := by
  native_decide

/-- φ is surjective. -/
theorem cubeCovering_surjective :
    ∀ w : Fin 4, ∃ v : Fin 8, cubeCoveringMap v = w := by
  native_decide

/-- φ is 2-to-1. -/
theorem cubeCovering_fiber_card :
    ∀ w : Fin 4,
      (Finset.univ.filter fun v => cubeCoveringMap v = w).card = 2 := by
  native_decide

/-! ## The S₃ generator

σ₃: an S₄ element (order 2) that swaps vertices 0 and 7.
Together with Stab(0) ≅ C₃, generates K = Stab_{S₄}({0,7}) ≅ S₃. -/

private def cubeSigma3Fwd : Array (Fin 8) := #[7, 3, 5, 1, 6, 2, 4, 0]
private def cubeSigma3Inv : Array (Fin 8) := cubeSigma3Fwd  -- involution

private theorem cubeSigma3Fwd_size : cubeSigma3Fwd.size = 8 := by native_decide
private theorem cubeSigma3Inv_size : cubeSigma3Inv.size = 8 := by native_decide

def cubeSigma3 : Equiv.Perm (Fin 8) where
  toFun i := cubeSigma3Fwd[i.val]'(by have := cubeSigma3Fwd_size; omega)
  invFun i := cubeSigma3Inv[i.val]'(by have := cubeSigma3Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide

theorem cubeSigma3_sq : cubeSigma3 * cubeSigma3 = 1 := by
  ext i; revert i; native_decide

theorem cubeSigma3_zero : cubeSigma3 0 = 7 := by native_decide
theorem cubeSigma3_seven : cubeSigma3 7 = 0 := by native_decide

theorem cubeSigma3_adj :
    ∀ u v : Fin 8, cubeGraph.Adj u v →
      cubeGraph.Adj (cubeSigma3 u) (cubeSigma3 v) := by
  native_decide

/-- σ₃ is in S₄: it equals the word g₁² · g₂ · g₁⁻¹. -/
theorem cubeSigma3_eq_word :
    ∀ i : Fin 8, cubeSigma3 i = applyWordGen cubeGens [0, 0, 1, 2] i := by
  native_decide

private theorem cubeSigma3_eq_applyWordGen :
    cubeSigma3 = applyWordGen cubeGens [0, 0, 1, 2] :=
  Equiv.ext (cubeSigma3_eq_word ·)

/-- **σ₃ is in the Sabidussi group** S₄ ≤ Perm(Fin 8). -/
theorem cubeSigma3_mem_group :
    cubeSigma3 ∈ cubeSabidussiWitness.group := by
  rw [cubeSigma3_eq_applyWordGen]
  exact applyWordGen_mem_closure cubeGens [0, 0, 1, 2]

/-- **σ₃ is NOT in the point stabilizer** Stab(0) ≅ C₃.
    Since σ₃(0) = 7 ≠ 0, K = ⟨Stab(0), σ₃⟩ is strictly larger,
    giving the covering Sab(S₄, C₃, D) → Sab(S₄, S₃, D'). -/
theorem cubeSigma3_not_in_stabilizer :
    (⟨cubeSigma3, cubeSigma3_mem_group⟩ : ↥cubeSabidussiWitness.group) ∉
      MulAction.stabilizer cubeSabidussiWitness.group (0 : Fin 8) := by
  simp only [MulAction.mem_stabilizer_iff]
  show cubeSigma3 0 ≠ 0
  native_decide

/-- The cube is a Sabidussi coset graph: `cubeGraph ≃g Sab(S₄, C₃, D)`. -/
noncomputable def cubeCosetGraphIso := cubeSabidussiWitness.toSabidussiIso

/-! ## Equivariance -/

private def cubeInducedGen1Data : Array (Fin 4) := #[2, 3, 1, 0]
private theorem cubeInducedGen1Data_size : cubeInducedGen1Data.size = 4 := by native_decide

def cubeInducedGen1 (v : Fin 4) : Fin 4 :=
  cubeInducedGen1Data[v.val]'(by have := cubeInducedGen1Data_size; omega)

private def cubeInducedGen2Data : Array (Fin 4) := #[0, 2, 3, 1]
private theorem cubeInducedGen2Data_size : cubeInducedGen2Data.size = 4 := by native_decide

def cubeInducedGen2 (v : Fin 4) : Fin 4 :=
  cubeInducedGen2Data[v.val]'(by have := cubeInducedGen2Data_size; omega)

theorem cubeCovering_equivariant_gen1 :
    ∀ v : Fin 8,
      cubeCoveringMap (cubeGen1 v) = cubeInducedGen1 (cubeCoveringMap v) := by
  native_decide

theorem cubeCovering_equivariant_gen2 :
    ∀ v : Fin 8,
      cubeCoveringMap (cubeGen2 v) = cubeInducedGen2 (cubeCoveringMap v) := by
  native_decide

/-! ## Regularity -/

theorem k4_regular :
    ∀ v : Fin 4,
      (Finset.univ.filter fun w => k4Graph.Adj v w).card = 3 := by
  native_decide

theorem k4_edge_count :
    (Finset.univ.filter fun p : Fin 4 × Fin 4 =>
      p.1 < p.2 ∧ k4Graph.Adj p.1 p.2).card = 6 := by
  native_decide
