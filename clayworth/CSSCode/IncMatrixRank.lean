/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Rank of the Boundary Operators

Upper and lower bounds on `rank(d1)` and `rank(d2)` over F₂ = ZMod 2.

**Upper bound** (proved): the transpose d1ᵀ has a nontrivial kernel (the all-ones
vector, since 2 = 0 in F₂). By rank-nullity, `rank(d1) ≤ V - 1`.

**Lower bound** (sorry'd): a spanning tree gives V - 1 linearly independent columns.
-/

import CSSCode.CellularSurface
import Mathlib.Algebra.Field.ZMod
import Mathlib.LinearAlgebra.Dimension.Finite
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.Combinatorics.SimpleGraph.Acyclic

open Matrix Finset

-- ZMod 2 is a field (2 is prime)
attribute [local instance] ZMod.instField

namespace CellularSurface

variable (S : CellularSurface)

/-! ### Upper bound: rank(d1) ≤ V - 1 -/

/-- Each column of d1 sums to 0 over F₂: every edge has exactly 2 endpoints. -/
theorem d1_col_sum_eq_zero (e : Fin S.nE) :
    ∑ v : Fin S.nV, S.d1 v e = 0 := by
  simp only [d1, Finset.sum_add_distrib]
  -- Each Kronecker sum equals 1 (exactly one vertex matches)
  have hsum : ∀ (a : Fin S.nV),
      ∑ v : Fin S.nV, (if a = v then (1 : ZMod 2) else 0) = 1 := fun a => by
    rw [Finset.sum_eq_single_of_mem a (Finset.mem_univ _)
      (fun b _ hne => if_neg (Ne.symm hne)), if_pos rfl]
  simp only [hsum]
  -- 1 + 1 = 0 in ZMod 2
  decide

/-- The all-ones vector lies in the kernel of d1ᵀ. -/
theorem one_mem_ker_d1T :
    (1 : Fin S.nV → ZMod 2) ∈ LinearMap.ker S.d1ᵀ.mulVecLin := by
  rw [LinearMap.mem_ker]
  ext e
  simp only [mulVecLin_apply, mulVec, dotProduct, transpose_apply, Pi.zero_apply,
    Pi.one_apply, mul_one]
  exact S.d1_col_sum_eq_zero e

/-- The all-ones vector in F₂^V is nonzero when V > 0. -/
theorem one_ne_zero_fin (hV : 0 < S.nV) :
    (1 : Fin S.nV → ZMod 2) ≠ 0 := by
  intro h; exact absurd (congr_fun h ⟨0, hV⟩) (by simp)

/-- **Upper bound**: `rank(d1) ≤ nV - 1`. -/
theorem d1_rank_le (hV : 0 < S.nV) : S.d1.rank ≤ S.nV - 1 := by
  rw [← rank_transpose S.d1]
  -- rank-nullity: rank(d1ᵀ) + dim(ker) = nV
  have hrn := LinearMap.finrank_range_add_finrank_ker S.d1ᵀ.mulVecLin
  simp only [rank, Module.finrank_pi_fintype, Module.finrank_self, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, smul_eq_mul, mul_one] at hrn ⊢
  -- ker(d1ᵀ) ≠ ⊥ since it contains the nonzero all-ones vector
  have hker_ne : LinearMap.ker S.d1ᵀ.mulVecLin ≠ ⊥ := by
    rw [ne_eq, LinearMap.ker_eq_bot']
    push_neg
    exact ⟨1, (LinearMap.mem_ker.mp S.one_mem_ker_d1T), S.one_ne_zero_fin hV⟩
  -- ker ≠ ⊥ implies dim(ker) ≥ 1
  have hker_pos : 1 ≤ Module.finrank (ZMod 2) (LinearMap.ker S.d1ᵀ.mulVecLin) :=
    Submodule.one_le_finrank_iff.mpr hker_ne
  omega

/-! ### Upper bound: rank(d2) ≤ F - 1 -/

/-- **Upper bound**: `rank(d2) ≤ nF - 1`, given that each edge borders an even
    number of face sides (`d2 *ᵥ 1 = 0`). -/
theorem d2_rank_le (hF : 0 < S.nF) (htwo : S.d2 *ᵥ 1 = 0) :
    S.d2.rank ≤ S.nF - 1 := by
  unfold rank
  have hrn := LinearMap.finrank_range_add_finrank_ker S.d2.mulVecLin
  simp only [Module.finrank_pi_fintype, Module.finrank_self, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, smul_eq_mul, mul_one] at hrn ⊢
  have hker_ne : LinearMap.ker S.d2.mulVecLin ≠ ⊥ := by
    rw [ne_eq, LinearMap.ker_eq_bot']
    push_neg
    refine ⟨1, ?_, ?_⟩
    · show S.d2.mulVecLin 1 = 0
      rwa [mulVecLin_apply]
    · intro h; exact absurd (congr_fun h ⟨0, hF⟩) (by simp)
  have hker_pos : 1 ≤ Module.finrank (ZMod 2) (LinearMap.ker S.d2.mulVecLin) :=
    Submodule.one_le_finrank_iff.mpr hker_ne
  omega

/-! ### Bridge to SimpleGraph -/

/-- The underlying simple graph of a CellularSurface. -/
def toSimpleGraph : SimpleGraph (Fin S.nV) where
  Adj u v := ∃ e : Fin S.nE,
    (S.edge_src e = u ∧ S.edge_tgt e = v) ∨ (S.edge_src e = v ∧ S.edge_tgt e = u)
  symm u v := fun ⟨e, h⟩ => ⟨e, h.symm⟩
  loopless := ⟨fun v ⟨e, h⟩ => by
    rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> exact S.edge_ne e (h1.trans h2.symm)⟩

/-- The dual graph: faces as vertices, adjacent when sharing a boundary edge. -/
def toDualSimpleGraph : SimpleGraph (Fin S.nF) where
  Adj f1 f2 := f1 ≠ f2 ∧ ∃ e : Fin S.nE,
    (∃ i, S.face_edge f1 i = e) ∧ (∃ j, S.face_edge f2 j = e)
  symm _ _ := fun ⟨hne, e, h1, h2⟩ => ⟨hne.symm, e, h2, h1⟩
  loopless := ⟨fun _ ⟨hne, _⟩ => hne rfl⟩

/-! ### Lower bound: rank(d1) ≥ V - 1 via kernel characterisation

The key idea: `ker(d1ᵀ) = span{1_V}` for connected graphs. Every element
of the kernel assigns equal values to adjacent vertices (since `d1ᵀ *ᵥ x = 0`
means `x(src e) + x(tgt e) = 0` for all edges, i.e., `x(src e) = x(tgt e)`).
By connectivity, this propagates along walks to give `x` constant on all
vertices. Over F₂, constant means `x = 0` or `x = 1_V`, so `ker = span{1_V}`
has dimension 1. Rank-nullity then gives `rank = V - 1`.
-/

/-- Each entry of `d1ᵀ *ᵥ x` equals `x(src e) + x(tgt e)`. -/
theorem d1T_mulVec_entry (x : Fin S.nV → ZMod 2) (e : Fin S.nE) :
    (S.d1ᵀ *ᵥ x) e = x (S.edge_src e) + x (S.edge_tgt e) := by
  simp only [mulVec, dotProduct, transpose_apply, d1, add_mul, Finset.sum_add_distrib]
  have hsum : ∀ (a : Fin S.nV),
      ∑ v : Fin S.nV, (if a = v then (1 : ZMod 2) else 0) * x v = x a := fun a => by
    rw [Finset.sum_eq_single_of_mem a (Finset.mem_univ _)
      (fun b _ hne => by simp [Ne.symm hne])]
    simp
  simp only [hsum]

/-- Elements of `ker(d1ᵀ)` assign equal values to adjacent vertices. -/
theorem ker_d1T_edge_eq {x : Fin S.nV → ZMod 2}
    (hx : x ∈ LinearMap.ker S.d1ᵀ.mulVecLin) (e : Fin S.nE) :
    x (S.edge_src e) = x (S.edge_tgt e) := by
  have h := congr_fun (LinearMap.mem_ker.mp hx) e
  simp only [mulVecLin_apply, Pi.zero_apply] at h
  rw [S.d1T_mulVec_entry x e] at h
  have : ∀ (a b : ZMod 2), a + b = 0 → a = b := by decide
  exact this _ _ h

/-- Walks in the graph preserve values of kernel elements. -/
theorem walk_preserves_ker {x : Fin S.nV → ZMod 2}
    (hx : x ∈ LinearMap.ker S.d1ᵀ.mulVecLin)
    {u v : Fin S.nV} (w : S.toSimpleGraph.Walk u v) :
    x u = x v := by
  induction w with
  | nil => rfl
  | cons hadj _ ih =>
    obtain ⟨e, h | h⟩ := hadj
    · -- src e = u, tgt e = w: x u = x(src e) = x(tgt e) = x w
      have := S.ker_d1T_edge_eq hx e
      rw [h.1, h.2] at this
      exact this.trans ih
    · -- src e = w, tgt e = u: x u = x(tgt e) = x(src e) = x w
      have := S.ker_d1T_edge_eq hx e
      rw [h.1, h.2] at this
      exact this.symm.trans ih

/-- For connected graphs, `ker(d1ᵀ) ≤ span{1_V}`: kernel elements are constant. -/
theorem ker_d1T_le_span_one (hconn : S.toSimpleGraph.Connected) :
    LinearMap.ker S.d1ᵀ.mulVecLin ≤ Submodule.span (ZMod 2) {1} := by
  intro x hx
  rw [Submodule.mem_span_singleton]
  have hV : 0 < S.nV := by obtain ⟨⟨_, hlt⟩⟩ := hconn.nonempty; omega
  -- x is constant on all vertices (by connectivity + walk preservation)
  have hconst : ∀ u v : Fin S.nV, x u = x v := fun u v =>
    S.walk_preserves_ker hx (hconn.preconnected u v).some
  -- x = x(0) • 1, so x ∈ span{1}
  exact ⟨x ⟨0, hV⟩, funext fun v => by
    simp [Pi.smul_apply, Pi.one_apply, smul_eq_mul, mul_one, hconst ⟨0, hV⟩ v]⟩

/-- For connected graphs, `ker(d1ᵀ)` has dimension exactly 1. -/
theorem ker_d1T_finrank_eq_one (hconn : S.toSimpleGraph.Connected) (hV : 0 < S.nV) :
    Module.finrank (ZMod 2) (LinearMap.ker S.d1ᵀ.mulVecLin) = 1 := by
  -- ker = span{1_V} by antisymmetry
  have heq : LinearMap.ker S.d1ᵀ.mulVecLin = Submodule.span (ZMod 2) {(1 : Fin S.nV → ZMod 2)} :=
    le_antisymm (S.ker_d1T_le_span_one hconn)
      ((Submodule.span_singleton_le_iff_mem _ _).mpr S.one_mem_ker_d1T)
  rw [heq]
  exact finrank_span_singleton (S.one_ne_zero_fin hV)

/-- **Lower bound**: `rank(d1) ≥ nV - 1` for connected graphs. -/
theorem d1_rank_ge (hconn : S.toSimpleGraph.Connected) (hV : 1 < S.nV) :
    S.nV - 1 ≤ S.d1.rank := by
  rw [← rank_transpose S.d1]
  have hrn := LinearMap.finrank_range_add_finrank_ker S.d1ᵀ.mulVecLin
  simp only [rank, Module.finrank_pi_fintype, Module.finrank_self, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, smul_eq_mul, mul_one] at hrn ⊢
  rw [S.ker_d1T_finrank_eq_one hconn (by omega)] at hrn
  omega

/-- **Rank theorem for d1**: `rank(d1) = nV - 1` for connected graphs. -/
theorem d1_rank_eq (hconn : S.toSimpleGraph.Connected) (hV : 1 < S.nV) :
    (S.d1.rank : ℤ) = S.nV - 1 := by
  have hle := S.d1_rank_le (by omega : 0 < S.nV)
  have hge := S.d1_rank_ge hconn hV
  omega

/-- **Lower bound**: `rank(d2) ≥ nF - 1` for connected dual graphs. -/
theorem d2_rank_ge (hconn_dual : S.toDualSimpleGraph.Connected) (hF : 1 < S.nF)
    (htwo : S.d2 *ᵥ 1 = 0) :
    S.nF - 1 ≤ S.d2.rank := by
  sorry

/-- **Rank theorem for d2**: `rank(d2) = nF - 1`. -/
theorem d2_rank_eq (hconn_dual : S.toDualSimpleGraph.Connected) (hF : 1 < S.nF)
    (htwo : S.d2 *ᵥ 1 = 0) :
    (S.d2.rank : ℤ) = S.nF - 1 := by
  have hle := S.d2_rank_le (by omega : 0 < S.nF) htwo
  have hge := S.d2_rank_ge hconn_dual hF htwo
  omega

end CellularSurface
