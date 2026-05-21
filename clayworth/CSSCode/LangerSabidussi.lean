/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Langer Graph is a Sabidussi Coset Graph

The **Langer graph** (63 vertices, 6-regular) is `Sab(G₂(2), H₁₉₂, D)`.
G₂(2) = ⟨σ₁, σ₂⟩ acts vertex-transitively on the 63 points of GH(2,2),
with point stabilizer H₁₉₂ of order 192 (maximal in G₂(2)).

Uses generators and transitivity witnesses from `G2Action.lean`.
-/

import CSSCode.SabidussiRepresentation
import CSSCode.G2Action

/-! ## The G₂(2) subgroup -/

/-- The two generators as a `Fin 2`-indexed family. -/
def langerGens : Fin 2 → Equiv.Perm (Fin 63)
  | 0 => g2gen1
  | 1 => g2gen2

/-- G₂(2) as a subgroup of Sym(63). -/
def langerG : Subgroup (Equiv.Perm (Fin 63)) :=
  Subgroup.closure (Set.range langerGens)

/-- Each `applyGen` value (from G2Action) is in G₂(2). -/
private theorem applyGen_mem (i : Fin 4) : applyGen i ∈ langerG := by
  unfold langerG
  match i with
  | 0 => exact Subgroup.subset_closure ⟨0, rfl⟩
  | 1 => exact Subgroup.inv_mem _ (Subgroup.subset_closure ⟨0, rfl⟩)
  | 2 => exact Subgroup.subset_closure ⟨1, rfl⟩
  | 3 => exact Subgroup.inv_mem _ (Subgroup.subset_closure ⟨1, rfl⟩)

/-- Every `applyWord` product is in G₂(2). -/
private theorem applyWord_mem (w : List (Fin 4)) : applyWord w ∈ langerG := by
  induction w with
  | nil => exact Subgroup.one_mem _
  | cons i rest ih =>
    -- applyWord (i :: rest) = (applyGen i).trans (applyWord rest) = (applyWord rest) * (applyGen i)
    exact langerG.mul_mem ih (applyGen_mem i)

/-! ## GraphAction and transitivity -/

noncomputable instance langerMulAction : MulAction langerG (Fin 63) :=
  MulAction.compHom (Fin 63) langerG.subtype

noncomputable instance langerGraphAction :
    GraphAction langerG (Fin 63) langerSimpleGraph where
  adj_smul := by
    intro ⟨σ, hσ⟩ u v hadj
    show langerSimpleGraph.Adj (σ u) (σ v)
    revert u v; change ∀ u v, langerSimpleGraph.Adj u v → langerSimpleGraph.Adj (σ u) (σ v)
    refine Subgroup.closure_induction
      (p := fun σ _ => ∀ u v, langerSimpleGraph.Adj u v → langerSimpleGraph.Adj (σ u) (σ v))
      ?_ ?_ ?_ ?_ hσ
    · intro x ⟨i, hi⟩; subst hi
      match i with
      | 0 => exact fun u v h => (g2gen1_langer_inv u v).mp h
      | 1 => exact fun u v h => (g2gen2_langer_inv u v).mp h
    · intro u v h; simpa
    · intro x y _ _ hx hy u v hadj
      simp only [Equiv.Perm.mul_apply]; exact hx _ _ (hy u v hadj)
    · intro x _ hx u v hadj
      let f : { p : Fin 63 × Fin 63 // langerSimpleGraph.Adj p.1 p.2 } →
              { p : Fin 63 × Fin 63 // langerSimpleGraph.Adj p.1 p.2 } :=
        fun ⟨⟨a, b⟩, hab⟩ => ⟨⟨x a, x b⟩, hx a b hab⟩
      have hf_surj : Function.Surjective f := Finite.surjective_of_injective (by
        intro ⟨⟨a₁, b₁⟩, _⟩ ⟨⟨a₂, b₂⟩, _⟩ h
        simp only [f, Subtype.mk.injEq, Prod.mk.injEq] at h
        exact Subtype.ext (Prod.ext (x.injective h.1) (x.injective h.2)))
      obtain ⟨⟨⟨a, b⟩, hab⟩, heq⟩ := hf_surj ⟨⟨u, v⟩, hadj⟩
      simp only [f, Subtype.mk.injEq, Prod.mk.injEq] at heq
      rw [show a = x⁻¹ u from by rw [← heq.1]; simp,
          show b = x⁻¹ v from by rw [← heq.2]; simp] at hab
      exact hab

noncomputable instance langerPretransitive :
    MulAction.IsPretransitive langerG (Fin 63) where
  exists_smul_eq := by
    intro x y
    let σ_x := applyWord (witnessWord x)
    let σ_y := applyWord (witnessWord y)
    have hmem : σ_x.symm.trans σ_y ∈ langerG :=
      langerG.mul_mem (applyWord_mem _) (langerG.inv_mem (applyWord_mem _))
    exact ⟨⟨σ_x.symm.trans σ_y, hmem⟩, by
      show (σ_x.symm.trans σ_y) x = y
      simp only [Equiv.trans_apply]
      rw [show σ_x.symm x = 0 from by
        rw [Equiv.symm_apply_eq]; exact (witnessWord_correct x).symm]
      exact witnessWord_correct y⟩

/-! ## The main theorem -/

/-- **The Langer graph is a Sabidussi coset graph.**

G₂(2) = ⟨σ₁, σ₂⟩ acts vertex-transitively on the 63 points of GH(2,2).
The stabilizer of vertex 0 is H₁₉₂ (order 192), giving:
  Langer ≃g Sab(G₂(2), H₁₉₂, D)

H₁₉₂ is maximal in G₂(2) (proved in `langer_action_primitive`),
so the Langer graph admits no non-trivial symmetric quotient. -/
noncomputable def langerSabidussiIso :
    langerSimpleGraph ≃g SimpleGraph.cosetGraph
      (MulAction.stabilizer langerG (0 : Fin 63))
      (connectionSet langerG langerSimpleGraph 0)
      (connectionSet.isConnectionSet 0) :=
  sabidussiIso 0
