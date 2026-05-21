/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Coset Graphs (Sabidussi Construction)

Local port of the `CosetGraph` construction from mathlib4 PR #39548.
Defines `IsConnectionSet` and `SimpleGraph.cosetGraph`.

Given a group G, a subgroup H ≤ G, and a connection set D ⊆ G, the
**coset graph** Sab(G, H, D) has vertices G/H and adjacency Hx ~ Hy iff x⁻¹y ∈ D.

## References

- Sabidussi, *Vertex-transitive graphs*, Monatsh. Math. 68 (1964)
- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.GroupTheory.QuotientGroup.Basic
import Mathlib.Tactic.Group

/-! ## Connection sets -/

/-- A connection set for H ≤ G: a subset D ⊆ G that is a union of (H,H)-double
    cosets, closed under inversion, and disjoint from H. -/
structure IsConnectionSet {G : Type*} [Group G] (H : Subgroup G) (D : Set G) : Prop where
  double_coset_stable : ∀ d ∈ D, ∀ h₁ ∈ (H : Set G), ∀ h₂ ∈ (H : Set G), h₁ * d * h₂ ∈ D
  inv_mem : ∀ d ∈ D, d⁻¹ ∈ D
  not_mem_of_mem_subgroup : ∀ d ∈ D, d ∉ (H : Set G)

namespace IsConnectionSet

variable {G : Type*} [Group G] {H : Subgroup G} {D : Set G}

theorem one_not_mem (hD : IsConnectionSet H D) : (1 : G) ∉ D :=
  fun h => hD.not_mem_of_mem_subgroup 1 h H.one_mem

/-- The coset relation is well-defined: if x₁ ≡ x₂ mod H and y₁ ≡ y₂ mod H,
    then x₁⁻¹y₁ ∈ D ↔ x₂⁻¹y₂ ∈ D. -/
theorem mem_iff_of_leftRel (hD : IsConnectionSet H D)
    {x₁ y₁ x₂ y₂ : G}
    (hx : (QuotientGroup.leftRel H).r x₁ x₂)
    (hy : (QuotientGroup.leftRel H).r y₁ y₂) :
    x₁⁻¹ * y₁ ∈ D ↔ x₂⁻¹ * y₂ ∈ D := by
  rw [QuotientGroup.leftRel_apply] at hx hy
  constructor
  · intro hadj
    -- x₂⁻¹y₂ = (x₂⁻¹x₁)(x₁⁻¹y₁)(y₁⁻¹y₂)
    have h1 : x₂⁻¹ * x₁ ∈ (H : Set G) := by
      show x₂⁻¹ * x₁ ∈ H
      have hmem : x₁⁻¹ * x₂ ∈ H := hx
      have hinv := H.inv_mem hmem
      simp only [mul_inv_rev, inv_inv] at hinv
      exact hinv
    have h2 : y₁⁻¹ * y₂ ∈ (H : Set G) := hy
    have key : x₂⁻¹ * y₂ = (x₂⁻¹ * x₁) * (x₁⁻¹ * y₁) * (y₁⁻¹ * y₂) := by group
    rw [key]
    exact hD.double_coset_stable _ hadj _ h1 _ h2
  · intro hadj
    have h1 : x₁⁻¹ * x₂ ∈ (H : Set G) := hx
    have h2 : y₂⁻¹ * y₁ ∈ (H : Set G) := by
      show y₂⁻¹ * y₁ ∈ H
      have hmem : y₁⁻¹ * y₂ ∈ H := hy
      have hinv := H.inv_mem hmem
      simp only [mul_inv_rev, inv_inv] at hinv
      exact hinv
    have key : x₁⁻¹ * y₁ = (x₁⁻¹ * x₂) * (x₂⁻¹ * y₂) * (y₂⁻¹ * y₁) := by group
    rw [key]
    exact hD.double_coset_stable _ hadj _ h1 _ h2

end IsConnectionSet

/-! ## The coset graph -/

/-- The **coset graph** Sab(G, H, D): vertices are right cosets G/H,
    Hx ~ Hy iff x⁻¹y ∈ D. -/
noncomputable def SimpleGraph.cosetGraph {G : Type*} [Group G]
    (H : Subgroup G) (D : Set G) (hD : IsConnectionSet H D) :
    SimpleGraph (G ⧸ H) where
  Adj q₁ q₂ := Quotient.liftOn₂ q₁ q₂ (fun x y => x⁻¹ * y ∈ D)
    (fun _ _ _ _ hx hy => propext (hD.mem_iff_of_leftRel hx hy))
  symm := by
    intro q₁ q₂
    refine Quotient.inductionOn₂ q₁ q₂ ?_
    intro x y
    simp only [Quotient.liftOn₂_mk]
    intro h
    have : y⁻¹ * x = (x⁻¹ * y)⁻¹ := by group
    rw [this]
    exact hD.inv_mem _ h
  loopless := ⟨by
    intro q
    refine Quotient.inductionOn q ?_
    intro x
    simp only [Quotient.liftOn₂_mk]
    intro h
    have : x⁻¹ * x = 1 := by group
    rw [this] at h
    exact hD.one_not_mem h⟩

namespace SimpleGraph.cosetGraph

variable {G : Type*} [Group G] {H : Subgroup G} {D : Set G} (hD : IsConnectionSet H D)

/-- Adjacency characterisation: ⟦x⟧ ~ ⟦y⟧ iff x⁻¹y ∈ D. -/
@[simp]
theorem adj_mk (x y : G) :
    (cosetGraph H D hD).Adj ⟦x⟧ ⟦y⟧ ↔ x⁻¹ * y ∈ D :=
  Iff.rfl

end SimpleGraph.cosetGraph
