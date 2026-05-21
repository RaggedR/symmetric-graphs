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

/-! ## Coset projection: H ≤ K gives Sab(G,H,D) → Sab(G,K,D)

When H ≤ K, the natural projection G/H → G/K sends gH ↦ gK. If D is
a connection set for both H and K, then this projection is a graph
homomorphism: adjacency in Sab(G,H,D) implies adjacency in Sab(G,K,D).

This is the Sabidussi analogue of "quotient by a block system": the
block system {gK/H : g ∈ G} on G/H has blocks of size [K:H], and the
quotient graph is Sab(G,K,D).

**Application**: Zhou-3 (Sab(PSL(2,13), S₃), 182v, cubic) projects to
Zhou-6 (Sab(PSL(2,13), D₁₂), 91v, 6-regular) via S₃ ≤ D₁₂ with [D₁₂:S₃] = 2.
Each Zhou-6 vertex has exactly 2 preimages, each of degree 3, whose 6
projected edges give the 6-regular structure. -/

/-- The natural projection G/H → G/K when H ≤ K, sending ⟦g⟧_H to ⟦g⟧_K. -/
def cosetProjection {G : Type*} [Group G] (H K : Subgroup G) (hle : H ≤ K) :
    G ⧸ H → G ⧸ K :=
  Quotient.map id fun a b hab => by
    change (QuotientGroup.leftRel K).r a b
    have hab' : (QuotientGroup.leftRel H).r a b := hab
    rw [QuotientGroup.leftRel_apply] at hab' ⊢
    exact hle hab'

@[simp]
theorem cosetProjection_mk {G : Type*} [Group G] (H K : Subgroup G) (hle : H ≤ K)
    (g : G) : cosetProjection H K hle ⟦g⟧ = ⟦g⟧ :=
  rfl

/-- **Coset projection is a graph homomorphism**: if H ≤ K and D is a
    connection set for both, then adjacency in Sab(G,H,D) implies
    adjacency in Sab(G,K,D).

    The proof is trivial: adjacency is `x⁻¹y ∈ D` in both graphs.
    The only difference is the coset equivalence, which is coarser for K. -/
theorem SimpleGraph.cosetGraph.proj_adj {G : Type*} [Group G]
    {H K : Subgroup G} (hle : H ≤ K) {D : Set G}
    (hDH : IsConnectionSet H D) (hDK : IsConnectionSet K D) :
    ∀ q₁ q₂ : G ⧸ H,
      (cosetGraph H D hDH).Adj q₁ q₂ →
      (cosetGraph K D hDK).Adj (cosetProjection H K hle q₁)
                                (cosetProjection H K hle q₂) := by
  intro q₁ q₂
  refine Quotient.inductionOn₂ q₁ q₂ ?_
  intro x y
  simp only [cosetProjection_mk]
  exact fun h => h

/-- The coset projection is surjective. -/
theorem cosetProjection_surjective {G : Type*} [Group G]
    (H K : Subgroup G) (hle : H ≤ K) :
    Function.Surjective (cosetProjection H K hle) := by
  intro q
  refine Quotient.inductionOn q ?_
  intro g
  exact ⟨⟦g⟧, by simp⟩
