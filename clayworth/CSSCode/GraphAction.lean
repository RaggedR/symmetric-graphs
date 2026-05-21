/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Group Actions on Simple Graphs

Local port of the `GraphAction` typeclass from mathlib4 PR #39530.

A group action on `V` is a *graph action* on `Γ : SimpleGraph V` if it
preserves adjacency: `Γ.Adj u v → Γ.Adj (g • u) (g • v)` for all `g : G`.
-/

import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.GroupTheory.GroupAction.Basic

/-- A group action on `V` is a *graph action* on `Γ : SimpleGraph V` if it preserves adjacency. -/
class GraphAction (G : Type*) (V : Type*) [Group G] [MulAction G V]
    (Γ : SimpleGraph V) : Prop where
  /-- The action preserves adjacency. -/
  adj_smul : ∀ (g : G) (u v : V), Γ.Adj u v → Γ.Adj (g • u) (g • v)

namespace GraphAction

variable {G V : Type*} [Group G] [MulAction G V] {Γ : SimpleGraph V} [GraphAction G V Γ]

/-- For group actions, adjacency is preserved in both directions. -/
theorem adj_smul_iff (g : G) (u v : V) : Γ.Adj (g • u) (g • v) ↔ Γ.Adj u v := by
  constructor
  · intro h
    have h' := adj_smul g⁻¹ (g • u) (g • v) h
    simp only [inv_smul_smul] at h'
    exact h'
  · exact adj_smul g u v

end GraphAction
