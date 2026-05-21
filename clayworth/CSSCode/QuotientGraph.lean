/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Quotient Graphs

The quotient of a simple graph by a surjection π : V → W. Two vertices
i, j in W are adjacent iff i ≠ j and there exist preimages u, v with
π u = i, π v = j, and u ~ v in the original graph.

Local port from mathlib4 PR #39551.

## References

- Lorimer, *Vertex-transitive graphs*, J. Graph Theory 8 (1984)
- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import Mathlib.Combinatorics.SimpleGraph.Basic

/-! ## Quotient graph construction -/

/-- The **quotient graph** of Γ by a map π : V → W.
    Two classes i, j : W are adjacent iff i ≠ j and some edge of Γ
    connects a vertex in π⁻¹(i) to a vertex in π⁻¹(j). -/
def SimpleGraph.quotientGraph {V W : Type*} (Γ : SimpleGraph V) (π : V → W)
    [DecidableEq W] : SimpleGraph W where
  Adj i j := i ≠ j ∧ ∃ u v, π u = i ∧ π v = j ∧ Γ.Adj u v
  symm := by
    intro i j ⟨hne, u, v, hu, hv, hadj⟩
    exact ⟨hne.symm, v, u, hv, hu, hadj.symm⟩
  loopless := ⟨fun i ⟨hne, _⟩ => hne rfl⟩

namespace SimpleGraph.quotientGraph

variable {V W : Type*} {Γ : SimpleGraph V} {π : V → W} [DecidableEq W]

/-- An edge in Γ projects to an edge in the quotient (when endpoints map to
    distinct vertices). -/
theorem adj_of_adj {u v : V} (hadj : Γ.Adj u v) (hne : π u ≠ π v) :
    (Γ.quotientGraph π).Adj (π u) (π v) :=
  ⟨hne, u, v, rfl, rfl, hadj⟩

end SimpleGraph.quotientGraph
