/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# SimpleGraph from a triangle group, with algebraic Sabidussi proof

Given `TriangleGroupData` for Δ(2,3,r) acting on a finite group G:
  - Vertices = left cosets of ⟨S⟩
  - Adjacency: v ~ w iff ∃ dart d at v whose T-image is at w

The graph is a coset graph by construction: the connection set is the
single double coset ⟨S⟩T⟨S⟩. Sabidussi's theorem is essentially
tautological — the group acts on cosets by left multiplication, and
left multiplication commutes with right multiplication (T, S).

## Key theorem

`DartLeftMult.preserves_adj`: if a dart permutation σ commutes with T
and induces a vertex permutation σ_v via vertexOf, then σ_v preserves
adjacency. Proof: apply σ to the witnessing dart. Cost: O(nDarts)
instead of O(nV²).
-/

import CSSCode.TriangleGroupSurface
import Mathlib.Combinatorics.SimpleGraph.Basic

namespace TriangleGroupData

variable {nDarts nV nE nF r : ℕ} (D : TriangleGroupData nDarts nV nE nF r)

/-! ## SimpleGraph from triangle group data -/

/-- The **triangle group graph**: vertices are ⟨S⟩-cosets, two vertices
are adjacent iff some dart at the first has its T-image at the second.

Symmetry is algebraic: T is an involution (T² = id in Δ(2,3,r)),
so if dart d witnesses v ~ w, then T(d) witnesses w ~ v. -/
def triangleGraph : SimpleGraph (Fin nV) where
  Adj v w := v ≠ w ∧ ∃ d : Fin nDarts, D.vertexOf d = v ∧ D.vertexOf (D.T d) = w
  symm v w := by
    intro ⟨hne, d, hv, hw⟩
    refine ⟨hne.symm, D.T d, hw, ?_⟩
    change D.vertexOf ((D.R.trans D.S) ((D.R.trans D.S) d)) = v
    rw [D.T_inv d, hv]
  loopless := ⟨fun v ⟨h, _⟩ => h rfl⟩

/-! ## Dart-level left multiplication -/

/-- A **dart-level left multiplication** permutation with the two key
properties needed for the algebraic Sabidussi proof:

1. `comm_T`: commutes with the edge involution T = RS
2. `vertex_compat`: induces a vertex permutation via vertexOf

These hold for any left multiplication in the regular representation,
since left and right multiplication commute. -/
structure DartLeftMult where
  /-- The dart-level permutation (left multiplication by some g ∈ G) -/
  σ : Equiv.Perm (Fin nDarts)
  /-- Commutes with T = RS (the edge involution) -/
  comm_T : ∀ d, σ (D.T d) = D.T (σ d)
  /-- The induced vertex permutation -/
  σ_v : Equiv.Perm (Fin nV)
  /-- Vertex induction: vertexOf(σ(d)) = σ_v(vertexOf(d)) -/
  vertex_compat : ∀ d, D.vertexOf (σ d) = σ_v (D.vertexOf d)

namespace DartLeftMult

variable {D : TriangleGroupData nDarts nV nE nF r}

/-- **Algebraic adjacency preservation**: if σ commutes with T and
induces σ_v on vertices, then σ_v preserves adjacency.

Proof: given dart d witnessing v ~ w, the dart σ(d) witnesses
σ_v(v) ~ σ_v(w). The key step: vertexOf(T(σ(d))) = vertexOf(σ(T(d)))
(by comm_T) = σ_v(vertexOf(T(d))) (by vertex_compat) = σ_v(w). -/
theorem preserves_adj (L : DartLeftMult D) (v w : Fin nV)
    (hadj : (D.triangleGraph).Adj v w) :
    (D.triangleGraph).Adj (L.σ_v v) (L.σ_v w) := by
  obtain ⟨hne, d, hv, hw⟩ := hadj
  refine ⟨L.σ_v.injective.ne hne, L.σ d, ?_, ?_⟩
  · rw [L.vertex_compat, hv]
  · rw [← L.comm_T, L.vertex_compat, hw]

/-- The inverse of a dart left multiplication is also a dart left
multiplication. Algebraic: σ⁻¹ commutes with T iff σ does. -/
def inv (L : DartLeftMult D) : DartLeftMult D where
  σ := L.σ.symm
  comm_T := by
    intro d
    apply L.σ.injective
    rw [Equiv.apply_symm_apply, L.comm_T, Equiv.apply_symm_apply]
  σ_v := L.σ_v.symm
  vertex_compat := by
    intro d
    apply L.σ_v.injective
    rw [Equiv.apply_symm_apply, ← L.vertex_compat, Equiv.apply_symm_apply]

/-- Inverse preserves adjacency (follows from `preserves_adj` + `inv`). -/
theorem inv_preserves_adj (L : DartLeftMult D) (v w : Fin nV)
    (hadj : (D.triangleGraph).Adj v w) :
    (D.triangleGraph).Adj (L.σ_v.symm v) (L.σ_v.symm w) :=
  L.inv.preserves_adj v w hadj

end DartLeftMult

end TriangleGroupData
