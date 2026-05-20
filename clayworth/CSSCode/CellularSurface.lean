/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Cellular Surfaces and ∂² = 0

A `CellularSurface` encodes the combinatorial data of a 2-cell embedding of a
graph on a closed surface: vertices, edges with endpoints, and faces whose
boundaries are closed directed trails.

We define the boundary operators ∂₁ (vertex-edge incidence) and ∂₂ (edge-face
boundary) as matrices over F₂ = ZMod 2 and prove the chain complex condition
∂₁ ∘ ∂₂ = 0.

## The proof

In a closed boundary walk, each vertex is visited an integer number of times.
Each visit uses exactly 2 incident edges (one entering, one leaving). Over F₂,
2k ≡ 0 for any k. The formal argument:

1. For each step i, `[src = v] + [tgt = v] = [stepStart = v] + [stepEnd = v]`
   (direction doesn't matter, by commutativity of addition).
2. The closed walk condition gives `stepEnd i = stepStart (next i)`.
3. Cyclic successor `next` is a bijection, so `∑ᵢ [stepEnd i = v] = ∑ᵢ [stepStart i = v]`.
4. Therefore `∑ᵢ d₁(v, eᵢ) = 2 · ∑ᵢ [stepStart i = v] = 0` in ZMod 2.
-/

import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

open Matrix Finset

/-! ### CellularSurface -/

/-- A `CellularSurface` is the combinatorial data of a 2-cell embedding of a
    graph on a closed surface. Faces are given as closed directed trails:
    ordered sequences of directed edges forming closed walks in the graph. -/
structure CellularSurface where
  /-- Number of vertices (0-cells) -/
  nV : ℕ
  /-- Number of edges (1-cells) -/
  nE : ℕ
  /-- Number of faces (2-cells) -/
  nF : ℕ
  /-- Source vertex of each edge -/
  edge_src : Fin nE → Fin nV
  /-- Target vertex of each edge -/
  edge_tgt : Fin nE → Fin nV
  /-- No loops: every edge has distinct endpoints -/
  edge_ne : ∀ e, edge_src e ≠ edge_tgt e
  /-- Number of edges bounding each face -/
  face_len : Fin nF → ℕ
  /-- Every face boundary is nonempty -/
  face_len_pos : ∀ f, 0 < face_len f
  /-- The sequence of edges traversed by face f's boundary -/
  face_edge : (f : Fin nF) → Fin (face_len f) → Fin nE
  /-- Traversal direction for each boundary step: `true` = src → tgt, `false` = tgt → src -/
  face_dir : (f : Fin nF) → Fin (face_len f) → Bool
  /-- Trail condition: no edge is traversed twice in a single face boundary -/
  face_inj : ∀ f, Function.Injective (face_edge f)
  /-- **Closed walk**: the end vertex of step `i` equals the start vertex of the
      cyclically next step `(i + 1) % face_len f`. This is the axiom from which
      ∂² = 0 follows. -/
  face_closed : ∀ (f : Fin nF) (i : Fin (face_len f)),
    let j : Fin (face_len f) :=
      ⟨(i.val + 1) % face_len f, Nat.mod_lt _ (face_len_pos f)⟩
    (if face_dir f i then edge_tgt (face_edge f i) else edge_src (face_edge f i)) =
    (if face_dir f j then edge_src (face_edge f j) else edge_tgt (face_edge f j))

namespace CellularSurface

variable (S : CellularSurface)

/-! ### Helper definitions -/

/-- Cyclic successor on face boundary indices. -/
def nextIdx (f : Fin S.nF) (i : Fin (S.face_len f)) : Fin (S.face_len f) :=
  ⟨(i.val + 1) % S.face_len f, Nat.mod_lt _ (S.face_len_pos f)⟩

/-- The starting vertex when traversing step `i` of face `f`'s boundary. -/
def stepStart (f : Fin S.nF) (i : Fin (S.face_len f)) : Fin S.nV :=
  if S.face_dir f i then S.edge_src (S.face_edge f i)
  else S.edge_tgt (S.face_edge f i)

/-- The ending vertex when traversing step `i` of face `f`'s boundary. -/
def stepEnd (f : Fin S.nF) (i : Fin (S.face_len f)) : Fin S.nV :=
  if S.face_dir f i then S.edge_tgt (S.face_edge f i)
  else S.edge_src (S.face_edge f i)

/-- The closed walk condition expressed via `stepEnd` and `stepStart`. -/
theorem stepEnd_eq_stepStart_next (f : Fin S.nF) (i : Fin (S.face_len f)) :
    S.stepEnd f i = S.stepStart f (S.nextIdx f i) :=
  S.face_closed f i

/-! ### Boundary operators over F₂ -/

/-- The boundary operator ∂₁: vertex–edge incidence matrix over F₂.
    Entry `(v, e)` is 1 if `v` is an endpoint of `e`, 0 otherwise.
    (Over F₂, a loop would give 1 + 1 = 0, but `edge_ne` forbids loops.) -/
def d1 : Matrix (Fin S.nV) (Fin S.nE) (ZMod 2) :=
  fun v e => (if S.edge_src e = v then 1 else 0) + (if S.edge_tgt e = v then 1 else 0)

/-- The boundary operator ∂₂: edge–face boundary matrix over F₂.
    Entry `(e, f)` counts mod 2 how many boundary steps of face `f` traverse edge `e`.
    For trails (which our faces are by `face_inj`), this is 0 or 1. -/
def d2 : Matrix (Fin S.nE) (Fin S.nF) (ZMod 2) :=
  fun e f => ∑ i : Fin (S.face_len f), if S.face_edge f i = e then 1 else 0

/-! ### The chain complex condition: ∂₁ ∘ ∂₂ = 0

The proof proceeds in five steps:

1. Expand the matrix product and swap summation order
2. Collapse the inner sum using `face_edge f i = e` as a Kronecker delta
3. Show `d₁(v, eᵢ) = [stepStart i = v] + [stepEnd i = v]` (direction-independent)
4. Reindex `∑ᵢ [stepEnd i = v]` to `∑ᵢ [stepStart i = v]` via cyclic bijection
5. Conclude `x + x = 0` in characteristic 2
-/

/-- `nextIdx` is injective: cyclic rotation by 1 preserves distinctness. -/
theorem nextIdx_injective (f : Fin S.nF) : Function.Injective (S.nextIdx f) := by
  intro a b h
  simp only [nextIdx, Fin.mk.injEq] at h
  have ha := a.isLt
  have hb := b.isLt
  ext
  -- Eliminate % by case analysis on whether a+1 and b+1 reach the boundary
  by_cases ha' : a.val + 1 < S.face_len f <;> by_cases hb' : b.val + 1 < S.face_len f
  · rw [Nat.mod_eq_of_lt ha', Nat.mod_eq_of_lt hb'] at h; omega
  · rw [Nat.mod_eq_of_lt ha', show b.val + 1 = S.face_len f from by omega, Nat.mod_self] at h
    omega
  · rw [show a.val + 1 = S.face_len f from by omega, Nat.mod_self, Nat.mod_eq_of_lt hb'] at h
    omega
  · omega

/-- `nextIdx` is bijective (injective on a finite type). -/
theorem nextIdx_bijective (f : Fin S.nF) : Function.Bijective (S.nextIdx f) :=
  (Finite.injective_iff_bijective.mp (S.nextIdx_injective f))

/-- `nextIdx` as a permutation (equivalence) on face boundary indices. -/
noncomputable def nextPerm (f : Fin S.nF) : Equiv.Perm (Fin (S.face_len f)) :=
  Equiv.ofBijective _ (S.nextIdx_bijective f)

/-- The incidence of vertex `v` with face edge `i` equals the sum of its
    start and end indicators. This holds regardless of traversal direction
    because `[src = v] + [tgt = v] = [start = v] + [end = v]` by commutativity. -/
theorem d1_eq_start_add_end (v : Fin S.nV) (f : Fin S.nF) (i : Fin (S.face_len f)) :
    S.d1 v (S.face_edge f i) =
      (if S.stepStart f i = v then 1 else 0) +
      (if S.stepEnd f i = v then 1 else 0) := by
  unfold d1 stepStart stepEnd
  rcases S.face_dir f i with _ | _
  · simp [add_comm]
  · simp

/-- **Main theorem**: ∂₁ ∘ ∂₂ = 0 over F₂.

    This is the chain complex condition for the cellular chain complex of a
    surface. It encodes the topological fact that the boundary of a boundary
    is empty: every face boundary is a cycle in the graph. -/
theorem d1_mul_d2_eq_zero : S.d1 * S.d2 = 0 := by
  ext v f
  simp only [Matrix.mul_apply, d2, Matrix.zero_apply]
  -- Step 1: Distribute the outer product into the inner sum, then swap
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  -- Step 2: For each i, collapse Σ_e d₁(v,e) · [face_edge f i = e] = d₁(v, face_edge f i)
  -- using the Kronecker delta property
  have step2 : ∀ i : Fin (S.face_len f),
      ∑ e : Fin S.nE, S.d1 v e * (if S.face_edge f i = e then 1 else 0) =
      S.d1 v (S.face_edge f i) := by
    intro i
    conv_lhs =>
      arg 2; ext e
      rw [show S.d1 v e * (if S.face_edge f i = e then 1 else 0) =
          (if e = S.face_edge f i then S.d1 v e else 0) from by
        split_ifs with h1 h2 h2 <;> simp_all [eq_comm]]
    rw [Finset.sum_ite_eq']
    simp [Finset.mem_univ]
  simp_rw [step2]
  -- Step 3: Replace d₁(v, eᵢ) with [start = v] + [end = v]
  simp_rw [S.d1_eq_start_add_end v f]
  -- Step 4: Split the sum
  rw [Finset.sum_add_distrib]
  -- Step 5: Reindex the stepEnd sum using the closed walk condition + cyclic bijection
  have step5 : ∀ i, (if S.stepEnd f i = v then (1 : ZMod 2) else 0) =
      (if S.stepStart f (S.nextIdx f i) = v then 1 else 0) := by
    intro i
    rw [← S.stepEnd_eq_stepStart_next]
  simp_rw [step5]
  -- Step 6: Reindex using cyclic bijection, then x + x = 0 in ZMod 2
  rw [show ∑ i, (if S.stepStart f (S.nextIdx f i) = v then (1 : ZMod 2) else 0) =
      ∑ i, (if S.stepStart f i = v then 1 else 0) from
    Equiv.sum_comp (S.nextPerm f) (fun j => if S.stepStart f j = v then 1 else 0)]
  have : ∀ (x : ZMod 2), x + x = 0 := by decide
  exact this _

end CellularSurface
