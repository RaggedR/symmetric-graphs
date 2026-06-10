/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# CellularSurface from a triangle group Δ(2,3,r)

Given a finite group G with elements R (order r), S (order 3), T = RS (order 2)
generating G, we construct a `CellularSurface` with:
  - Vertices = left cosets g⟨S⟩  (|G|/3 vertices)
  - Edges = left cosets g⟨T⟩    (|G|/2 edges)
  - Faces = left cosets g⟨R⟩    (|G|/r faces, each with r edges)

The key insight: `face_closed` follows from `T = RS` and `S` preserving vertex
cosets: `vertexOf(T(d)) = vertexOf(S(R(d))) = vertexOf(R(d))`.
-/

import CSSCode.CellularSurface

/-! ## arrayToFin helper -/

/-- Convert an Array to a function on Fin, given a size proof. O(1) lookup. -/
def arrayToFin {α : Type*} {n : ℕ} (a : Array α) (h : a.size = n) (i : Fin n) : α :=
  a[i.val]'(by omega)

/-! ## TriangleGroupData -/

/-- The data of a triangle group Δ(2,3,r) acting on darts, plus coset maps.

Darts = Fin nDarts (group elements in the regular representation).
Right multiplication by R, S permutes darts. Left cosets of ⟨S⟩, ⟨T⟩, ⟨R⟩
give vertices, edges, faces. -/
structure TriangleGroupData (nDarts nV nE nF r : ℕ) where
  /-- Face rotation (order r) -/
  R : Equiv.Perm (Fin nDarts)
  /-- Vertex rotation (order 3) -/
  S : Equiv.Perm (Fin nDarts)
  /-- T = R.trans S (edge involution: T(d) = S(R(d))) -/
  T_eq_RS : ∀ d : Fin nDarts, (R.trans S) d = (R.trans S) d  -- trivially true; used as marker
  /-- Vertex coset projection: dart → vertex -/
  vertexOf : Fin nDarts → Fin nV
  /-- Edge coset projection: dart → edge -/
  edgeOf : Fin nDarts → Fin nE
  /-- Face coset projection: dart → face -/
  faceOf : Fin nDarts → Fin nF
  /-- Vertex coset representative -/
  vertexRep : Fin nV → Fin nDarts
  /-- Edge coset representative -/
  edgeRep : Fin nE → Fin nDarts
  /-- Face coset representative -/
  faceRep : Fin nF → Fin nDarts
  /-- S preserves vertex cosets -/
  vertexOf_S : ∀ d, vertexOf (S d) = vertexOf d
  /-- vertexRep round-trips -/
  vertexOf_rep : ∀ v, vertexOf (vertexRep v) = v
  /-- T preserves edge cosets -/
  edgeOf_T : ∀ d, edgeOf ((R.trans S) d) = edgeOf d
  /-- edgeRep round-trips -/
  edgeOf_rep : ∀ e, edgeOf (edgeRep e) = e
  /-- R preserves face cosets -/
  faceOf_R : ∀ d, faceOf (R d) = faceOf d
  /-- faceRep round-trips -/
  faceOf_rep : ∀ f, faceOf (faceRep f) = f
  /-- No loops -/
  edge_ne : ∀ e : Fin nE,
    vertexOf (edgeRep e) ≠ vertexOf ((R.trans S) (edgeRep e))
  /-- Trail: face boundaries visit distinct edges -/
  face_inj : ∀ (f : Fin nF) (i j : Fin r),
    edgeOf ((R ^ i.val) (faceRep f)) = edgeOf ((R ^ j.val) (faceRep f)) → i = j
  /-- T = RS is an involution: T(T(d)) = d. Holds because (RS)² = 1 in Δ(2,3,r). -/
  T_inv : ∀ d : Fin nDarts, (R.trans S) ((R.trans S) d) = d

namespace TriangleGroupData

variable {nDarts nV nE nF r : ℕ} (D : TriangleGroupData nDarts nV nE nF r)

/-- The edge involution T = R.trans S. -/
abbrev T : Equiv.Perm (Fin nDarts) := D.R.trans D.S

/-! ## Key algebraic lemma -/

/-- `vertexOf(T(d)) = vertexOf(R(d))`: the "other end" of a dart has the same
vertex as the R-successor. This is the triangle group identity:
  T(d) = S(R(d)), and S preserves vertex cosets. -/
theorem vertexOf_T_eq_R (d : Fin nDarts) :
    D.vertexOf (D.T d) = D.vertexOf (D.R d) := by
  change D.vertexOf ((D.R.trans D.S) d) = D.vertexOf (D.R d)
  simp only [Equiv.trans_apply]
  exact D.vertexOf_S (D.R d)

/-! ## CellularSurface construction -/

/-- Build a `CellularSurface` from triangle group data.

Edges: `edge_src e = vertexOf(edgeRep e)`, `edge_tgt e = vertexOf(T(edgeRep e))`.
Faces: step i of face f uses dart `R^i(faceRep f)`.
Direction: `true` iff the dart is the canonical dart of its edge.

The `face_closed` axiom is proved by `native_decide` for now.
The algebraic argument: endVertex(i) = vertexOf(T(R^i g)) = vertexOf(R^(i+1) g)
= startVertex(i+1), by `vertexOf_T_eq_R`.
-/
def toCellularSurface
    (hr : 0 < r)
    (h_face_closed : ∀ (f : Fin nF) (i : Fin r),
      let d := (D.R ^ i.val) (D.faceRep f)
      let e := D.edgeOf d
      let dir := decide (D.edgeRep e = d)
      let j : Fin r := ⟨(i.val + 1) % r, Nat.mod_lt _ hr⟩
      let d' := (D.R ^ j.val) (D.faceRep f)
      let e' := D.edgeOf d'
      let dir' := decide (D.edgeRep e' = d')
      (if dir then D.vertexOf (D.T (D.edgeRep e))
       else D.vertexOf (D.edgeRep e)) =
      (if dir' then D.vertexOf (D.edgeRep e')
       else D.vertexOf (D.T (D.edgeRep e'))))
    : CellularSurface where
  nV := nV
  nE := nE
  nF := nF
  edge_src e := D.vertexOf (D.edgeRep e)
  edge_tgt e := D.vertexOf (D.T (D.edgeRep e))
  edge_ne := D.edge_ne
  face_len _ := r
  face_len_pos _ := hr
  face_edge f i := D.edgeOf ((D.R ^ i.val) (D.faceRep f))
  face_dir f i :=
    decide (D.edgeRep (D.edgeOf ((D.R ^ i.val) (D.faceRep f))) =
            (D.R ^ i.val) (D.faceRep f))
  face_inj f := by
    intro i j h
    exact D.face_inj f i j h
  face_closed f i := by
    exact h_face_closed f i

end TriangleGroupData
