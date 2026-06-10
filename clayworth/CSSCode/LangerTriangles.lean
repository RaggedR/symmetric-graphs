/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Triangles in the Langer graph and the induced G₂(2) line-action

The Langer graph (point-collinearity of GH(2,2)) has exactly 63 triangles,
one for each line of GH(2,2). Each triangle {p₁, p₂, p₃} consists of the
3 points incident to a line.

G₂(2) acts on points (G2Action.lean) preserving adjacency. Since adjacency-
preserving maps send triangles to triangles, the point-action **induces** an
action on lines. This file:

1. Defines triangles and the line↔triangle bijection
2. Defines the induced line-action generators (precomputed, verified equivariant)
3. Proves algebraically that the line-action preserves the "share a vertex"
   relation — the key step for the non-isomorphism proof in DualLanger.lean

## The algebraic argument

Two lines l₁, l₂ share a vertex iff `linePointSet l₁ ∩ linePointSet l₂` is
nonempty. If σ is a point-automorphism with induced line-action τ, then
`linePointSet (τ l) = (linePointSet l).image σ`. Since injective maps commute
with intersection (`Finset.image_inter`) and preserve non-emptiness
(`Finset.image_nonempty`), the share-vertex relation is invariant.
-/

import CSSCode.G2Action
import CSSCode.Tutte12Cage

/-! ## Section 1: Triangles in the Langer graph -/

/-- The set of all triangles (ordered triples i < j < k) in the Langer graph. -/
def langerTriangles : Finset (Fin 63 × Fin 63 × Fin 63) :=
  Finset.univ.filter fun (i, j, k) =>
    i < j ∧ j < k ∧
    langerSimpleGraph.Adj i j ∧ langerSimpleGraph.Adj j k ∧ langerSimpleGraph.Adj i k

/-- The Langer graph has exactly 63 triangles — one for each line of GH(2,2). -/
theorem langer_triangles_count : langerTriangles.card = 63 := by native_decide

/-- Each vertex of the Langer graph is in exactly 3 triangles.
    Geometrically: each point of GH(2,2) lies on exactly 3 lines. -/
theorem langer_triangles_per_vertex :
    ∀ v : Fin 63, (langerTriangles.filter fun (i, j, k) =>
      v = i ∨ v = j ∨ v = k).card = 3 := by
  native_decide

/-! ## Section 2: The line-triangle correspondence -/

/-- For a line vertex l (0-indexed as Fin 63), its 3 point neighbors in
    the Tutte 12-cage. These are the points incident to line l in GH(2,2). -/
def linePointNeighbors (l : Fin 63) : Finset (Fin 63) :=
  Finset.univ.filter fun p =>
    tutte12AdjBool ⟨p.val, by omega⟩ ⟨l.val + 63, by omega⟩

/-- Each line of GH(2,2) is incident to exactly 3 points. -/
theorem line_has_three_points :
    ∀ l : Fin 63, (linePointNeighbors l).card = 3 := by
  native_decide

/-- Precomputed line-to-triangle data: for each line index l ∈ Fin 63,
    the sorted triple (p₁, p₂, p₃) of its 3 point neighbors with p₁ < p₂ < p₃.
    Computed by GAP from the Tutte 12-cage edge list. -/
private def lineTriangleData : Array (Fin 63 × Fin 63 × Fin 63) := #[
  (⟨0, by omega⟩, ⟨15, by omega⟩, ⟨16, by omega⟩),
  (⟨0, by omega⟩, ⟨31, by omega⟩, ⟨32, by omega⟩),
  (⟨0, by omega⟩, ⟨47, by omega⟩, ⟨48, by omega⟩),
  (⟨1, by omega⟩, ⟨7, by omega⟩, ⟨8, by omega⟩),
  (⟨1, by omega⟩, ⟨31, by omega⟩, ⟨33, by omega⟩),
  (⟨1, by omega⟩, ⟨39, by omega⟩, ⟨40, by omega⟩),
  (⟨2, by omega⟩, ⟨23, by omega⟩, ⟨24, by omega⟩),
  (⟨2, by omega⟩, ⟨31, by omega⟩, ⟨34, by omega⟩),
  (⟨2, by omega⟩, ⟨55, by omega⟩, ⟨56, by omega⟩),
  (⟨3, by omega⟩, ⟨7, by omega⟩, ⟨9, by omega⟩),
  (⟨3, by omega⟩, ⟨15, by omega⟩, ⟨17, by omega⟩),
  (⟨3, by omega⟩, ⟨23, by omega⟩, ⟨25, by omega⟩),
  (⟨4, by omega⟩, ⟨15, by omega⟩, ⟨18, by omega⟩),
  (⟨4, by omega⟩, ⟨39, by omega⟩, ⟨41, by omega⟩),
  (⟨4, by omega⟩, ⟨55, by omega⟩, ⟨57, by omega⟩),
  (⟨5, by omega⟩, ⟨7, by omega⟩, ⟨10, by omega⟩),
  (⟨5, by omega⟩, ⟨47, by omega⟩, ⟨49, by omega⟩),
  (⟨5, by omega⟩, ⟨55, by omega⟩, ⟨58, by omega⟩),
  (⟨6, by omega⟩, ⟨23, by omega⟩, ⟨26, by omega⟩),
  (⟨6, by omega⟩, ⟨39, by omega⟩, ⟨42, by omega⟩),
  (⟨6, by omega⟩, ⟨47, by omega⟩, ⟨50, by omega⟩),
  (⟨8, by omega⟩, ⟨35, by omega⟩, ⟨46, by omega⟩),
  (⟨8, by omega⟩, ⟨37, by omega⟩, ⟨45, by omega⟩),
  (⟨9, by omega⟩, ⟨19, by omega⟩, ⟨30, by omega⟩),
  (⟨9, by omega⟩, ⟨21, by omega⟩, ⟨28, by omega⟩),
  (⟨10, by omega⟩, ⟨51, by omega⟩, ⟨61, by omega⟩),
  (⟨10, by omega⟩, ⟨53, by omega⟩, ⟨60, by omega⟩),
  (⟨11, by omega⟩, ⟨17, by omega⟩, ⟨29, by omega⟩),
  (⟨11, by omega⟩, ⟨33, by omega⟩, ⟨44, by omega⟩),
  (⟨11, by omega⟩, ⟨49, by omega⟩, ⟨62, by omega⟩),
  (⟨12, by omega⟩, ⟨22, by omega⟩, ⟨25, by omega⟩),
  (⟨12, by omega⟩, ⟨33, by omega⟩, ⟨43, by omega⟩),
  (⟨12, by omega⟩, ⟨54, by omega⟩, ⟨58, by omega⟩),
  (⟨13, by omega⟩, ⟨17, by omega⟩, ⟨27, by omega⟩),
  (⟨13, by omega⟩, ⟨38, by omega⟩, ⟨40, by omega⟩),
  (⟨13, by omega⟩, ⟨52, by omega⟩, ⟨58, by omega⟩),
  (⟨14, by omega⟩, ⟨20, by omega⟩, ⟨25, by omega⟩),
  (⟨14, by omega⟩, ⟨36, by omega⟩, ⟨40, by omega⟩),
  (⟨14, by omega⟩, ⟨49, by omega⟩, ⟨59, by omega⟩),
  (⟨16, by omega⟩, ⟨35, by omega⟩, ⟨54, by omega⟩),
  (⟨16, by omega⟩, ⟨36, by omega⟩, ⟨53, by omega⟩),
  (⟨18, by omega⟩, ⟨43, by omega⟩, ⟨61, by omega⟩),
  (⟨18, by omega⟩, ⟨45, by omega⟩, ⟨59, by omega⟩),
  (⟨19, by omega⟩, ⟨32, by omega⟩, ⟨52, by omega⟩),
  (⟨19, by omega⟩, ⟨41, by omega⟩, ⟨62, by omega⟩),
  (⟨20, by omega⟩, ⟨32, by omega⟩, ⟨51, by omega⟩),
  (⟨20, by omega⟩, ⟨46, by omega⟩, ⟨57, by omega⟩),
  (⟨21, by omega⟩, ⟨38, by omega⟩, ⟨48, by omega⟩),
  (⟨21, by omega⟩, ⟨44, by omega⟩, ⟨57, by omega⟩),
  (⟨22, by omega⟩, ⟨37, by omega⟩, ⟨48, by omega⟩),
  (⟨22, by omega⟩, ⟨41, by omega⟩, ⟨60, by omega⟩),
  (⟨24, by omega⟩, ⟨35, by omega⟩, ⟨62, by omega⟩),
  (⟨24, by omega⟩, ⟨38, by omega⟩, ⟨61, by omega⟩),
  (⟨26, by omega⟩, ⟨44, by omega⟩, ⟨53, by omega⟩),
  (⟨26, by omega⟩, ⟨45, by omega⟩, ⟨52, by omega⟩),
  (⟨27, by omega⟩, ⟨34, by omega⟩, ⟨60, by omega⟩),
  (⟨27, by omega⟩, ⟨46, by omega⟩, ⟨50, by omega⟩),
  (⟨28, by omega⟩, ⟨34, by omega⟩, ⟨59, by omega⟩),
  (⟨28, by omega⟩, ⟨42, by omega⟩, ⟨54, by omega⟩),
  (⟨29, by omega⟩, ⟨37, by omega⟩, ⟨56, by omega⟩),
  (⟨29, by omega⟩, ⟨42, by omega⟩, ⟨51, by omega⟩),
  (⟨30, by omega⟩, ⟨36, by omega⟩, ⟨56, by omega⟩),
  (⟨30, by omega⟩, ⟨43, by omega⟩, ⟨50, by omega⟩)
]

private theorem lineTriangleData_size : lineTriangleData.size = 63 := by native_decide

/-- Map from line index to sorted triangle: l ↦ (p₁, p₂, p₃) where
    p₁ < p₂ < p₃ are the 3 points incident to line l. -/
def lineToTriangle (l : Fin 63) : Fin 63 × Fin 63 × Fin 63 :=
  lineTriangleData[l.val]'(by have := lineTriangleData_size; omega)

/-- The lineToTriangle data matches the Tutte 12-cage adjacency. -/
theorem lineToTriangle_matches_tutte12 :
    ∀ l : Fin 63,
      let (p₁, p₂, p₃) := lineToTriangle l
      tutte12AdjBool ⟨p₁.val, by omega⟩ ⟨l.val + 63, by omega⟩ ∧
      tutte12AdjBool ⟨p₂.val, by omega⟩ ⟨l.val + 63, by omega⟩ ∧
      tutte12AdjBool ⟨p₃.val, by omega⟩ ⟨l.val + 63, by omega⟩ := by
  native_decide

/-- Every lineToTriangle output is a valid Langer triangle. -/
theorem lineToTriangle_mem :
    ∀ l : Fin 63, lineToTriangle l ∈ langerTriangles := by
  native_decide

/-- The line-to-triangle map is injective. -/
theorem lineToTriangle_injective :
    ∀ l₁ l₂ : Fin 63, lineToTriangle l₁ = lineToTriangle l₂ → l₁ = l₂ := by
  native_decide

/-- The line-to-triangle map is surjective: every Langer triangle is a line. -/
theorem lineToTriangle_surjective :
    ∀ t ∈ langerTriangles, ∃ l : Fin 63, lineToTriangle l = t := by
  native_decide

/-! ## Section 3: LinePointSet — Finset representation of line points -/

/-- The 3-element point set of line l, as a Finset. -/
def linePointSet (l : Fin 63) : Finset (Fin 63) :=
  let (p₁, p₂, p₃) := lineToTriangle l
  {p₁, p₂, p₃}

/-- linePointSet has cardinality 3 for every line. -/
theorem linePointSet_card :
    ∀ l : Fin 63, (linePointSet l).card = 3 := by native_decide

/-- Two lines share a vertex: their point sets have nonempty intersection. -/
def linesShareVertex (l₁ l₂ : Fin 63) : Prop :=
  l₁ ≠ l₂ ∧ (linePointSet l₁ ∩ linePointSet l₂).Nonempty

instance decLinesShareVertex : DecidableRel (fun l₁ l₂ => linesShareVertex l₁ l₂) :=
  fun _ _ => inferInstanceAs (Decidable (_ ∧ _))

/-! ## Section 4: G₂(2) line-action generators

The point-action of G₂(2) on the 63 points of GH(2,2) induces an action on
the 63 lines via the triangle bijection. These arrays are precomputed from
the point-action generators; Section 5 verifies they match the induced action.
-/

/-- G₂(2) line-action generator σ₁ (order 7): forward permutation. -/
def g2lineGen1Fwd : Array (Fin 63) := #[
  4, 3, 5, 11, 9, 10, 17, 15, 16, 6, 7, 8, 1, 0, 2,
  18, 19, 20, 14, 12, 13, 30, 36, 51, 52, 54, 53, 57, 24, 58, 61, 23, 62, 55,
  33, 56, 59, 27, 60, 31, 28, 43, 45, 21, 39, 22, 49, 34, 47, 37, 40, 32, 35,
  48, 46, 26, 50, 25, 41, 38, 42, 29, 44]

/-- G₂(2) line-action generator σ₁ (order 7): inverse permutation. -/
def g2lineGen1Inv : Array (Fin 63) := #[
  13, 12, 14, 1, 0, 2, 9, 10, 11, 4, 5, 3, 19, 20, 18, 7, 8, 6, 15, 16, 17,
  43, 45, 31, 28, 57, 55, 37, 40, 61, 21, 39, 51, 34, 47, 52, 22, 49, 59, 44,
  50, 58, 60, 41, 62, 42, 54, 48, 53, 46, 56, 23, 24, 26, 25, 33, 35, 27, 29,
  36, 38, 30, 32]

/-- G₂(2) line-action generator σ₂ (order 6): forward permutation. -/
def g2lineGen2Fwd : Array (Fin 63) := #[
  15, 25, 26, 36, 45, 46, 27, 60, 59, 38, 16, 29, 17, 14, 8, 37, 40, 61, 28,
  48, 53, 11, 30, 42, 57, 34, 5, 2, 1, 0, 44, 43, 23, 20, 56, 62, 51, 21, 39,
  9, 3, 35, 32, 41, 12, 52, 6, 55, 7, 50, 13, 10, 33, 4, 31, 19, 18, 58, 24,
  49, 47, 22, 54]

/-- G₂(2) line-action generator σ₂ (order 6): inverse permutation. -/
def g2lineGen2Inv : Array (Fin 63) := #[
  29, 28, 27, 40, 53, 26, 46, 48, 14, 39, 51, 21, 44, 50, 13, 0, 10, 12, 56,
  55, 33, 37, 61, 32, 58, 1, 2, 6, 18, 11, 22, 54, 42, 52, 25, 41, 3, 15, 9,
  38, 16, 43, 23, 31, 30, 4, 5, 60, 19, 59, 49, 36, 45, 20, 62, 47, 34, 24,
  57, 8, 7, 17, 35]

private theorem g2lineGen1Fwd_size : g2lineGen1Fwd.size = 63 := by decide
private theorem g2lineGen1Inv_size : g2lineGen1Inv.size = 63 := by decide
private theorem g2lineGen2Fwd_size : g2lineGen2Fwd.size = 63 := by decide
private theorem g2lineGen2Inv_size : g2lineGen2Inv.size = 63 := by decide

/-- G₂(2) line-action generator σ₁ (order 7). -/
def g2lineGen1 : Equiv.Perm (Fin 63) where
  toFun i := g2lineGen1Fwd[i.val]'(by have := g2lineGen1Fwd_size; omega)
  invFun i := g2lineGen1Inv[i.val]'(by have := g2lineGen1Inv_size; omega)
  left_inv := by intro i; revert i; decide
  right_inv := by intro i; revert i; decide

/-- G₂(2) line-action generator σ₂ (order 6). -/
def g2lineGen2 : Equiv.Perm (Fin 63) where
  toFun i := g2lineGen2Fwd[i.val]'(by have := g2lineGen2Fwd_size; omega)
  invFun i := g2lineGen2Inv[i.val]'(by have := g2lineGen2Inv_size; omega)
  left_inv := by intro i; revert i; decide
  right_inv := by intro i; revert i; decide

/-! ## Section 5: Equivariance — the line-action IS the induced action

The precomputed line-action permutations are verified to be equivariant
with the point-action: applying the line-action to a line l gives the line
whose point set is the image of l's point set under the point-action.

This is the KEY bridge: it connects the precomputed arrays to the point-action
algebraically, enabling all subsequent proofs to be structural. -/

/-- σ₁ line-action is equivariant with σ₁ point-action on point sets. -/
theorem linePointSet_gen1_equivariant :
    ∀ l : Fin 63, linePointSet (g2lineGen1 l) = (linePointSet l).image g2gen1 := by
  native_decide

/-- σ₂ line-action is equivariant with σ₂ point-action on point sets. -/
theorem linePointSet_gen2_equivariant :
    ∀ l : Fin 63, linePointSet (g2lineGen2 l) = (linePointSet l).image g2gen2 := by
  native_decide

/-- **Triangle-sharing invariance under equivariant action.**

    If σ is an injective point-map and τ is a line-map such that
    `linePointSet (τ l) = (linePointSet l).image σ`, then τ preserves
    the "share a vertex" relation.

    Proof: injective functions commute with intersection (Finset.image_inter)
    and preserve non-emptiness (Finset.image_nonempty). -/
theorem share_vertex_equivariant (σ : Equiv.Perm (Fin 63)) (τ : Equiv.Perm (Fin 63))
    (h : ∀ l : Fin 63, linePointSet (τ l) = (linePointSet l).image σ)
    (l₁ l₂ : Fin 63) :
    (linePointSet l₁ ∩ linePointSet l₂).Nonempty ↔
    (linePointSet (τ l₁) ∩ linePointSet (τ l₂)).Nonempty := by
  rw [h l₁, h l₂, ← Finset.image_inter _ _ σ.injective, Finset.image_nonempty]

/-- σ₁ line-action preserves the share-vertex relation. -/
theorem g2lineGen1_share_inv :
    ∀ l₁ l₂ : Fin 63, linesShareVertex l₁ l₂ ↔
      linesShareVertex (g2lineGen1 l₁) (g2lineGen1 l₂) := by
  intro l₁ l₂
  simp only [linesShareVertex, g2lineGen1.injective.ne_iff]
  exact and_congr_right' (share_vertex_equivariant g2gen1 g2lineGen1
    linePointSet_gen1_equivariant l₁ l₂)

/-- σ₂ line-action preserves the share-vertex relation. -/
theorem g2lineGen2_share_inv :
    ∀ l₁ l₂ : Fin 63, linesShareVertex l₁ l₂ ↔
      linesShareVertex (g2lineGen2 l₁) (g2lineGen2 l₂) := by
  intro l₁ l₂
  simp only [linesShareVertex, g2lineGen2.injective.ne_iff]
  exact and_congr_right' (share_vertex_equivariant g2gen2 g2lineGen2
    linePointSet_gen2_equivariant l₁ l₂)

/-! ## Section 6: Line-action word machinery

Mirroring G2Action.lean: a word over {σ₁, σ₁⁻¹, σ₂, σ₂⁻¹} composes to a
permutation on lines. By induction on word length, every such permutation
preserves the share-vertex relation. NO native_decide in this section. -/

/-- Line-action generator alphabet: σ₁, σ₁⁻¹, σ₂, σ₂⁻¹. -/
def applyLineGen : Fin 4 → Equiv.Perm (Fin 63)
  | 0 => g2lineGen1
  | 1 => g2lineGen1.symm
  | 2 => g2lineGen2
  | 3 => g2lineGen2.symm

/-- Apply a line-action generator word. -/
def applyLineWord : List (Fin 4) → Equiv.Perm (Fin 63)
  | [] => Equiv.refl _
  | g :: gs => (applyLineGen g).trans (applyLineWord gs)

/-- If σ preserves share-vertex, so does σ⁻¹. -/
private theorem symm_preserves_share (σ : Equiv.Perm (Fin 63))
    (h : ∀ l₁ l₂ : Fin 63, linesShareVertex l₁ l₂ ↔
      linesShareVertex (σ l₁) (σ l₂))
    (l₁ l₂ : Fin 63) :
    linesShareVertex l₁ l₂ ↔ linesShareVertex (σ.symm l₁) (σ.symm l₂) := by
  have := h (σ.symm l₁) (σ.symm l₂)
  simp only [Equiv.apply_symm_apply] at this
  exact this.symm

/-- Each individual line-action generator preserves share-vertex. -/
private theorem applyLineGen_share_inv (g : Fin 4) (l₁ l₂ : Fin 63) :
    linesShareVertex l₁ l₂ ↔ linesShareVertex (applyLineGen g l₁) (applyLineGen g l₂) := by
  fin_cases g <;> simp only [applyLineGen]
  · exact g2lineGen1_share_inv l₁ l₂
  · exact symm_preserves_share g2lineGen1 (g2lineGen1_share_inv) l₁ l₂
  · exact g2lineGen2_share_inv l₁ l₂
  · exact symm_preserves_share g2lineGen2 (g2lineGen2_share_inv) l₁ l₂

/-- Any line-action word preserves the share-vertex relation (structural induction). -/
theorem applyLineWord_share_inv (w : List (Fin 4)) :
    ∀ l₁ l₂ : Fin 63, linesShareVertex l₁ l₂ ↔
      linesShareVertex (applyLineWord w l₁) (applyLineWord w l₂) := by
  induction w with
  | nil => intro l₁ l₂; simp [applyLineWord]
  | cons g gs ih =>
    intro l₁ l₂
    simp only [applyLineWord, Equiv.trans_apply]
    rw [applyLineGen_share_inv g l₁ l₂]
    exact ih (applyLineGen g l₁) (applyLineGen g l₂)

/-! ## Section 7: Transitivity witnesses

BFS-computed words mapping line 0 → line v, for each line v.
These are the same data as in the original DualLanger.lean. -/

/-- BFS witness words for line transitivity: maps line 0 → line v. -/
private def lineWitnessWordData : Array (List (Fin 4)) := #[
  [], [2,0,2,2], [1,2,0], [3,3,1], [0], [0,3,1,2], [0,0,0], [2,1], [1,2,2],
  [0,0], [2,1,1], [3,3], [1,1,1,3], [1], [1,2], [2], [1,2,2,0], [1,1,1],
  [2,0], [0,3,0,3], [1,1], [2,2,2], [0,2,0], [0,2,1,3], [3,0,2], [3,0,3,0],
  [0,3,1], [2,2,0], [2,0,2], [3], [3,1,1], [0,0,3,0], [0,2,1,2], [1,1,3],
  [1,1,3,1], [0,2,2,0], [0,2,0,0], [2,2], [0,0,2], [0,0,3], [1,3,0], [3,0,0],
  [0,2,1], [2,2,2,1], [0,0,3,1], [0,2], [0,0,0,3], [0,3,0,0], [0,3,0],
  [1,3,3], [1,3], [2,1,1,3], [0,2,2], [0,3], [0,0,0,3,1], [0,3,1,1], [1,3,1],
  [3,0,3], [3,0], [0,0,2,1], [2,1,2], [3,1], [0,0,3,1,1]
]

private theorem lineWitnessWordData_size : lineWitnessWordData.size = 63 := by decide

/-- Line witness word: maps line 0 to line v. -/
def lineWitnessWord (v : Fin 63) : List (Fin 4) :=
  lineWitnessWordData[v.val]'(by have := lineWitnessWordData_size; omega)

/-- Each line witness word correctly maps 0 to its target line. -/
theorem lineWitnessWord_correct :
    ∀ v : Fin 63, applyLineWord (lineWitnessWord v) 0 = v := by
  decide
