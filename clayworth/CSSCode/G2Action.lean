/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# G₂(2) action on Q(6,2) — structural proof of langer_eq_tutte12_distance2

The Langer graph is a Sabidussi graph: G₂(2) acts vertex-transitively on
its 63 vertices. We prove `langer_eq_tutte12_distance2'` structurally:

1. Two generators σ₁ (order 7), σ₂ (order 6) of G₂(2) ≤ Sym(63)
2. Each generator preserves both `langerSimpleGraph.Adj` and `tutte12Distance2Bool`
3. By induction on word length, any product of generators preserves both
4. Transitivity: for each vertex v, a BFS witness word maps 0 → v
5. Both graphs agree at vertex 0 (63 cases)
6. Transport via witness: for any pair (p,q), reduce to the vertex-0 case

The structural glue (word induction + transport) uses no `native_decide`.
Each `native_decide` encodes a named mathematical claim about the generators.
-/

import CSSCode.Q62
import CSSCode.Tutte12Cage

/-! ## Section 1: G₂(2) generators as Equiv.Perm (Fin 63)

Computed by `compute_g2_generators.py` using graph automorphism search
in the same vertex ordering as Lean's `q62Indices`.
- σ₁: order 7, from the companion matrix of x³ + x + 1 acting on Q(6,2)
- σ₂: order 6, a graph automorphism mapping vertex 0 → 10
Together they generate G₂(2) of order 12096 = |Aut(GH(2,2))|. -/

private def g2gen1Fwd : Array (Fin 63) := #[
  ⟨1, by omega⟩, ⟨3, by omega⟩, ⟨5, by omega⟩, ⟨2, by omega⟩, ⟨0, by omega⟩,
  ⟨6, by omega⟩, ⟨4, by omega⟩, ⟨23, by omega⟩, ⟨25, by omega⟩, ⟨24, by omega⟩,
  ⟨26, by omega⟩, ⟨28, by omega⟩, ⟨30, by omega⟩, ⟨27, by omega⟩, ⟨29, by omega⟩,
  ⟨31, by omega⟩, ⟨33, by omega⟩, ⟨34, by omega⟩, ⟨32, by omega⟩, ⟨35, by omega⟩,
  ⟨37, by omega⟩, ⟨38, by omega⟩, ⟨36, by omega⟩, ⟨55, by omega⟩, ⟨58, by omega⟩,
  ⟨56, by omega⟩, ⟨57, by omega⟩, ⟨60, by omega⟩, ⟨61, by omega⟩, ⟨59, by omega⟩,
  ⟨62, by omega⟩, ⟨7, by omega⟩, ⟨8, by omega⟩, ⟨9, by omega⟩, ⟨10, by omega⟩,
  ⟨12, by omega⟩, ⟨11, by omega⟩, ⟨14, by omega⟩, ⟨13, by omega⟩, ⟨15, by omega⟩,
  ⟨17, by omega⟩, ⟨16, by omega⟩, ⟨18, by omega⟩, ⟨19, by omega⟩, ⟨21, by omega⟩,
  ⟨20, by omega⟩, ⟨22, by omega⟩, ⟨39, by omega⟩, ⟨40, by omega⟩, ⟨42, by omega⟩,
  ⟨41, by omega⟩, ⟨45, by omega⟩, ⟨46, by omega⟩, ⟨44, by omega⟩, ⟨43, by omega⟩,
  ⟨47, by omega⟩, ⟨49, by omega⟩, ⟨48, by omega⟩, ⟨50, by omega⟩, ⟨51, by omega⟩,
  ⟨53, by omega⟩, ⟨52, by omega⟩, ⟨54, by omega⟩]

private def g2gen1Inv : Array (Fin 63) := #[
  ⟨4, by omega⟩, ⟨0, by omega⟩, ⟨3, by omega⟩, ⟨1, by omega⟩, ⟨6, by omega⟩,
  ⟨2, by omega⟩, ⟨5, by omega⟩, ⟨31, by omega⟩, ⟨32, by omega⟩, ⟨33, by omega⟩,
  ⟨34, by omega⟩, ⟨36, by omega⟩, ⟨35, by omega⟩, ⟨38, by omega⟩, ⟨37, by omega⟩,
  ⟨39, by omega⟩, ⟨41, by omega⟩, ⟨40, by omega⟩, ⟨42, by omega⟩, ⟨43, by omega⟩,
  ⟨45, by omega⟩, ⟨44, by omega⟩, ⟨46, by omega⟩, ⟨7, by omega⟩, ⟨9, by omega⟩,
  ⟨8, by omega⟩, ⟨10, by omega⟩, ⟨13, by omega⟩, ⟨11, by omega⟩, ⟨14, by omega⟩,
  ⟨12, by omega⟩, ⟨15, by omega⟩, ⟨18, by omega⟩, ⟨16, by omega⟩, ⟨17, by omega⟩,
  ⟨19, by omega⟩, ⟨22, by omega⟩, ⟨20, by omega⟩, ⟨21, by omega⟩, ⟨47, by omega⟩,
  ⟨48, by omega⟩, ⟨50, by omega⟩, ⟨49, by omega⟩, ⟨54, by omega⟩, ⟨53, by omega⟩,
  ⟨51, by omega⟩, ⟨52, by omega⟩, ⟨55, by omega⟩, ⟨57, by omega⟩, ⟨56, by omega⟩,
  ⟨58, by omega⟩, ⟨59, by omega⟩, ⟨61, by omega⟩, ⟨60, by omega⟩, ⟨62, by omega⟩,
  ⟨23, by omega⟩, ⟨25, by omega⟩, ⟨26, by omega⟩, ⟨24, by omega⟩, ⟨29, by omega⟩,
  ⟨27, by omega⟩, ⟨28, by omega⟩, ⟨30, by omega⟩]

private def g2gen2Fwd : Array (Fin 63) := #[
  ⟨10, by omega⟩, ⟨20, by omega⟩, ⟨29, by omega⟩, ⟨49, by omega⟩, ⟨55, by omega⟩,
  ⟨36, by omega⟩, ⟨44, by omega⟩, ⟨14, by omega⟩, ⟨25, by omega⟩, ⟨59, by omega⟩,
  ⟨40, by omega⟩, ⟨0, by omega⟩, ⟨19, by omega⟩, ⟨50, by omega⟩, ⟨35, by omega⟩,
  ⟨5, by omega⟩, ⟨7, by omega⟩, ⟨47, by omega⟩, ⟨58, by omega⟩, ⟨18, by omega⟩,
  ⟨24, by omega⟩, ⟨34, by omega⟩, ⟨41, by omega⟩, ⟨11, by omega⟩, ⟨17, by omega⟩,
  ⟨62, by omega⟩, ⟨33, by omega⟩, ⟨6, by omega⟩, ⟨28, by omega⟩, ⟨48, by omega⟩,
  ⟨45, by omega⟩, ⟨51, by omega⟩, ⟨61, by omega⟩, ⟨32, by omega⟩, ⟨42, by omega⟩,
  ⟨3, by omega⟩, ⟨8, by omega⟩, ⟨22, by omega⟩, ⟨27, by omega⟩, ⟨57, by omega⟩,
  ⟨46, by omega⟩, ⟨4, by omega⟩, ⟨21, by omega⟩, ⟨52, by omega⟩, ⟨31, by omega⟩,
  ⟨12, by omega⟩, ⟨23, by omega⟩, ⟨53, by omega⟩, ⟨60, by omega⟩, ⟨16, by omega⟩,
  ⟨26, by omega⟩, ⟨38, by omega⟩, ⟨43, by omega⟩, ⟨1, by omega⟩, ⟨9, by omega⟩,
  ⟨56, by omega⟩, ⟨37, by omega⟩, ⟨2, by omega⟩, ⟨30, by omega⟩, ⟨54, by omega⟩,
  ⟨39, by omega⟩, ⟨13, by omega⟩, ⟨15, by omega⟩]

private def g2gen2Inv : Array (Fin 63) := #[
  ⟨11, by omega⟩, ⟨53, by omega⟩, ⟨57, by omega⟩, ⟨35, by omega⟩, ⟨41, by omega⟩,
  ⟨15, by omega⟩, ⟨27, by omega⟩, ⟨16, by omega⟩, ⟨36, by omega⟩, ⟨54, by omega⟩,
  ⟨0, by omega⟩, ⟨23, by omega⟩, ⟨45, by omega⟩, ⟨61, by omega⟩, ⟨7, by omega⟩,
  ⟨62, by omega⟩, ⟨49, by omega⟩, ⟨24, by omega⟩, ⟨19, by omega⟩, ⟨12, by omega⟩,
  ⟨1, by omega⟩, ⟨42, by omega⟩, ⟨37, by omega⟩, ⟨46, by omega⟩, ⟨20, by omega⟩,
  ⟨8, by omega⟩, ⟨50, by omega⟩, ⟨38, by omega⟩, ⟨28, by omega⟩, ⟨2, by omega⟩,
  ⟨58, by omega⟩, ⟨44, by omega⟩, ⟨33, by omega⟩, ⟨26, by omega⟩, ⟨21, by omega⟩,
  ⟨14, by omega⟩, ⟨5, by omega⟩, ⟨56, by omega⟩, ⟨51, by omega⟩, ⟨60, by omega⟩,
  ⟨10, by omega⟩, ⟨22, by omega⟩, ⟨34, by omega⟩, ⟨52, by omega⟩, ⟨6, by omega⟩,
  ⟨30, by omega⟩, ⟨40, by omega⟩, ⟨17, by omega⟩, ⟨29, by omega⟩, ⟨3, by omega⟩,
  ⟨13, by omega⟩, ⟨31, by omega⟩, ⟨43, by omega⟩, ⟨47, by omega⟩, ⟨59, by omega⟩,
  ⟨4, by omega⟩, ⟨55, by omega⟩, ⟨39, by omega⟩, ⟨18, by omega⟩, ⟨9, by omega⟩,
  ⟨48, by omega⟩, ⟨32, by omega⟩, ⟨25, by omega⟩]

private theorem g2gen1Fwd_size : g2gen1Fwd.size = 63 := by native_decide
private theorem g2gen1Inv_size : g2gen1Inv.size = 63 := by native_decide
private theorem g2gen2Fwd_size : g2gen2Fwd.size = 63 := by native_decide
private theorem g2gen2Inv_size : g2gen2Inv.size = 63 := by native_decide

/-- First generator of G₂(2): order 7, from companion matrix of x³ + x + 1. -/
def g2gen1 : Equiv.Perm (Fin 63) where
  toFun i := g2gen1Fwd[i.val]'(by have := g2gen1Fwd_size; omega)
  invFun i := g2gen1Inv[i.val]'(by have := g2gen1Inv_size; omega)
  left_inv := by intro i; revert i; native_decide
  right_inv := by intro i; revert i; native_decide

/-- Second generator of G₂(2): order 6, maps vertex 0 → 10. -/
def g2gen2 : Equiv.Perm (Fin 63) where
  toFun i := g2gen2Fwd[i.val]'(by have := g2gen2Fwd_size; omega)
  invFun i := g2gen2Inv[i.val]'(by have := g2gen2Inv_size; omega)
  left_inv := by intro i; revert i; native_decide
  right_inv := by intro i; revert i; native_decide

/-! ## Section 2: Generator invariance

Each generator preserves both the Langer graph adjacency and the
Tutte 12-cage distance-2 adjacency. -/

/-- σ₁ preserves Langer adjacency. -/
theorem g2gen1_langer_inv :
    ∀ i j : Fin 63,
      langerSimpleGraph.Adj i j ↔ langerSimpleGraph.Adj (g2gen1 i) (g2gen1 j) := by
  native_decide

/-- σ₂ preserves Langer adjacency. -/
theorem g2gen2_langer_inv :
    ∀ i j : Fin 63,
      langerSimpleGraph.Adj i j ↔ langerSimpleGraph.Adj (g2gen2 i) (g2gen2 j) := by
  native_decide

/-- σ₁ preserves distance-2 adjacency. -/
theorem g2gen1_dist2_inv :
    ∀ p q : Fin 63,
      tutte12Distance2Bool p q = tutte12Distance2Bool (g2gen1 p) (g2gen1 q) := by
  native_decide

/-- σ₂ preserves distance-2 adjacency. -/
theorem g2gen2_dist2_inv :
    ∀ p q : Fin 63,
      tutte12Distance2Bool p q = tutte12Distance2Bool (g2gen2 p) (g2gen2 q) := by
  native_decide

/-! ## Section 3: Word action — structural induction

A generator word is a `List (Fin 4)` encoding a product of generators:
  0 = σ₁, 1 = σ₁⁻¹, 2 = σ₂, 3 = σ₂⁻¹.

We prove by induction on word length that any such product preserves
both graphs. This section uses NO `native_decide`. -/

/-- The four generator actions: σ₁, σ₁⁻¹, σ₂, σ₂⁻¹. -/
def applyGen : Fin 4 → Equiv.Perm (Fin 63)
  | 0 => g2gen1
  | 1 => g2gen1.symm
  | 2 => g2gen2
  | 3 => g2gen2.symm

/-- Apply a generator word (left-to-right):
    `applyWord [a, b, c] v` = `applyGen c (applyGen b (applyGen a v))`. -/
def applyWord : List (Fin 4) → Equiv.Perm (Fin 63)
  | [] => Equiv.refl _
  | g :: gs => (applyGen g).trans (applyWord gs)

/-- If σ preserves adjacency, so does σ⁻¹. -/
private theorem symm_preserves_langer (σ : Equiv.Perm (Fin 63))
    (h : ∀ i j : Fin 63,
      langerSimpleGraph.Adj i j ↔ langerSimpleGraph.Adj (σ i) (σ j))
    (i j : Fin 63) :
    langerSimpleGraph.Adj i j ↔ langerSimpleGraph.Adj (σ.symm i) (σ.symm j) := by
  have := h (σ.symm i) (σ.symm j)
  simp only [Equiv.apply_symm_apply] at this
  exact this.symm

/-- Each individual generator preserves Langer adjacency. -/
private theorem applyGen_langer_inv (g : Fin 4) (i j : Fin 63) :
    langerSimpleGraph.Adj i j ↔ langerSimpleGraph.Adj (applyGen g i) (applyGen g j) := by
  fin_cases g <;> simp only [applyGen]
  · exact g2gen1_langer_inv i j
  · exact symm_preserves_langer g2gen1 (g2gen1_langer_inv) i j
  · exact g2gen2_langer_inv i j
  · exact symm_preserves_langer g2gen2 (g2gen2_langer_inv) i j

/-- Any generator word preserves Langer adjacency (structural induction). -/
theorem applyWord_langer_inv (w : List (Fin 4)) :
    ∀ i j : Fin 63, langerSimpleGraph.Adj i j ↔
      langerSimpleGraph.Adj (applyWord w i) (applyWord w j) := by
  induction w with
  | nil => intro i j; simp [applyWord]
  | cons g gs ih =>
    intro i j
    simp only [applyWord, Equiv.trans_apply]
    rw [applyGen_langer_inv g i j]
    exact ih (applyGen g i) (applyGen g j)

/-- If σ preserves Bool equality, so does σ⁻¹. -/
private theorem symm_preserves_dist2 (σ : Equiv.Perm (Fin 63))
    (h : ∀ p q : Fin 63, tutte12Distance2Bool p q = tutte12Distance2Bool (σ p) (σ q))
    (p q : Fin 63) :
    tutte12Distance2Bool p q = tutte12Distance2Bool (σ.symm p) (σ.symm q) := by
  have := h (σ.symm p) (σ.symm q)
  simp only [Equiv.apply_symm_apply] at this
  exact this.symm

/-- Each individual generator preserves distance-2 adjacency. -/
private theorem applyGen_dist2_inv (g : Fin 4) (p q : Fin 63) :
    tutte12Distance2Bool p q = tutte12Distance2Bool (applyGen g p) (applyGen g q) := by
  fin_cases g <;> simp only [applyGen]
  · exact g2gen1_dist2_inv p q
  · exact symm_preserves_dist2 g2gen1 (g2gen1_dist2_inv) p q
  · exact g2gen2_dist2_inv p q
  · exact symm_preserves_dist2 g2gen2 (g2gen2_dist2_inv) p q

/-- Any generator word preserves distance-2 adjacency (structural induction). -/
theorem applyWord_dist2_inv (w : List (Fin 4)) :
    ∀ p q : Fin 63, tutte12Distance2Bool p q =
      tutte12Distance2Bool (applyWord w p) (applyWord w q) := by
  induction w with
  | nil => intro p q; simp [applyWord]
  | cons g gs ih =>
    intro p q
    simp only [applyWord, Equiv.trans_apply]
    rw [applyGen_dist2_inv g p q]
    exact ih (applyGen g p) (applyGen g q)

/-! ## Section 4: Transitivity witnesses

For each vertex v ∈ Fin 63, a BFS-computed word mapping 0 → v.
Computed by `compute_g2_generators.py`. Stored as an Array for
efficient lookup and decidability. -/

/-- BFS witness words stored as an Array. Each entry is the generator word
    (encoded as `List (Fin 4)`) that maps vertex 0 to vertex v. -/
private def witnessWordData : Array (List (Fin 4)) := #[
  [],           [0],          [0,0,0],      [0,0],        [1],
  [1,1,1],      [1,1],        [3,3,1],      [3,1,2],      [0,2,2,1],
  [2],          [3],          [0,0,3,0],    [1,1,3,1],    [0,0,3,3],
  [1,1,1,3],    [1,3,0],      [2,1,1],      [2,1,2,0],    [0,0,3,1],
  [0,2],        [2,1,3],      [1,3,3],      [3,3],        [0,2,2],
  [1,2,2,1],    [2,0],        [1,1,3],      [3,0],        [0,0,0,2],
  [0,2,1,3],    [0,3,0,2],    [2,0,2,2],    [2,0,2],      [2,1],
  [0,0,3],      [3,1],        [0,2,0],      [1,1,3,3],    [0,3,1,2],
  [2,2],        [1,3],        [2,1,2],      [0,0,3,1,1],  [0,3,0],
  [0,2,1],      [2,2,2],      [0,3,3],      [2,2,1],      [0,0,2],
  [1,3,1],      [0,2,1,1],    [2,2,2,1],    [0,3],        [0,2,2,1,3],
  [1,2],        [1,2,2],      [2,0,0],      [0,2,2,0],    [0,0,0,2,0],
  [0,3,1],      [3,0,0],      [0,2,1,3,0]]

private theorem witnessWordData_size : witnessWordData.size = 63 := by native_decide

/-- BFS witness words: `witnessWord v` is a generator word mapping 0 → v. -/
def witnessWord (v : Fin 63) : List (Fin 4) :=
  witnessWordData[v.val]'(by have := witnessWordData_size; omega)

/-- Each witness word correctly maps 0 to its target vertex. -/
theorem witnessWord_correct :
    ∀ v : Fin 63, applyWord (witnessWord v) 0 = v := by
  native_decide

/-! ## Section 5: Agreement at vertex 0

Both graphs agree on the neighborhood of vertex 0.
This is the base case — only 63 pairs to check. -/

/-- The Langer graph and distance-2 graph agree at vertex 0. -/
theorem graphs_agree_at_zero :
    ∀ v : Fin 63,
      langerSimpleGraph.Adj 0 v ↔ (tutte12Distance2Bool 0 v = true) := by
  native_decide

/-! ## Section 6: The main theorem — structural proof

For any pair (p, q), we transport to the vertex-0 case via the
transitivity witness for p. The proof is pure structural reasoning:
no `native_decide`, just applying the lemmas from Sections 3–5. -/

/-- Transport lemma: Langer adjacency at (σ(0), q) reduces to (0, σ⁻¹(q)). -/
private theorem transport_langer (w : List (Fin 4)) (q : Fin 63) :
    langerSimpleGraph.Adj (applyWord w 0) q ↔
      langerSimpleGraph.Adj 0 ((applyWord w).symm q) := by
  conv_lhs => rw [← Equiv.apply_symm_apply (applyWord w) q]
  exact ((applyWord_langer_inv w) 0 ((applyWord w).symm q)).symm

/-- Transport lemma: distance-2 adjacency at (σ(0), q) reduces to (0, σ⁻¹(q)). -/
private theorem transport_dist2 (w : List (Fin 4)) (q : Fin 63) :
    tutte12Distance2Bool (applyWord w 0) q =
      tutte12Distance2Bool 0 ((applyWord w).symm q) := by
  conv_lhs => rw [← Equiv.apply_symm_apply (applyWord w) q]
  exact ((applyWord_dist2_inv w) 0 ((applyWord w).symm q)).symm

/-- **The Langer graph equals the distance-2 graph of the Tutte 12-cage
    restricted to points** — proved structurally via G₂(2) transitivity.

    The argument:
    1. Let σ = applyWord (witnessWord p), so σ(0) = p
    2. Adj p q ↔ Adj 0 (σ⁻¹ q)  [by transport_langer]
    3. dist2 p q = dist2 0 (σ⁻¹ q)  [by transport_dist2]
    4. Adj 0 (σ⁻¹ q) ↔ dist2 0 (σ⁻¹ q)  [by graphs_agree_at_zero] -/
theorem langer_eq_tutte12_distance2' :
    ∀ p q : Fin 63,
      langerSimpleGraph.Adj p q ↔ (tutte12Distance2Bool p q = true) := by
  intro p q
  -- Rewrite p as σ(0) where σ is the witness permutation
  rw [← witnessWord_correct p]
  -- Transport both sides to vertex 0
  rw [transport_langer, transport_dist2]
  -- Apply the base case
  exact graphs_agree_at_zero _

/-! ## Section 7: λ=1 for the Langer graph — by G₂(2) transport

Each pair of adjacent vertices has exactly 1 common neighbor.
Proved by verifying at vertex 0 (`decide`) and transporting via the
G₂(2) action. The transport uses the fact that any graph automorphism
bijects common-neighbor sets, preserving their cardinality. -/

/-- A Langer graph automorphism bijects common-neighbor sets. -/
private theorem card_commonNeighbors_perm (σ : Equiv.Perm (Fin 63))
    (hσ : ∀ i j : Fin 63,
      langerSimpleGraph.Adj i j ↔ langerSimpleGraph.Adj (σ i) (σ j))
    (p q : Fin 63) :
    (Finset.univ.filter fun w =>
      langerSimpleGraph.Adj p w ∧ langerSimpleGraph.Adj q w).card =
    (Finset.univ.filter fun w =>
      langerSimpleGraph.Adj (σ p) w ∧ langerSimpleGraph.Adj (σ q) w).card := by
  -- The image of the p,q filter set under σ equals the σp,σq filter set
  suffices h : (Finset.univ.filter fun w =>
      langerSimpleGraph.Adj (σ p) w ∧ langerSimpleGraph.Adj (σ q) w) =
    (Finset.univ.filter fun w =>
      langerSimpleGraph.Adj p w ∧ langerSimpleGraph.Adj q w).image σ by
    rw [h, Finset.card_image_of_injective _ σ.injective]
  ext w
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_image]
  constructor
  · intro ⟨h1, h2⟩
    exact ⟨σ.symm w,
      ⟨(hσ p (σ.symm w)).mpr (by rwa [σ.apply_symm_apply]),
       (hσ q (σ.symm w)).mpr (by rwa [σ.apply_symm_apply])⟩,
      σ.apply_symm_apply w⟩
  · rintro ⟨u, ⟨h1, h2⟩, rfl⟩
    exact ⟨(hσ p u).mp h1, (hσ q u).mp h2⟩

/-- λ=1 at vertex 0: each neighbor of 0 has exactly 1 common neighbor with 0. -/
private theorem langer_lambda_at_zero :
    ∀ q : Fin 63, langerSimpleGraph.Adj 0 q →
      (Finset.univ.filter fun w =>
        langerSimpleGraph.Adj 0 w ∧ langerSimpleGraph.Adj q w).card = 1 := by
  decide

/-- **λ=1 for the Langer graph**: each pair of adjacent vertices has exactly 1
    common neighbor. This is the defining property of the collinearity graph of
    a generalized hexagon: two collinear points lie on a unique line, which
    contains exactly one other point.

    Proved by G₂(2) transport: verify at vertex 0 (6 neighbors, `decide`),
    then transport to all vertices via the group action. -/
theorem langer_lambda :
    ∀ p q : Fin 63, langerSimpleGraph.Adj p q →
      (Finset.univ.filter fun w =>
        langerSimpleGraph.Adj p w ∧ langerSimpleGraph.Adj q w).card = 1 := by
  intro p q hadj
  -- Let σ = applyWord(witnessWord p), so σ(0) = p
  -- σ⁻¹ preserves adjacency: rewrite using applyWord_langer_inv
  have hσ_adj : ∀ i j : Fin 63, langerSimpleGraph.Adj i j ↔
      langerSimpleGraph.Adj ((applyWord (witnessWord p)).symm i)
        ((applyWord (witnessWord p)).symm j) := by
    intro i j
    have := applyWord_langer_inv (witnessWord p)
      ((applyWord (witnessWord p)).symm i) ((applyWord (witnessWord p)).symm j)
    simp only [Equiv.apply_symm_apply] at this
    exact this.symm
  -- Transport: common-neighbor count is invariant under σ⁻¹
  rw [card_commonNeighbors_perm (applyWord (witnessWord p)).symm hσ_adj p q]
  -- σ⁻¹(p) = 0 since σ(0) = p
  have h0 : (applyWord (witnessWord p)).symm p = 0 := by
    rw [Equiv.symm_apply_eq]; exact (witnessWord_correct p).symm
  -- Adjacency transports to vertex 0
  have hadj0 := h0 ▸ (hσ_adj p q).mp hadj
  rw [h0]
  exact langer_lambda_at_zero _ hadj0

/-! ## Section 8: Primitivity of the G₂(2) action

The action of G₂(2) on the 63 points of GH(2,2) is **primitive**: the point
stabiliser H₁₉₂ (order 192 = 12096/63) is a maximal subgroup of G₂(2).

We verify this directly using Atkinson's algorithm (1975): for each v ≠ 0,
compute the smallest block containing {0, v} by propagating equivalences
through generators. A block B satisfies: ∀ g ∈ G, gB = B ∨ gB ∩ B = ∅.
Equivalently: if i ~ j (same block) then g(i) ~ g(j). Starting from 0 ~ v,
we propagate via a queue of newly-merged pairs until stable. If the block
is always all of Ω, no non-trivial block system exists.

The queue-based approach processes only the O(n) actual merges rather than
scanning all elements at each round, keeping the computation small. -/

/-- All four generator permutations: σ₁, σ₁⁻¹, σ₂, σ₂⁻¹. -/
private def g2AllGens : List (Equiv.Perm (Fin 63)) :=
  [g2gen1, g2gen1.symm, g2gen2, g2gen2.symm]

/-- Merge two equivalence classes: replace all occurrences of `hi` with `lo`. -/
private def mergeRep (f : Fin 63 → Fin 63) (lo hi : Fin 63) : Fin 63 → Fin 63 :=
  fun i => let fi := f i; if fi = hi then lo else fi

/-- Atkinson's algorithm: compute the smallest block containing {0, v}.

    State: `(queue, partition)` where `queue` is a list of newly-merged pairs
    and `partition : Fin 63 → Fin 63` maps each element to its class representative.

    At each step, pop a pair (x, y) from the queue. For each generator g,
    check if g(x) and g(y) have the same representative. If not, merge their
    classes and enqueue the new pair. Terminate when the queue is empty.

    `fuel` bounds the number of steps (at most 62 merges × 4 generators = 248
    enqueued pairs, so fuel = 300 is generous). -/
private def atkinson (gens : List (Equiv.Perm (Fin 63)))
    (queue : List (Fin 63 × Fin 63)) (part : Fin 63 → Fin 63) :
    Nat → Fin 63 → Fin 63
  | 0 => part
  | fuel + 1 =>
    match queue with
    | [] => part
    | (x, y) :: rest =>
      -- For each generator, check if g(x) and g(y) need merging
      let Acc := List (Fin 63 × Fin 63) × (Fin 63 → Fin 63)
      let (newPairs, part') := gens.foldl (fun (acc : Acc) g =>
        let (pairs, p) := acc
        let rgx := p (g x)
        let rgy := p (g y)
        if rgx = rgy then (pairs, p)
        else
          let lo := min rgx rgy
          let hi := max rgx rgy
          ((lo, hi) :: pairs, mergeRep p lo hi)) ([], part)
      atkinson gens (rest ++ newPairs) part' fuel

/-- Compute the block class of vertex 0 after Atkinson's algorithm starting
    from {0, v}. Returns `true` iff all vertices have representative 0
    (i.e., the block is all of Ω). -/
private def blockIsFullBool (v : Fin 63) : Bool :=
  let initPart : Fin 63 → Fin 63 := fun i => if i = v then 0 else i
  let finalPart := atkinson g2AllGens [(0, v)] initPart 300
  decide (∀ w : Fin 63, finalPart w = 0)

/-- **Primitivity of the G₂(2) action on the 63 points of GH(2,2).**

    For every vertex v ≠ 0, Atkinson's algorithm starting from 0 ~ v
    merges all 63 vertices into a single equivalence class. This means
    no non-trivial block system exists, so the point stabiliser H₁₉₂ ≤ G₂(2)
    is a **maximal subgroup** — there is no subgroup K with H₁₉₂ < K < G₂(2).

    The Langer graph is primitive and admits no cubic double cover
    (see `tutte12_points_independent`). -/
theorem langer_action_primitive :
    ∀ v : Fin 63, v ≠ 0 → blockIsFullBool v = true := by
  native_decide
