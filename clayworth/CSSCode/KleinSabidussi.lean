/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Klein: algebraic Sabidussi coset graph proof

The **Klein quartic graph** (56 vertices, cubic) is Sab(PSL(2,7), C₃).
|PSL(2,7)| = 168, |C₃| = 3, 168/3 = 56 vertices.

The proof is algebraic: adjacency preservation for vertex generators
follows from the fact that the corresponding dart-level left
multiplications commute with T = RS (the edge involution). This
reduces the verification from O(56² = 3136) vertex pairs to O(168)
dart checks.
-/

import CSSCode.TriangleGroupGraph
import CSSCode.KleinTriangleData
import CSSCode.SabidussiWitness

set_option linter.style.longLine false
set_option linter.style.nativeDecide false

/-! ## Graph definition -/

/-- The Klein quartic graph, defined from the triangle group Δ(2,3,7).

Vertices are the 56 left cosets of ⟨S⟩ in PSL(2,7). Two vertices are
adjacent iff the edge involution T = RS maps a dart at one to a dart
at the other. -/
def kleinGraph : SimpleGraph (Fin 56) := kleinTriangleData.triangleGraph

/-! ## Dart-level left multiplication -/

/-- Left multiplication by R on darts commutes with T and induces the
vertex-action generator `kleinVertexGens 0`. -/
def kleinDartLR : kleinTriangleData.DartLeftMult where
  σ := kleinDartGenR
  comm_T := by native_decide
  σ_v := kleinVertexGens 0
  vertex_compat := by native_decide

/-- Left multiplication by S on darts commutes with T and induces the
vertex-action generator `kleinVertexGens 1`. -/
def kleinDartLS : kleinTriangleData.DartLeftMult where
  σ := kleinDartGenS
  comm_T := by native_decide
  σ_v := kleinVertexGens 1
  vertex_compat := by native_decide

/-! ## Sabidussi witness -/

/-- The Klein graph is a **Sabidussi coset graph**: Sab(PSL(2,7), C₃, D).

Adjacency preservation is proved algebraically via `DartLeftMult.preserves_adj`:
left multiplication commutes with T, so the witnessing dart transfers directly.
No O(nV²) enumeration needed. -/
def kleinWit : SabidussiWitness kleinGraph where
  numGens := 2
  gens := kleinVertexGens
  gens_preserve_adj := by
    intro ⟨i, hi⟩ u v hadj
    match i with
    | 0 => exact kleinDartLR.preserves_adj u v hadj
    | 1 => exact kleinDartLS.preserves_adj u v hadj
  gens_inv_preserve_adj := by
    intro ⟨i, hi⟩ u v hadj
    match i with
    | 0 => exact kleinDartLR.inv_preserves_adj u v hadj
    | 1 => exact kleinDartLS.inv_preserves_adj u v hadj
  base := 0
  witnessWords := kleinWitness
  witness_correct := by native_decide

/-- **The Klein graph is a Sabidussi coset graph.**

`kleinGraph ≃g Sab(PSL(2,7), C₃, D)` with 56 = 168/3 vertices.
The connection set D = ⟨S⟩T⟨S⟩ is a single double coset — the graph
is arc-transitive. -/
noncomputable def kleinSabidussiIso := kleinWit.toSabidussiIso
