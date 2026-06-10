/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Clayworth: algebraic Sabidussi coset graph proof

The **Clayworth graph** (4032 vertices, cubic) is Sab(G₂(2), C₃).
|G₂(2)| = 12096, |C₃| = 3, 12096/3 = 4032 vertices.

The proof is algebraic: adjacency preservation for vertex generators
follows from the fact that the corresponding dart-level left
multiplications commute with T = RS (the edge involution). This
reduces the verification from O(4032² = 16M) vertex pairs to O(12096)
dart checks — the difference between CI timeout and a fast build.

## Context

This is one of four CSS surface codes:
  - Heawood (F014A), Möbius-Kantor (F016A): Sabidussi via voltage
    graphs on K₂ (mathlib4 PR #39651)
  - Klein (F056A): Sabidussi via Δ(2,3,7) triangle group
  - Clayworth (4032v): Sabidussi via Δ(2,3,12) triangle group (this file)
-/

import CSSCode.TriangleGroupGraph
import CSSCode.ClayworthTriangleData
import CSSCode.SabidussiWitness

set_option linter.style.longLine false
set_option linter.style.nativeDecide false

/-! ## Graph definition -/

/-- The **Clayworth graph**, defined from the triangle group Δ(2,3,12).
Marked `@[irreducible]` so the kernel does not unfold the existential
adjacency (over Fin 12096) when type-checking downstream structures. -/
@[irreducible] def clayworthGraph : SimpleGraph (Fin 4032) :=
  clayworthTriangleData.triangleGraph

/-- Opaque bridge: unfold `clayworthGraph.Adj` to the existential form.
This is an opaque `theorem`, so the kernel sees the TYPE but not the
proof term (which contains the irreducible unfolding). -/
private theorem clayworthAdj_iff (u v : Fin 4032) :
    clayworthGraph.Adj u v ↔
    (u ≠ v ∧ ∃ d : Fin 12096,
      clayworthTriangleData.vertexOf d = u ∧
      clayworthTriangleData.vertexOf (clayworthTriangleData.T d) = v) := by
  unfold clayworthGraph TriangleGroupData.triangleGraph
  exact Iff.rfl

/-! ## Dart-level algebraic proofs (all O(12096) native_decide) -/

private theorem clayworthLR_comm_T :
    ∀ d, clayworthDartGenR (clayworthTriangleData.T d) =
         clayworthTriangleData.T (clayworthDartGenR d) := by native_decide

private theorem clayworthLR_vertex_compat :
    ∀ d, clayworthTriangleData.vertexOf (clayworthDartGenR d) =
         (clayworthVertexGens 0) (clayworthTriangleData.vertexOf d) := by native_decide

private theorem clayworthLS_comm_T :
    ∀ d, clayworthDartGenS (clayworthTriangleData.T d) =
         clayworthTriangleData.T (clayworthDartGenS d) := by native_decide

private theorem clayworthLS_vertex_compat :
    ∀ d, clayworthTriangleData.vertexOf (clayworthDartGenS d) =
         (clayworthVertexGens 1) (clayworthTriangleData.vertexOf d) := by native_decide

/-! ## Adjacency preservation (algebraic, via opaque bridge) -/

/-- Helper: prove adjacency preservation from dart-level commutation.
All access to the graph internals goes through `clayworthAdj_iff` (opaque). -/
private theorem preserves_adj_of_dart
    (σ_d : Equiv.Perm (Fin 12096)) (σ_v : Equiv.Perm (Fin 4032))
    (comm_T : ∀ d, σ_d (clayworthTriangleData.T d) = clayworthTriangleData.T (σ_d d))
    (compat : ∀ d, clayworthTriangleData.vertexOf (σ_d d) = σ_v (clayworthTriangleData.vertexOf d))
    (u v : Fin 4032) (hadj : clayworthGraph.Adj u v) :
    clayworthGraph.Adj (σ_v u) (σ_v v) := by
  rw [clayworthAdj_iff] at hadj ⊢
  obtain ⟨hne, d, hu, hv⟩ := hadj
  exact ⟨σ_v.injective.ne hne, σ_d d, by rw [compat, hu],
    by change clayworthTriangleData.vertexOf (clayworthTriangleData.T (σ_d d)) = _
       rw [← comm_T, compat, hv]⟩

private theorem clayworthR_preserves_adj (u v : Fin 4032) (hadj : clayworthGraph.Adj u v) :
    clayworthGraph.Adj ((clayworthVertexGens 0) u) ((clayworthVertexGens 0) v) :=
  preserves_adj_of_dart clayworthDartGenR (clayworthVertexGens 0)
    clayworthLR_comm_T clayworthLR_vertex_compat u v hadj

private theorem clayworthS_preserves_adj (u v : Fin 4032) (hadj : clayworthGraph.Adj u v) :
    clayworthGraph.Adj ((clayworthVertexGens 1) u) ((clayworthVertexGens 1) v) :=
  preserves_adj_of_dart clayworthDartGenS (clayworthVertexGens 1)
    clayworthLS_comm_T clayworthLS_vertex_compat u v hadj

/-- Inverse commutation: algebraic, no native_decide. -/
private theorem inv_comm_T (σ : Equiv.Perm (Fin 12096))
    (h : ∀ d, σ (clayworthTriangleData.T d) = clayworthTriangleData.T (σ d)) :
    ∀ d, σ.symm (clayworthTriangleData.T d) = clayworthTriangleData.T (σ.symm d) := by
  intro d; apply σ.injective; rw [Equiv.apply_symm_apply, h, Equiv.apply_symm_apply]

private theorem inv_compat (σ_d : Equiv.Perm (Fin 12096)) (σ_v : Equiv.Perm (Fin 4032))
    (h : ∀ d, clayworthTriangleData.vertexOf (σ_d d) = σ_v (clayworthTriangleData.vertexOf d)) :
    ∀ d, clayworthTriangleData.vertexOf (σ_d.symm d) = σ_v.symm (clayworthTriangleData.vertexOf d) := by
  intro d; apply σ_v.injective
  rw [Equiv.apply_symm_apply, ← h, Equiv.apply_symm_apply]

private theorem clayworthRi_preserves_adj (u v : Fin 4032) (hadj : clayworthGraph.Adj u v) :
    clayworthGraph.Adj ((clayworthVertexGens 0).symm u) ((clayworthVertexGens 0).symm v) :=
  preserves_adj_of_dart clayworthDartGenR.symm (clayworthVertexGens 0).symm
    (inv_comm_T _ clayworthLR_comm_T) (inv_compat _ _ clayworthLR_vertex_compat) u v hadj

private theorem clayworthSi_preserves_adj (u v : Fin 4032) (hadj : clayworthGraph.Adj u v) :
    clayworthGraph.Adj ((clayworthVertexGens 1).symm u) ((clayworthVertexGens 1).symm v) :=
  preserves_adj_of_dart clayworthDartGenS.symm (clayworthVertexGens 1).symm
    (inv_comm_T _ clayworthLS_comm_T) (inv_compat _ _ clayworthLS_vertex_compat) u v hadj

private theorem clayworthWitness_correct :
    ∀ v : Fin 4032, applyWordGen clayworthVertexGens (clayworthWitness v) 0 = v := by
  native_decide

/-! ## Sabidussi witness -/

/-- The Clayworth graph is a **Sabidussi coset graph**: Sab(G₂(2), C₃, D).

Adjacency preservation is proved algebraically: left multiplication
commutes with T, so the witnessing dart transfers directly. No O(nV²)
enumeration — all native_decide checks are O(nDarts) = O(12096). -/
def clayworthWit : SabidussiWitness clayworthGraph where
  numGens := 2
  gens := clayworthVertexGens
  gens_preserve_adj i := by
    match i with
    | ⟨0, _⟩ => exact clayworthR_preserves_adj
    | ⟨1, _⟩ => exact clayworthS_preserves_adj
  gens_inv_preserve_adj i := by
    match i with
    | ⟨0, _⟩ => exact clayworthRi_preserves_adj
    | ⟨1, _⟩ => exact clayworthSi_preserves_adj
  base := 0
  witnessWords := clayworthWitness
  witness_correct := clayworthWitness_correct

/-- **The Clayworth graph is a Sabidussi coset graph.** -/
noncomputable def clayworthSabidussiIso := clayworthWit.toSabidussiIso
