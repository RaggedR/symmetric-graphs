/-
Test: use attribute [irreducible] to seal witness data AFTER native_decide
but BEFORE structure assembly.
-/

import CSSCode.TriangleGroupGraph
import CSSCode.ClayworthTriangleData
import CSSCode.SabidussiWitness

set_option linter.style.longLine false
set_option linter.style.nativeDecide false
set_option linter.style.maxHeartbeats false

@[irreducible] def clayworthGraph' : SimpleGraph (Fin 4032) :=
  clayworthTriangleData.triangleGraph

private theorem clayworthAdj_iff' (u v : Fin 4032) :
    clayworthGraph'.Adj u v ↔
    (u ≠ v ∧ ∃ d : Fin 12096,
      clayworthTriangleData.vertexOf d = u ∧
      clayworthTriangleData.vertexOf (clayworthTriangleData.T d) = v) := by
  unfold clayworthGraph' TriangleGroupData.triangleGraph
  exact Iff.rfl

-- Dart-level theorems
private theorem test_LR_comm_T :
    ∀ d, clayworthDartGenR (clayworthTriangleData.T d) =
         clayworthTriangleData.T (clayworthDartGenR d) := by native_decide

private theorem test_LR_vertex :
    ∀ d, clayworthTriangleData.vertexOf (clayworthDartGenR d) =
         (clayworthVertexGens 0) (clayworthTriangleData.vertexOf d) := by native_decide

private theorem test_LS_comm_T :
    ∀ d, clayworthDartGenS (clayworthTriangleData.T d) =
         clayworthTriangleData.T (clayworthDartGenS d) := by native_decide

private theorem test_LS_vertex :
    ∀ d, clayworthTriangleData.vertexOf (clayworthDartGenS d) =
         (clayworthVertexGens 1) (clayworthTriangleData.vertexOf d) := by native_decide

-- Adjacency preservation (algebraic)
private theorem preserves_adj_of_dart'
    (σ_d : Equiv.Perm (Fin 12096)) (σ_v : Equiv.Perm (Fin 4032))
    (comm_T : ∀ d, σ_d (clayworthTriangleData.T d) = clayworthTriangleData.T (σ_d d))
    (compat : ∀ d, clayworthTriangleData.vertexOf (σ_d d) = σ_v (clayworthTriangleData.vertexOf d))
    (u v : Fin 4032) (hadj : clayworthGraph'.Adj u v) :
    clayworthGraph'.Adj (σ_v u) (σ_v v) := by
  rw [clayworthAdj_iff'] at hadj ⊢
  obtain ⟨hne, d, hu, hv⟩ := hadj
  exact ⟨σ_v.injective.ne hne, σ_d d, by rw [compat, hu],
    by change clayworthTriangleData.vertexOf (clayworthTriangleData.T (σ_d d)) = _
       rw [← comm_T, compat, hv]⟩

private theorem R_preserves_adj' (u v : Fin 4032) (hadj : clayworthGraph'.Adj u v) :
    clayworthGraph'.Adj ((clayworthVertexGens 0) u) ((clayworthVertexGens 0) v) :=
  preserves_adj_of_dart' clayworthDartGenR (clayworthVertexGens 0)
    test_LR_comm_T test_LR_vertex u v hadj

private theorem S_preserves_adj' (u v : Fin 4032) (hadj : clayworthGraph'.Adj u v) :
    clayworthGraph'.Adj ((clayworthVertexGens 1) u) ((clayworthVertexGens 1) v) :=
  preserves_adj_of_dart' clayworthDartGenS (clayworthVertexGens 1)
    test_LS_comm_T test_LS_vertex u v hadj

private theorem inv_comm_T' (σ : Equiv.Perm (Fin 12096))
    (h : ∀ d, σ (clayworthTriangleData.T d) = clayworthTriangleData.T (σ d)) :
    ∀ d, σ.symm (clayworthTriangleData.T d) = clayworthTriangleData.T (σ.symm d) := by
  intro d; apply σ.injective; rw [Equiv.apply_symm_apply, h, Equiv.apply_symm_apply]

private theorem inv_compat' (σ_d : Equiv.Perm (Fin 12096)) (σ_v : Equiv.Perm (Fin 4032))
    (h : ∀ d, clayworthTriangleData.vertexOf (σ_d d) = σ_v (clayworthTriangleData.vertexOf d)) :
    ∀ d, clayworthTriangleData.vertexOf (σ_d.symm d) = σ_v.symm (clayworthTriangleData.vertexOf d) := by
  intro d; apply σ_v.injective
  rw [Equiv.apply_symm_apply, ← h, Equiv.apply_symm_apply]

private theorem Ri_preserves_adj' (u v : Fin 4032) (hadj : clayworthGraph'.Adj u v) :
    clayworthGraph'.Adj ((clayworthVertexGens 0).symm u) ((clayworthVertexGens 0).symm v) :=
  preserves_adj_of_dart' clayworthDartGenR.symm (clayworthVertexGens 0).symm
    (inv_comm_T' _ test_LR_comm_T) (inv_compat' _ _ test_LR_vertex) u v hadj

private theorem Si_preserves_adj' (u v : Fin 4032) (hadj : clayworthGraph'.Adj u v) :
    clayworthGraph'.Adj ((clayworthVertexGens 1).symm u) ((clayworthVertexGens 1).symm v) :=
  preserves_adj_of_dart' clayworthDartGenS.symm (clayworthVertexGens 1).symm
    (inv_comm_T' _ test_LS_comm_T) (inv_compat' _ _ test_LS_vertex) u v hadj

-- Step 1: Prove witness_correct with native_decide (data still reducible)
set_option maxHeartbeats 800000 in
private theorem witness_correct' :
    ∀ v : Fin 4032, applyWordGen clayworthVertexGens (clayworthWitness v) 0 = v := by native_decide

-- Step 2: SEAL the witness data so the kernel doesn't normalize it
-- during structure assembly
attribute [irreducible] clayworthWitness

-- Step 3: Build the structure (kernel won't unfold clayworthWitness)
def clayworthWit' : SabidussiWitness clayworthGraph' where
  numGens := 2
  gens := clayworthVertexGens
  gens_preserve_adj i := by
    match i with
    | ⟨0, _⟩ => exact R_preserves_adj'
    | ⟨1, _⟩ => exact S_preserves_adj'
  gens_inv_preserve_adj i := by
    match i with
    | ⟨0, _⟩ => exact Ri_preserves_adj'
    | ⟨1, _⟩ => exact Si_preserves_adj'
  base := 0
  witnessWords := clayworthWitness
  witness_correct := witness_correct'
