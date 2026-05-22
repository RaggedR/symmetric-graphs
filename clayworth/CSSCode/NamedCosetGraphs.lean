/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Named Graphs as Coset / Voltage Graphs

Algebraic definitions of the named graphs from the Foster census,
using voltage graph and coset graph constructions.

## Voltage graphs

A voltage graph on a base graph B with voltage group Γ (abelian) is a
regular Γ-covering. Vertices are V(B) × Γ, and edges lift via the voltage
assignment. Every cubic arc-transitive graph of small order arises as
a voltage graph on K₂ (two base vertices, three parallel edges with voltages).

  | Graph | Base | Γ | Voltages | V |
  |-------|------|---|----------|---|
  | Heawood (F014A) | K₂ | Z₇ | {0, 4, 6} | 14 |
  | Möbius-Kantor (F016A) | K₂ | Z₈ | {0, 1, 3} | 16 |

## References

- Gross & Tucker, *Topological Graph Theory*, 1987
- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import CSSCode.CosetGraph
import CSSCode.Heawood
import CSSCode.HeawoodSabidussi
import CSSCode.MkSabidussi
import CSSCode.QuotientGraph
import Mathlib.Data.ZMod.Basic

/-! ## Voltage graph on K₂

  A voltage graph on K₂ with cyclic group Zₘ and three voltages {v₁, v₂, v₃}
  gives a cubic graph on 2m vertices:
    - Vertices: Fin 2 × ZMod m  (point i, fibre element g)
    - Edges: (0, g) ~ (1, g + vⱼ) for j = 1, 2, 3 -/

/-- A cubic voltage graph on K₂ with voltage group Zₘ.
    Three voltages v₁, v₂, v₃ ∈ Zₘ give a cubic graph on 2m vertices.
    Edge from fibre 0 to fibre 1: (0, g) ~ (1, g + vⱼ).
    Edge from fibre 1 to fibre 0: (1, h) ~ (0, h - vⱼ). -/
def voltageGraphK2 (m : ℕ) [NeZero m]
    (v₁ v₂ v₃ : ZMod m) : SimpleGraph (Fin 2 × ZMod m) where
  Adj p q :=
    -- From fibre 0 to fibre 1: difference q.2 - p.2 is a voltage
    (p.1 = 0 ∧ q.1 = 1 ∧ q.2 - p.2 ∈ ({v₁, v₂, v₃} : Set (ZMod m))) ∨
    -- From fibre 1 to fibre 0: difference p.2 - q.2 is a voltage
    (p.1 = 1 ∧ q.1 = 0 ∧ p.2 - q.2 ∈ ({v₁, v₂, v₃} : Set (ZMod m)))
  symm := by
    intro p q hpq
    rcases hpq with ⟨hp, hq, hv⟩ | ⟨hp, hq, hv⟩
    · exact Or.inr ⟨hq, hp, hv⟩
    · exact Or.inl ⟨hq, hp, hv⟩
  loopless := ⟨fun p hp => by
    rcases hp with ⟨hp, hq, _⟩ | ⟨hp, hq, _⟩ <;> simp [hp] at hq⟩

/-! ## Heawood graph (F014A) -/

/-- The **Heawood graph** as a voltage graph on K₂ with Z₇ voltages {0, 4, 6}.

    This is the Levi graph (incidence graph) of the Fano plane PG(2,2).
    14 vertices, cubic, girth 6, Aut ≅ PGL(2,7), |Aut| = 336.

    Sabidussi: Sab(G₄₂, C₃, a) where G₄₂ = Z₇ ⋊ Z₆. -/
def heawoodVoltage : SimpleGraph (Fin 2 × ZMod 7) :=
  voltageGraphK2 7 0 4 6

instance : DecidableRel heawoodVoltage.Adj := by
  intro p q; unfold heawoodVoltage voltageGraphK2; simp only; exact instDecidableOr

/-- Heawood is 3-regular. -/
theorem heawood_voltage_regular :
    ∀ v : Fin 2 × ZMod 7,
      (Finset.univ.filter fun w => heawoodVoltage.Adj v w).card = 3 := by
  native_decide

/-- The Heawood graph has 42 directed edges (= 21 undirected edges × 2). -/
theorem heawood_voltage_directed_edges :
    (Finset.univ.filter fun p : (Fin 2 × ZMod 7) × (Fin 2 × ZMod 7) =>
      heawoodVoltage.Adj p.1 p.2).card = 42 := by
  native_decide

/-! ## Möbius-Kantor graph (F016A) -/

/-- The **Möbius-Kantor graph** as a voltage graph on K₂ with Z₈ voltages {0, 1, 3}.

    This is GP(8,3) — the generalised Petersen graph with outer 8-cycle and
    inner star with step 3.
    16 vertices, cubic, girth 6, Aut ≅ GL(2,3) : C₂, |Aut| = 96.

    Sabidussi: Sab(GL(2,3), C₃, a). -/
def mobiusKantorVoltage : SimpleGraph (Fin 2 × ZMod 8) :=
  voltageGraphK2 8 0 1 3

instance : DecidableRel mobiusKantorVoltage.Adj := by
  intro p q; unfold mobiusKantorVoltage voltageGraphK2; simp only; exact instDecidableOr

/-- Möbius-Kantor is 3-regular. -/
theorem mk_voltage_regular :
    ∀ v : Fin 2 × ZMod 8,
      (Finset.univ.filter fun w => mobiusKantorVoltage.Adj v w).card = 3 := by
  native_decide

/-- Möbius-Kantor has 48 directed edges (= 24 undirected × 2). -/
theorem mk_voltage_directed_edges :
    (Finset.univ.filter fun p : (Fin 2 × ZMod 8) × (Fin 2 × ZMod 8) =>
      mobiusKantorVoltage.Adj p.1 p.2).card = 48 := by
  native_decide

/-! ## Heawood isomorphism: voltage graph ≅ CellularSurface

The Heawood graph has two descriptions in this project:
  1. **Combinatorial** (Heawood.lean): a CellularSurface with explicit edge_src/edge_tgt
  2. **Algebraic** (above): a voltage graph on K₂ with Z₇ voltages {0, 4, 6}

We prove these define the same graph via the bijection v ↦ (v/7, v mod 7),
which maps Fano points 0–6 to fibre 0 and Fano lines 7–13 to fibre 1. -/

/-- The bijection Fin 14 → Fin 2 × ZMod 7 mapping vertex v to (v/7, v mod 7). -/
def heawoodBij (v : Fin 14) : Fin 2 × ZMod 7 :=
  (⟨v.val / 7, by omega⟩, (v.val : ZMod 7))

/-- The Heawood voltage graph is isomorphic to the CellularSurface graph.

    An edge of the CellularSurface (given by edge_src/edge_tgt) exists between
    vertices u and v iff they are adjacent in the voltage graph under the
    bijection v ↦ (v/7, v mod 7). This closes the loop between the algebraic
    and combinatorial descriptions of the Heawood graph. -/
theorem heawood_surface_eq_voltage :
    ∀ u v : Fin 14,
      (∃ e : Fin 21,
        (heawoodSurface.edge_src e = u ∧ heawoodSurface.edge_tgt e = v) ∨
        (heawoodSurface.edge_src e = v ∧ heawoodSurface.edge_tgt e = u)) ↔
      heawoodVoltage.Adj (heawoodBij u) (heawoodBij v) := by
  native_decide

/-! ## Heawood isomorphism: voltage graph ≅ Sabidussi coset graph

The Heawood graph also has a Sabidussi description:
  3. **Sabidussi** (HeawoodSabidussi.lean): `heawGraph ≃g Sab(G₄₂, H₂₄, D)`

We prove `heawGraph ≃g heawoodVoltage` via a concrete bijection,
connecting the voltage construction to the Sabidussi coset graph. -/

/-- The bijection Fin 14 → Fin 2 × ZMod 7 that makes the Sabidussi graph
    (`heawGraph`) isomorphic to the voltage graph (`heawoodVoltage`).

    heawGraph 0 → (0,0), 1 → (1,1), 2 → (0,3), 3 → (1,4), 4 → (0,1),
    5 → (1,3), 6 → (1,2), 7 → (0,5), 8 → (1,0), 9 → (0,6),
    10 → (1,6), 11 → (1,5), 12 → (0,4), 13 → (0,2). -/
private def heawSabBijData : Array (Fin 2 × ZMod 7) := #[
  (0, 0), (1, 1), (0, 3), (1, 4), (0, 1), (1, 3), (1, 2),
  (0, 5), (1, 0), (0, 6), (1, 6), (1, 5), (0, 4), (0, 2)]

private theorem heawSabBijData_size : heawSabBijData.size = 14 := by native_decide

private def heawSabBijFwd (v : Fin 14) : Fin 2 × ZMod 7 :=
  heawSabBijData[v.val]'(by have := heawSabBijData_size; omega)

/-- `heawGraph ≃g heawoodVoltage`: the Sabidussi edge-list graph is isomorphic
    to the voltage graph on K₂ with Z₇ voltages {0, 4, 6}.

    This proves adjacency is preserved in both directions. -/
theorem heawGraph_adj_iff_voltage :
    ∀ u v : Fin 14,
      heawGraph.Adj u v ↔ heawoodVoltage.Adj (heawSabBijFwd u) (heawSabBijFwd v) := by
  native_decide

/-- **The Heawood voltage graph is a Sabidussi coset graph.**

    Composing `heawGraph_adj_iff_voltage` with `heawoodSabidussiIso` from
    `HeawoodSabidussi.lean`, we obtain:
      `heawoodVoltage ≃g heawGraph ≃g Sab(G₄₂, H₂₄, D)`
    connecting all three descriptions of the Heawood graph. -/
noncomputable def heawoodVoltage_isSabidussi := heawoodSabidussiWitness.toSabidussiIso

/-! ## Möbius-Kantor isomorphism: voltage graph ≅ Sabidussi coset graph

Same treatment for the Möbius-Kantor graph: connect the voltage description
(`mobiusKantorVoltage` on `Fin 2 × ZMod 8`) to the Sabidussi description
(`mkSabGraph ≃g Sab(GL(2,3), C₃, D)` from `MkSabidussi.lean`). -/

private def mkSabBijData : Array (Fin 2 × ZMod 8) := #[
  (0, 0), (1, 2), (0, 4), (1, 6), (0, 3), (1, 1), (0, 7), (1, 5),
  (0, 2), (1, 4), (0, 6), (1, 0), (0, 5), (1, 3), (0, 1), (1, 7)]

private theorem mkSabBijData_size : mkSabBijData.size = 16 := by native_decide

private def mkSabBijFwd (v : Fin 16) : Fin 2 × ZMod 8 :=
  mkSabBijData[v.val]'(by have := mkSabBijData_size; omega)

/-- `mkSabGraph ≃g mobiusKantorVoltage`: the Sabidussi edge-list graph is
    isomorphic to the voltage graph on K₂ with Z₈ voltages {0, 1, 3}. -/
theorem mkGraph_adj_iff_voltage :
    ∀ u v : Fin 16,
      mkSabGraph.Adj u v ↔ mobiusKantorVoltage.Adj (mkSabBijFwd u) (mkSabBijFwd v) := by
  native_decide

/-- **The Möbius-Kantor voltage graph is a Sabidussi coset graph.**

    Composing `mkGraph_adj_iff_voltage` with `mkSabidussiIso`:
      `mobiusKantorVoltage ≃g mkSabGraph ≃g Sab(GL(2,3), C₃, D)` -/
noncomputable def mkVoltage_isSabidussi := mkSabidussiWitness.toSabidussiIso

/-! ## Voltage graph layout quotients

Every voltage graph on K₂ with cyclic group Zₘ has a natural quotient: the fibre
projection (fibre, g) ↦ fibre collapses each orbit to a single vertex, recovering
the base graph K₂. This is the quotient used for the concentric-ring layout. -/

instance mkQuotientDecAdj :
    DecidableRel (mobiusKantorVoltage.quotientGraph
      (Prod.fst : Fin 2 × ZMod 8 → Fin 2)).Adj := by
  intro i j; unfold SimpleGraph.quotientGraph; simp only; exact instDecidableAnd

/-- The Möbius-Kantor voltage graph quotients to K₂ under the fibre projection.
    Every pair of distinct fibres is connected (the quotient is the complete
    graph on 2 vertices). -/
theorem mk_voltage_quotient_complete :
    ∀ i j : Fin 2, i ≠ j →
    (mobiusKantorVoltage.quotientGraph
      (Prod.fst : Fin 2 × ZMod 8 → Fin 2)).Adj i j := by
  native_decide

instance heawoodQuotientDecAdj :
    DecidableRel (heawoodVoltage.quotientGraph
      (Prod.fst : Fin 2 × ZMod 7 → Fin 2)).Adj := by
  intro i j; unfold SimpleGraph.quotientGraph; simp only; exact instDecidableAnd

/-- The Heawood voltage graph quotients to K₂ under the fibre projection. -/
theorem heawood_voltage_quotient_complete :
    ∀ i j : Fin 2, i ≠ j →
    (heawoodVoltage.quotientGraph
      (Prod.fst : Fin 2 × ZMod 7 → Fin 2)).Adj i j := by
  native_decide
