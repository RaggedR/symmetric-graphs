/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Named Graphs from the G₂(2) Programme

Graphs from the research programme connecting exceptional Lie groups to
symmetric graph theory. These are not CSS codes but are included for context.

  - **Zhou (F182A)**: Sab(PSL(2,13), S₃), 182 vertices, cubic, girth 7.
    Vertex stabiliser S₃, so no triangle group / surface tiling.
    **Primitive**: Aut ≅ PSL(2,13) acts primitively on vertices (no
    non-trivial block system giving a regular quotient graph).

  - **Langer**: 63 vertices, 6-regular, distance-regular.
    Collinearity graph of the generalized hexagon GH(2,2) of G₂(2).
    **Primitive**: Aut ≅ G₂(2) acts primitively on the 63 points (no
    non-trivial block system). See `Q62.lean` for the algebraic definition.

  - **Klein (F056A)**: Sab(PSL(2,7), C₃), 56 vertices, cubic, girth 7.
    **Primitive**: Aut ≅ PGL(2,7) acts primitively on vertices (no
    non-trivial block system giving a regular quotient graph).
    Surface tiling as `KleinSurface.lean` (genus 3).

  - **Tutte 12-cage**: 126 vertices (63 points + 63 lines), cubic, girth 12.
    Incidence graph of GH(2,2), the split Cayley hexagon over 𝔽₂.
    **Semisymmetric**: edge-transitive but NOT vertex-transitive (not Sabidussi).
    Aut ≅ G₂(2) (order 12096), two orbits of 63 (points and lines).
    No polarity exists for GH(2,2) over 𝔽₂.
    Distance-2 graph on points = Langer graph (proved in `Tutte12Cage.lean`).

  ## Generalized polygon hierarchy over 𝔽₂

  | Geometry | Type | Incidence graph | Girth | Aut | Polarity? |
  |----------|------|-----------------|-------|-----|-----------|
  | PG(2,2) (Fano) | gen. 3-gon | Heawood (14v) | 6 | PGL(2,7) | Yes |
  | GH(2,2) | gen. 6-gon | Tutte 12-cage (126v) | 12 | G₂(2) | No |

  PSL(2,7) ≅ GL(3,2) (Aut of Fano plane) is a subgroup of G₂(2).
  See `Tutte12Cage.lean` for the graph definition and proofs.
-/

import CSSCode.CSSCode

/-! ### Zhou graph (F182A) -/

/-- The Zhou graph: Foster census F182A, Sab(PSL(2,13), S₃).
    182 vertices, cubic, girth 7, |Aut| = PSL(2,13), order 1092.
    No surface tiling (vertex stabiliser S₃ is not cyclic). -/
structure ZhouGraph where
  /-- Number of vertices -/
  nV : ℕ := 182
  /-- Degree (3-regular) -/
  degree : ℕ := 3
  /-- Number of edges = nV × degree / 2 -/
  nE : ℕ := 273
  /-- Girth -/
  girth : ℕ := 7
  /-- |Aut| = |PSL(2,13)| -/
  autOrder : ℕ := 1092

/-- The Zhou graph parameters. -/
def zhouGraph : ZhouGraph := {}

/-! ### Langer graph -/

/-- The Langer graph: collinearity graph of the generalized hexagon GH(2,2).
    63 vertices, 6-regular, distance-regular.
    Aut = G₂(2), order 12096. Vertex stabiliser of order 192,
    arc stabiliser of order 32. Primitive (no non-trivial block system). -/
structure LangerGraph where
  /-- Number of vertices = |GH(2,2) points| -/
  nV : ℕ := 63
  /-- Degree (6-regular) -/
  degree : ℕ := 6
  /-- Number of edges = nV × degree / 2 -/
  nE : ℕ := 189
  /-- |Aut| = |G₂(2)| -/
  autOrder : ℕ := 12096

/-- The Langer graph parameters. -/
def langerGraph : LangerGraph := {}

