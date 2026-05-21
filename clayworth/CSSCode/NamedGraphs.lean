/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Named Graphs from the G₂(2) Programme

Graphs from the research programme connecting exceptional Lie groups to
symmetric graph theory. These are not CSS codes but are included for context.

  - **Zhou (F182A)**: Sab(PSL(2,13), S₃), 182 vertices, cubic, girth 7.
    Vertex stabiliser S₃, so no triangle group / surface tiling.

  - **Langer**: 63 vertices, 6-regular, distance-regular.
    Collinearity graph of the generalized hexagon GH(2,2) of G₂(2).
    Quotient of the Clayworth graph.

  - **Tutte 12-cage (F126A)**: 126 vertices, cubic, girth 12.
    Incidence graph of GH(2,2). Aut = G₂(2), order 12096.
    Parent of the Clayworth construction but has no {p,3} surface tiling
    (12F = 2E = 378 gives non-integer F).
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
    Aut = G₂(2), order 12096. Vertex stabiliser S₄, edge stabiliser D₆.
    Quotient of the Clayworth graph by GH(2,2). -/
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

/-! ### Tutte 12-cage (F126A) -/

/-- The Tutte 12-cage: incidence graph of the generalized hexagon GH(2,2).
    Foster census F126A. 126 vertices (63 points + 63 lines), cubic, girth 12.
    Aut = G₂(2), order 12096.
    No {p,3} surface tiling: 12F = 2E = 378 gives F = 31.5. -/
structure Tutte12Cage where
  /-- Number of vertices = 2 × |GH(2,2) points| -/
  nV : ℕ := 126
  /-- Degree (cubic) -/
  degree : ℕ := 3
  /-- Number of edges -/
  nE : ℕ := 189
  /-- Girth -/
  girth : ℕ := 12
  /-- |Aut| = |G₂(2)| -/
  autOrder : ℕ := 12096

/-- The Tutte 12-cage parameters. -/
def tutte12Cage : Tutte12Cage := {}
