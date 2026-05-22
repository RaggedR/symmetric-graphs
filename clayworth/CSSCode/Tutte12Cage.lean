/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# The Tutte 12-cage: incidence graph of GH(2,2)

The Tutte 12-cage is the incidence graph of the split Cayley hexagon GH(2,2),
the unique generalized 6-gon over 𝔽₂. It has 126 vertices (63 points + 63 lines),
189 edges, and girth 12.

It is **semisymmetric**: edge-transitive but not vertex-transitive.
Aut ≅ G₂(2) (order 12096) preserves the point/line bipartition with
two orbits of 63 vertices. No polarity exists for GH(2,2) over 𝔽₂.

This is NOT a Sabidussi graph. It is the incidence graph of a rank-2
coset geometry (Tits building of type G₂), with vertices from two coset
spaces G₂(2)/P₁ ∪ G₂(2)/P₂ joined by incidence.

## Relationship to other named graphs

  - **Heawood** = incidence graph of PG(2,2) (Fano plane, gen. 3-gon over 𝔽₂), girth 6
  - **Tutte 12-cage** = incidence graph of GH(2,2) (gen. 6-gon over 𝔽₂), girth 12
  - **Langer** = collinearity graph of GH(2,2) = distance-2 graph on the point orbit

## References

- Benson, *Minimal regular graphs of girth 8 and 12*, Canad. J. Math. 18 (1966)
- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

import CSSCode.Q62

/-! ## The Tutte 12-cage graph -/

/-- Edge list for the Tutte 12-cage (189 edges, 0-indexed).
    Vertices 0–62 are points of GH(2,2), vertices 63–125 are lines.
    Every edge connects a point to a line (bipartite). -/
private def tutte12Edges : Array (Fin 126 × Fin 126) := #[
  (0,63), (0,64), (0,65), (1,66), (1,67), (1,68),
  (2,69), (2,70), (2,71), (3,72), (3,73), (3,74),
  (4,75), (4,76), (4,77), (5,78), (5,79), (5,80),
  (6,81), (6,82), (6,83), (7,66), (7,72), (7,78),
  (8,66), (8,84), (8,85), (9,72), (9,86), (9,87),
  (10,78), (10,88), (10,89), (11,90), (11,91), (11,92),
  (12,93), (12,94), (12,95), (13,96), (13,97), (13,98),
  (14,99), (14,100), (14,101), (15,63), (15,73), (15,75),
  (16,63), (16,102), (16,103), (17,73), (17,90), (17,96),
  (18,75), (18,104), (18,105), (19,86), (19,106), (19,107),
  (20,99), (20,108), (20,109), (21,87), (21,110), (21,111),
  (22,93), (22,112), (22,113), (23,69), (23,74), (23,81),
  (24,69), (24,114), (24,115), (25,74), (25,93), (25,99),
  (26,81), (26,116), (26,117), (27,96), (27,118), (27,119),
  (28,87), (28,120), (28,121), (29,90), (29,122), (29,123),
  (30,86), (30,124), (30,125), (31,64), (31,67), (31,70),
  (32,64), (32,106), (32,108), (33,67), (33,91), (33,94),
  (34,70), (34,118), (34,120), (35,84), (35,102), (35,114),
  (36,100), (36,103), (36,124), (37,85), (37,112), (37,122),
  (38,97), (38,110), (38,115), (39,68), (39,76), (39,82),
  (40,68), (40,97), (40,100), (41,76), (41,107), (41,113),
  (42,82), (42,121), (42,123), (43,94), (43,104), (43,125),
  (44,91), (44,111), (44,116), (45,85), (45,105), (45,117),
  (46,84), (46,109), (46,119), (47,65), (47,79), (47,83),
  (48,65), (48,110), (48,112), (49,79), (49,92), (49,101),
  (50,83), (50,119), (50,125), (51,88), (51,108), (51,123),
  (52,98), (52,106), (52,117), (53,89), (53,103), (53,116),
  (54,95), (54,102), (54,121), (55,71), (55,77), (55,80),
  (56,71), (56,122), (56,124), (57,77), (57,109), (57,111),
  (58,80), (58,95), (58,98), (59,101), (59,105), (59,120),
  (60,89), (60,113), (60,118), (61,88), (61,104), (61,115),
  (62,92), (62,107), (62,114)
]

/-- Adjacency test for the Tutte 12-cage. Exposed for the dual Langer
    graph construction (distance-2 on lines). -/
def tutte12AdjBool (u v : Fin 126) : Bool :=
  tutte12Edges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

/-- The **Tutte 12-cage** as a SimpleGraph on 126 vertices. -/
def tutte12CageGraph : SimpleGraph (Fin 126) where
  Adj u v := tutte12AdjBool u v
  symm u v := by simp only [tutte12AdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [tutte12AdjBool]; revert u; native_decide⟩

instance tutte12DecAdj : DecidableRel tutte12CageGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (tutte12AdjBool u v))

/-! ## Basic properties -/

/-- The Tutte 12-cage is 3-regular. -/
theorem tutte12_regular :
    ∀ v : Fin 126, (Finset.univ.filter fun w => tutte12CageGraph.Adj v w).card = 3 := by
  native_decide

/-- The Tutte 12-cage has 189 edges. -/
theorem tutte12_edges :
    (Finset.univ.filter fun p : Fin 126 × Fin 126 =>
      p.1 < p.2 ∧ tutte12CageGraph.Adj p.1 p.2).card = 189 := by
  native_decide

/-! ## Bipartiteness

The Tutte 12-cage is bipartite: vertices 0–62 (points of GH(2,2)) form one
part, vertices 63–125 (lines) form the other. Every edge crosses the cut. -/

/-- The bipartition: point (< 63) or line (≥ 63). -/
def tutte12Side (v : Fin 126) : Fin 2 :=
  if v.val < 63 then 0 else 1

/-- The Tutte 12-cage is bipartite: every edge connects a point to a line. -/
theorem tutte12_bipartite :
    ∀ u v : Fin 126, tutte12CageGraph.Adj u v → tutte12Side u ≠ tutte12Side v := by
  native_decide

/-- Each side of the bipartition has exactly 63 vertices. -/
theorem tutte12_bipartition_sizes :
    (Finset.univ.filter fun v : Fin 126 => tutte12Side v = 0).card = 63 ∧
    (Finset.univ.filter fun v : Fin 126 => tutte12Side v = 1).card = 63 := by
  native_decide

/-! ## The Langer graph is the distance-2 graph on the point orbit

Two points p, q of GH(2,2) are collinear (adjacent in the Langer graph)
iff they are at distance 2 in the Tutte 12-cage (connected through a
shared line vertex). We prove this matches the algebraic Langer graph
definition from Q62.lean. -/

/-- Distance-2 adjacency on points: p ~ q iff they share a line.
    Equivalently, ∃ line vertex l with p ~ l and q ~ l in the Tutte 12-cage. -/
def tutte12Distance2Bool (p q : Fin 63) : Bool :=
  p != q &&
  -- Check if any line vertex (index 63..125) is adjacent to both p and q
  let pv : Fin 126 := ⟨p.val, by omega⟩
  let qv : Fin 126 := ⟨q.val, by omega⟩
  (List.finRange 63).any fun l =>
    let lv : Fin 126 := ⟨l.val + 63, by omega⟩
    tutte12AdjBool pv lv && tutte12AdjBool qv lv

/-- The Langer graph (from Q62.lean) equals the distance-2 graph
    of the Tutte 12-cage restricted to points.

    This connects the algebraic definition (Zorn product on Q(6,2))
    to the geometric definition (collinearity in GH(2,2)).

    See `G2Action.langer_eq_tutte12_distance2'` for the structural proof
    via G₂(2) transitivity: both graphs are G₂(2)-invariant and agree
    at vertex 0, so by the vertex-transitive action they agree everywhere. -/
theorem langer_eq_tutte12_distance2 :
    ∀ p q : Fin 63,
      langerSimpleGraph.Adj p q ↔ tutte12Distance2Bool p q := by
  native_decide

/-! ## No cubic double cover of the Langer graph

The Zhou graph (F182A, 182 vertices, cubic) is a Z₂ double cover of its
quotient Zhou-6 (91 vertices, 6-regular) via the block system S₃ < D₁₂.
One might ask: does the same hold for the Langer graph? That is, does there
exist a 126-vertex cubic graph covering the Langer graph via a 2:1 map?

The natural candidate is the Tutte 12-cage (126 vertices, cubic, same G₂(2)
symmetry). However, it **cannot** cover the Langer graph:

1. The Tutte 12-cage is semisymmetric (not vertex-transitive), while any
   covering of the vertex-transitive Langer graph must be vertex-transitive.

2. More fundamentally, the point stabiliser H₁₉₂ in G₂(2) (structure
   ((C₄ × C₄) ⋊ C₃) ⋊ C₂) ⋊ C₂) has three index-2 subgroups of order 96,
   but **none** yield a size-3 suborbit in the coset action on 126 points.
   The smallest suborbits are size 6 and 12 — no cubic Sabidussi graph exists.

This contrasts with PSL(2,13), where S₃ < D₁₂ gives the cubic Zhou-3 → 6-regular
Zhou-6 double cover, because S₃ has a suborbit of size 3 in the coset action.

Verification: no surjection from the Tutte 12-cage to the Langer graph
maps edges to edges (each Langer edge would need a Tutte 12-cage edge in
its preimage). We verify the weaker necessary condition: the Tutte 12-cage
point subgraph (induced on vertices 0–62) has no edges, so no point-point
edge can cover a Langer edge. -/

/-- The Tutte 12-cage has no edges between points: the induced subgraph
    on vertices 0–62 is edgeless. This obstructs any covering map to
    the Langer graph, which has 189 edges between its 63 vertices. -/
theorem tutte12_points_independent :
    ∀ p q : Fin 63,
      ¬ tutte12CageGraph.Adj ⟨p.val, by omega⟩ ⟨q.val, by omega⟩ := by
  native_decide
