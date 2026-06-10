# Group Categories — Sabidussi Coset Graphs

The graphs in this project fall into categories according to their **Sabidussi
construction** Sab(G, H): vertices are the cosets G/H of a vertex stabilizer
H ≤ G, edges given by a self-paired connecting set. The category is determined
by the stabilizer H, not by the group G — graphs with the same H form a family
sharing genus formulas, CSS code parameters, and covering relationships.

Source data: `lean/named_graphs/README.md`, `clayworth/named_graphs/README.md`,
`docs/CATALOGUE.md` §1–§6. Reference: Langer, *Symmetric Graphs and their
Quotients*, [arXiv:1306.4798](https://arxiv.org/abs/1306.4798).

---

## Category 1 — Sab(G, C₃): the genus tower (with CSS codes)

Stabilizer C₃ gives cubic graphs of size |G|/3 that carry {p,3} surface
tilings, hence a genus and a CSS quantum code (Breuckmann–Terhal,
[arXiv:1506.04029](https://arxiv.org/abs/1506.04029)).

| Graph | Foster | V | Sab(G, H) | genus | CSS |
|-------|:------:|:-:|-----------|:-----:|-----|
| Heawood | F014A | 14 | Sab(G₄₂, C₃) | 1 | [[21, 2]] |
| Möbius–Kantor | F016A | 16 | Sab(GL(2,3), C₃) | 2 | [[24, 4]] |
| Double cover of Petersen | — | 20 | Sab(PSL(2,5), C₃) | — | — |
| Klein graph | F056A | 56 | Sab(PSL(2,7), C₃) | 3 | [[84, 6]] |
| **Clayworth graph** | — | 4032 | Sab(G₂(2), C₃) | 505 | [[6048, 1010]] |

Clayworth is not drawn (too large); its quotients are (see Category 4).

## Category 2 — Sab(G, S₃): the PSL(2,q) cubic census

Stabilizer S₃ gives cubic graphs of size |G|/6.

| Graph | Foster | V | Sab(G, H) |
|-------|:------:|:-:|-----------|
| Petersen | F010A | 10 | Sab(PSL(2,5), S₃) |
| **Zhou-3 graph** | F182A | 182 | Sab(PSL(2,13), S₃) |

## Category 3 — Sab(G, S₄): the deep quotients

Stabilizer S₄ gives cubic graphs of size |G|/24 — the bottom of the
C₃ ⊂ S₃ ⊂ S₄ tower (covers q³/6 → q³/12 → q³/48, all Ramanujan; see
`docs/RESEARCH_PROGRAM.md`).

| Graph | Foster | V | Sab(G, H) |
|-------|:------:|:-:|-----------|
| Heawood | F014A | 14 | Sab(PGL(2,7), S₄) |
| Biggs–Smith | F102A | 102 | Sab(PSL(2,17), S₄) |

Note: Heawood appears in **two** categories — Sab(G₄₂, C₃) and
Sab(PGL(2,7), S₄) are different coset constructions of the same graph.

## Category 4 — G₂(2) family: hexagons and quotients

Built from the generalized hexagon GH(2,2) rather than a single Sab(G, H);
the Tutte 12-cage is **not** vertex-transitive (point/line orbits 63 + 63),
so it is not itself a Sabidussi graph.

| Graph | Foster | V | Construction |
|-------|:------:|:-:|--------------|
| Tutte 12-cage | F126A | 126 | incidence graph of GH(2,2), from G₂(2) |
| GH(2,2) collinearity graph | — | 63 | G₂(2) on 63 points |
| **Langer graph** | — | 63 (6-regular) | Sab(G₂(2), H₁₉₂) = distance-2 graph on the points of GH(2,2) |
| Dual Langer graph | — | 63 | dual of the above (`lean/named_graphs/dual-langer-63v.jpg`) |
| Meinhold-72 | — | 72 | orbital graph on G₂(2)/PSL(3,2) cosets |

The Tutte 12-cage has no {p,3} tiling (12F = 2E = 378 ⇒ F = 31.5, not an
integer) — included as parent context for the Clayworth construction.

**The Langer graph is NOT a quotient of Clayworth** (disproved in
`clayworth/PR_NOTE.md`): the Clayworth C₃ fixes zero points of GH(2,2), so
no coset map gC₃ ↦ gH₁₉₂ exists, and all small quotients of Clayworth are
degenerate. The two graphs are algebraic siblings via G₂(2), not
quotient-related. Lean proof `langer_eq_tutte12_distance2` identifies the
Langer graph with collinearity-at-distance-2 in the Tutte 12-cage.

## Category 5 — larger stabilizers: higher valency

Larger stabilizers give non-cubic Sabidussi graphs.

| Graph | Foster | V | Sab(G, H) | Valency |
|-------|:------:|:-:|-----------|:-------:|
| **Zhou-6 graph** | — | 91 | Sab(PSL(2,13), D₁₂) | 6 |
| Langer graph | — | 63 | (see Category 4) | 6 |

---

## Named graphs at a glance

| Name | Definition | Where drawn |
|------|-----------|-------------|
| **Langer graph** | Sab(G₂(2), H₁₉₂), 6-regular, 63 vertices; primitive, arc-transitive | `lean/named_graphs/langer-63v.jpg` |
| **Zhou-3 graph** | F182A = Sab(PSL(2,13), S₃), cubic, 182 vertices | `lean/named_graphs/zhou3-F182A.jpg`, `individual_graphs/auto-psl2_13.*` |
| **Zhou-6 graph** | Sab(PSL(2,13), D₁₂), 6-regular, 91 vertices | `lean/named_graphs/zhou6-91v.jpg` |
| **Clayworth graph** | Sab(G₂(2), C₃), 4032 vertices, genus 505, [[6048, 1010]] | not drawn (quotients in `clayworth/`) |
| **Meinhold-72** | orbital graph on 72 cosets G₂(2)/PSL(3,2); rank-8 Hecke algebra, B² = 13I | `individual_graphs/auto-meinhold72.*` |
| **Klein graph** | F056A = Sab(PSL(2,7), C₃), skeleton of the Klein quartic tiling | `lean/named_graphs/klein-F056A.jpg`, `individual_graphs/F056A.*` |
