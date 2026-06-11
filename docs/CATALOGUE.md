# Catalogue: Groups Studied, Algorithms, Artifacts

Companion to `RESEARCH_PROGRAM.md` ("Drawable Shadows of Lie Groups over Finite
Fields") and `HECKE_RESEARCH.md` (Hecke algebra H(G,H) + irreps for every entry
below). For each group we record: the graph(s) drawn, the script that
builds/draws them, and where the rendered output lives. Last updated 2026-06-10.

All paths relative to `~/git/symmetric-graphs/` unless noted. The standard pipeline is
`src/draw_graph.py` (layout engine, GAP interface, TikZ emission) +
`src/symmetric-layout.g` (Aut, orbits, quotients); anything named `auto-*` in
`individual_graphs/` was produced by it. Gallery: `gallery.pdf`.

**Legend:** entries marked **⚑ FUTURE** are open work — either computed in GAP
but never rendered, or identified gaps awaiting any computation. They are
collected in §9 (Future Research Directions).

---

## 1. A₁ family — PSL(2,q) / PGL(2,q)  (the complete programme)

Census scripts: `src/census/psl2_cubic_census.g`, `psl2_cubic_search.g`,
`psl2_dickson_census.g`, `psl2_fast_census.g`, `psl2_census_v3.g`,
`psl2_graphs.g`. Spectral checks: `src/spectral_analysis.py`.

| Group | Stabilizer H | Graph (|V|) | Build/draw | Output |
|---|---|---|---|---|
| PSL(2,5) ≅ A₅ | S₃ | Petersen (10) | `src/draw_graph.py petersen` | `individual_graphs/auto-petersen.{tex,pdf}` (+ `auto-petersen_voltage`) |
| PSL(2,5) | C₃ | double cover of Petersen (20) | census scripts | `individual_graphs/auto-psl2_5_c3.{tex,pdf}` |
| PSL(2,7) | C₃ | 56-vertex coset graph | census scripts | `individual_graphs/auto-psl2_7_c3.{tex,pdf,jpg}` |
| PGL(2,7) = SO(3,F₇) | S₄ | Heawood (14) = Levi of PG(2,2) | `draw_graph.py heawood` | `individual_graphs/auto-heawood.{tex,pdf}`, voltage version `auto-heawood_voltage` |
| PGL(2,7) (Aut) | — | Coxeter graph (28) | `draw_graph.py coxeter` | `individual_graphs/auto-coxeter.{tex,pdf}`, `quotients/coxeter.*` |
| PSL(2,13) | S₃ | Foster 182 / coset graph | census scripts | `individual_graphs/auto-psl2_13.{tex,pdf,jpg}` |
| PSL(2,17) | S₄ | **Biggs–Smith** (102) | census scripts + pipeline | `individual_graphs/auto-biggssmith.{tex,pdf}`, `auto-biggssmith2.pdf` |
| PSL(2,19) | (cubic Ramanujan) | cr19 | census scripts | `individual_graphs/auto-cr19.{tex,pdf}` |
| PSL(2,31) | (cubic Ramanujan) | cr31 | census scripts | `individual_graphs/auto-cr31.{tex,pdf}` |
| PSL(2,7) (surface form) | — | Klein quartic tilings/quotients | `src/compute_klein_quotients.g` (in `clayworth/`), `klein-quartic/`, `src/census/klein_quartic_codes.g` | `individual_graphs/auto-klein-c7.*`, `auto-klein-minimal.*`, `klein-z2-quotient.{tex,pdf,jpg}`, `klein-quartic/` |
| GL(2,3) (Aut of Bolza) | — | Bolza surface graphs | `src/census/bolza_graphs.g`, `bolza-surface/*.py` | `bolza-surface/bolza_disk.{tex,pdf,jpg}`, `octagon_cut.*`, HTML visualizations |
| PSL(2,8) (= ²G₂(3)′!) | — | Macbeath surface graph | `src/census/macbeath_graph.g` | **⚑ FUTURE** — GAP output exists, drawing never rendered |

Tower structure: C₃ ⊂ S₃ ⊂ S₄ gives covers q³/6 → q³/12 → q³/48 (see
RESEARCH_PROGRAM.md). All confirmed Ramanujan.

## 2. A₂ family — PSL(3,q) (projective planes, generalized 3-gons)

| Group | Graph | Build/draw | Output |
|---|---|---|---|
| PSL(3,2) (⋊ duality: PGL(2,7)) | Heawood = Levi of PG(2,2) | pipeline | `individual_graphs/auto-heawood.*` |
| PSL(3,3) | Levi of PG(2,3) (26) | pipeline | `individual_graphs/auto-pg23.{tex,pdf}` |

Valency is q+1, so only q=2 is cubic; higher q drawn as incidence geometry.

## 3. B₂ = C₂ family — Sp(4,2) (generalized quadrangles)

| Group | Graph | Build/draw | Output |
|---|---|---|---|
| Sp(4,2) ≅ S₆ (Aut order 1440) | Tutte–Coxeter (30) = Levi of GQ(2,2), the (3,8)-cage | `draw_graph.py tuttecoxeter` | `individual_graphs/auto-tuttecoxeter.*`, `tuttecoxeter-clean.{tex,pdf,jpg}`, `cage-pair.*` |
| (related) | Clebsch graph computations | `hecke/hecke_clebsch.g`, notes `hecke/Clebsch.md` | **⚑ FUTURE** — GAP/notes exist, drawing never rendered |

## 4. G₂ family — G₂(2) (generalized hexagons) — ACTIVE

| Group | Graph | Build/draw | Output |
|---|---|---|---|
| G₂(2) | Tutte 12-cage (126) = incidence graph of GH(2,2) | `src/tutte12cage.json`, `src/census/tutte12cage.g`, check `clayworth/check_126_sabidussi.g` | `individual_graphs/auto-tuttecage12.{tex,pdf,jpg}`, `tutte12cage-clean.*` |
| G₂(2) on 63 points | GH(2,2) collinearity/incidence drawing | `src/build_hexagon.py`, `src/g2_cubic_search.g`, `clayworth/compute_g2_generators.py` | `individual_graphs/auto-gh22.{tex,pdf,jpg}` |
| G₂(2) / PSL(3,2), 72 cosets | Meinhold-72 orbital graph | GAP scripts in `~/git/hecke/` (`meinhold_hecke.g`, `why13_*.g`) | `individual_graphs/auto-meinhold72.{tex,pdf}` |
| — theory — | rank-8 Hecke algebra, B² = 13I, sub-hexagons of order (2,1) | `~/git/hecke/why13_fano.g`, `why13_orbitals.g`, `why13_subhexagon.g`, `why13_classes.g` | **`~/git/hecke/why13.pdf`** (expository paper), `~/git/hecke/MEINHOLD-HECKE.md` |

This answers Key Question 5 of RESEARCH_PROGRAM.md: 13 = 7 + 2·3 via
NNᵀ = 7I + 2·SRG(36,21,12,12) on the 36+36 chiral Fano structures
(= sub-hexagons of order (2,1) = distance-2 ovoids).

## 5. Sporadic groups

| Group | Status | Script | Output |
|---|---|---|---|
| M₁₁ | cubic Sabidussi search attempted | `src/census/m11_cubic_search.g`, `m11_component.g` | **⚑ FUTURE** — GAP search only; no drawing, search outcome undocumented |

## 6. Small / classical named graphs (Foster census baseline, ≤ 30 vertices)

All drawn by `src/draw_graph.py <name>` / `src/draw_all.py`; outputs in
`individual_graphs/auto-<name>.{tex,pdf}`, quotient versions in `quotients/`,
TikZ sources mirrored in `named_graphs_tikz/`.

| Graph | Aut group | Notes |
|---|---|---|
| K₄ (4) | S₄ | F004A |
| K₃,₃ (6) | S₃≀C₂ (72) | F006A; Levi of the "digon" geometry |
| Cube (8) | S₄×C₂ | F008A |
| Petersen (10) | S₅ ≅ PGL(2,5) | F010A; Kneser K(5,2) |
| Heawood (14) | PGL(2,7) | F014A; PG(2,2) |
| Möbius–Kantor (16) | order 96 | F016A; voltage version `auto-mk_voltage` |
| Pappus (18) | order 216 | F018A; AG(2,3)-configuration |
| Desargues (20) | S₅×C₂ | F020A; Levi of Desargues configuration |
| Dodecahedron (20) | A₅×C₂ | F020B |
| Nauru (24) | order 144 | F024A; generalized Petersen GP(12,5) |
| F26A (26) | order 78 | |
| Coxeter (28) | PGL(2,7) | F028A |
| Tutte–Coxeter (30) | order 1440 (Aut S₆) | F030A; GQ(2,2) |
| F056A (56) | — | `F056A.{pdf,jpg}`, `F056A-minimal.*` |
| Biggs–Smith (102) | PSL(2,17) | F102A; distance-transitive |
| Tutte 12-cage (126) | G₂(2) (order 12096), **not vertex-transitive**: point/line orbits 63+63 (no duality at q=2; H(q) self-dual only for q=3^h) | F126A; GH(2,2); corrected 2026-06-10 via `~/git/hecke/hecke_census.g` |

## 7. Quotients & surfaces side-programme

| Object | Script | Output |
|---|---|---|
| Orbit quotients of all gallery graphs | `src/draw_quotients.py`, `src/enumerate-quotients.g` | `quotients/*.pdf`, `quotients/all_quotients.pdf` |
| Zhou ℤ₂-quotients | (clayworth) | `individual_graphs/zhou-z2-quotient.*`, `zhou-quotient-spectral.*` |
| Clayworth cellular surfaces / codes | `src/clayworth_build.g`, `src/clayworth_codes.py`, `clayworth/` | `graphviz-clayworth-quotients/individual_graphs/` (worktree) |
| Error-correcting codes from graphs | `src/graph_codes.py` | console tables |
| Bolza / Klein / Macbeath surfaces | `bolza-surface/`, `klein-quartic/`, `src/census/macbeath_graph.g` | `bolza-surface/*.{pdf,jpg,html}` |

---

## 8. Gap analysis: what the classification of finite simple groups says we've missed

CFSG: cyclic Cₚ, alternating Aₙ (n ≥ 5), 16 families of Lie type, 26 sporadics.
Our coverage and the natural missing targets:

### Covered (at least one graph)
- **A₁**: PSL(2,q) for q = 5, 7, 13, 17, 19, 31 (+ PGL forms; 8 via Macbeath); q = 9, 11, 23, 25, 27, 29 are gaps *within* the family — note PSL(2,9) ≅ A₆ ≅ Sp₄(2)′ and PSL(2,11) (the second-smallest gap, Foster 110/220 per RESEARCH_PROGRAM tables but no drawing in `individual_graphs/`).
- **A₂**: q = 2, 3 only.
- **B₂**: Sp₄(2) only.
- **G₂**: q = 2 only.
- **Alternating**: A₅, A₆, S₅, S₆ implicitly via the exceptional isomorphisms above; nothing for A₇+ as such.
- **Sporadic**: M₁₁ attempted, nothing rendered.

### Missing — ranked by feasibility × interest

1. **²A₂ = PSU(3,q), Hermitian unitals.** PSU(3,3) we have *implicitly*
   (= G₂(2)′) but no graph is drawn in its own right; the unital S(2,4,28)
   for q=3 is small and drawable. PSU(3,2)... is solvable (not simple) —
   worth a remark, not a graph. Natural first target: **PSU(3,3) on 28
   points** (Hermitian unital, rank 2).
2. **C₂ at odd q: PSp₄(3) ≅ PSU(4,2)** (order 25920, the smallest exceptional
   coincidence between a symplectic and a unitary group). Levi graph of
   GQ(3,3) has 40+40 = 80 vertices, valency 4 — beyond cubic but within the
   incidence-drawing pipeline (like `auto-pg23`). Also the unique simple
   group with an exceptional Schur multiplier story at 80 points.
3. **²B₂ = Sz(8)** (order 29120, smallest Suzuki). Acts on the Suzuki ovoid
   (65 points). Cubic Sabidussi search is a finite GAP computation —
   exactly the shape of `m11_cubic_search.g`. No generalized polygon
   (rank 1 BN-pair), so this tests the programme outside rank 2.
4. **²G₂(3)′ ≅ PSL(2,8)** — we already have Macbeath; drawing it would
   close the smallest Ree case "for free." The first genuine Ree group
   R(27) is order ~10¹⁰: out of reach.
5. **³D₄(2)** (order 211 341 312) — acts on the **twisted triality hexagon
   GH(2,8)**: 819 points, the only other hexagon with s=2. The natural
   sequel to the why-13 story (its Hecke/orbital structure vs ours). GAP
   can handle the 819-point permutation representation; the 12-cage
   analogue has 819+2457 vertices — drawable only as quotients.
6. **A₂ at q = 4: PSL(3,4)** — the projective plane PG(2,4), 21+21 Levi
   graph, central to the Mathieu groups (S(5,8,24) is built from it). A
   bridge from the Lie-type programme to the sporadics.
7. **Sporadics adjacent to what we have**:
   - **M₁₂** (and finish M₁₁): both have cubic arc-transitive candidates;
     M₁₂'s outer automorphism mirrors our PSU(3,3) chirality story.
   - **J₂ ⊂ G₂(4)**: *the* sporadic of the G₂ thread; acts rank-3 on 100
     points (SRG(100,36,14,12)). Drawing its orbital graph extends §4.
   - **M₂₄ / M₂₂**: Steiner systems S(5,8,24), S(3,6,22); incidence
     drawings, not cubic.
   - **HS** (Higman–Sims graph SRG(100,22,0,6)) — iconic, drawable with
     the orbit-ring method.
8. **Higher-rank Aₙ: PSL(4,2) ≅ A₈** (order 20160) on the 15 points of
   PG(3,2) — ties the alternating family to Lie type; the 15+35 Levi
   graph of points/planes is small.
9. **F₄, E₆, E₇, E₈, ²F₄, ²E₆**: smallest cases (²F₄(2)′ Tits group, order
   17 971 200) are at the edge of feasibility; ²F₄(2)′ acts on the
   generalized octagon GO(2,4) — **the only thick octagon with s=2**, and
   Feit–Higman says octagons are the last family. A long-term jewel:
   it would complete a "generalized polygon trilogy" (Heawood 6-cage /
   GH(2,2) 12-cage / GO(2,4) octagon) — 3-gon, 6-gon, 8-gon at s=2.
10. **Cyclic/abelian shadows**: circulants are the degenerate base case —
    deliberately out of scope.

---

## 9. ⚑ Future Research Directions

Consolidated open work, in suggested order of attack.

### 9.1 Finish what's started (GAP done, nothing rendered)
- **Macbeath surface graph** (PSL(2,8) = ²G₂(3)′): run the drawing pipeline
  on the output of `src/census/macbeath_graph.g`. Bonus: closes the smallest
  Ree case for free.
- **Clebsch graph** (`hecke/hecke_clebsch.g` + `Clebsch.md`): render it, and
  fold its Hecke-algebra notes into the §4 story.
- **M₁₁ cubic search**: document the outcome of `m11_cubic_search.g`
  (did any (H, a) give a connected cubic graph?) and render the result or
  record the negative.
- **PSL(2,11)**: in the RESEARCH_PROGRAM tables (110/220 vertices) but
  absent from `individual_graphs/` — fill the smallest gap in the A₁ census.

### 9.2 New computations (gap analysis of §8, ranked)
1. **PSU(3,3) on the 28-point Hermitian unital** — our G₂(2)′ workhorse,
   never drawn in its own right.
2. **Sz(8)** cubic Sabidussi search (65-point ovoid action) — first test of
   the programme outside rank-2 BN-pairs.
3. **PSp₄(3) ≅ PSU(4,2)**: GQ(3,3) Levi graph (80 vertices, valency 4).
4. **PSL(3,4)**: PG(2,4) Levi graph — the bridge toward the Mathieu groups.
5. **J₂ ⊂ G₂(4)**: rank-3 orbital graph on 100 points — the sporadic
   continuation of the G₂ thread.
6. **M₁₂** (after M₁₁), exploiting its outer automorphism (mirrors the
   PSU(3,3) chirality story of `~/git/hecke/why13.tex` §6).
7. **³D₄(2)** on the twisted triality hexagon GH(2,8) (819 points): the
   sequel to why-13 — does its Hecke algebra have an analogous B² = c·I?
8. **²F₄(2)′ (Tits group) on the generalized octagon GO(2,4)** — the
   capstone: completes the Feit–Higman family n = 3, 4, 6, 8 at s = 2
   (Heawood, Tutte–Coxeter, Tutte 12-cage, octagon).

### 9.3 Housekeeping
- Verify the Aut attributions in §6 currently given only as group orders
  (Möbius–Kantor 96, Pappus 216, Nauru 144, F26A 78, F056A).
- Cross-link RESEARCH_PROGRAM.md "Beyond A₁" to `~/git/hecke/why13.tex`
  (its Key Question 5 is now answered).

---

*Provenance: assembled from `RESEARCH_PROGRAM.md`, directory listings of
`individual_graphs/`, `src/`, `src/census/`, `quotients/`, `bolza-surface/`,
`clayworth/`, and `~/git/hecke/`. Group attributions for the Foster-census
entries are standard; entries marked TBD have GAP output but no rendered
drawing yet.*
