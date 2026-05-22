# Clayworth Graph: Quotient Classification

The Clayworth graph is Sab(G₂(2), C₃) — 4032 vertices, cubic, arc-transitive, embeddable on a genus-505 surface. The C₃ subgroup belongs to the conjugacy class with 0 fixed points on GH(2,2) (28 conjugates).

## The Langer graph is NOT a quotient

**Theorem.** The Langer graph (63 vertices, 6-regular, primitive) is not a quotient of any Clayworth graph.

**Proof.** The Langer graph is Sab(G₂(2), H₁₉₂) where H₁₉₂ = Stab_G(p) is the point stabiliser of GH(2,2), with structure (SL(2,3) : C₄) : C₂. For a coset projection G₂(2)/C₃ → G₂(2)/H₁₉₂^g to exist, we need C₃ ≤ H₁₉₂^g.

Every conjugate of H₁₉₂ is a point stabiliser Stab_G(p) for some p ∈ {1,...,63}. A subgroup C₃ ≤ Stab_G(p) requires the C₃ generator to fix p. The Clayworth C₃ generator is a product of 21 disjoint 3-cycles with **zero fixed points** on the 63 GH(2,2) points. Therefore C₃ ⊄ Stab_G(p) for any p.

Furthermore, the three order-192 overgroups of C₃ in G₂(2) have structure (((C₄×C₄):C₂):C₃):C₂, which is NOT conjugate to H₁₉₂ ≅ (SL(2,3):C₄):C₂. The resulting 63-vertex block quotient is 24-regular and verified NOT isomorphic to the 6-regular Langer graph. □

**Contrast:** The OTHER C₃ conjugacy class (336 conjugates, 9 fixed points) IS contained in point stabilisers. Its coset graph is disconnected — the 56-vertex component is the Klein graph.

## Conjugacy classes of C₃ in G₂(2)

| Class | Conjugates | Fixed points | Coset graph |
|-------|------------|-------------|-------------|
| 1 | 28 | 9 | Disconnected (Klein component, 56v) |
| 2 | 336 | **0** | **Clayworth** (connected, 4032v) |

## The quotient lattice

G₂(2) has **85 subgroups** between C₃ and G₂(2). These define block systems on the 4032 Clayworth vertices. Quotient structure depends on the choice of Clayworth orbital (= choice of double coset / connection set).

**Only 2 quotients are degenerate** (for any orbital choice):
- K₂₈ via K = ((C₃×C₃):C₃):QD₁₆ (order 432) — always complete
- K₂ via K = PSU(3,3) (order 6048) — always bipartite

**The other 83 quotients are connected non-complete vertex-transitive graphs.**

## Star quotients

The most notable quotients (sparse, high girth, full G₂(2) symmetry):

| V | Degree | Girth | |Aut| | Edge mult | K | |K| | Block size |
|-----|--------|-------|-----------|-----------|-----|------|------------|
| **2016** | **6** | **8** | ? | 1 | C₆ | 6 | 2 |
| 2016 | 6 | 7 | ? | 1 | C₆ | 6 | 2 |
| **1344** | **9** | **6** | ? | 1 | C₃×C₃ | 9 | 3 |
| 1008 | 12 | 6 | ? | 1 | A₄ | 12 | 4 |
| **576** | **7** | **6** | G₂(2) | 3 | C₇:C₃ | 21 | 7 |
| **504** | **6** | **6** | G₂(2) | 4 | S₄ | 24 | 8 |
| **72** | **7** | **4** | G₂(2) | 24 | PSL(3,2) | 168 | 56 |
| 126 | 16 | 4 | G₂(2) | 6 | ((C₄×C₄):C₃):C₂ | 96 | 32 |

### The 2016-vertex girth-8 graph

The block-size-2 quotient via K ≅ C₆ produces a **6-regular graph on 2016 vertices with girth 8** and edge multiplicity 1 (every inter-block edge is simple). All 6048 edges cross between blocks; no internal edges exist. The Moore bound for (6,8) is 187, so this is not a cage, but it is a natural algebraic construction.

The seven K=6 overgroups give 2016-vertex quotients with girth spectrum: **[3, 4, 4, 6, 7, 7, 8]**. The C₆ overgroups give higher girth (7, 7, 8) than the S₃ ones (3, 4, 4, 6).

### The 72-vertex quotient via PSL(3,2)

G₂(2) contains PSL(3,2) = PSL(2,7) as a maximal subgroup (order 168, index 72). The 72-block quotient is a 7-regular graph with girth 4 and Aut = G₂(2). Each pair of adjacent blocks shares 24 crossing edges. This is a G₂(2)-orbital graph on a representation NOT isomorphic to the GH(2,2) point action.

### The 504-vertex girth-6 graph

Via K = S₄ (order 24, blocks of 8): a 6-regular graph on 504 vertices with girth 6 and Aut = G₂(2). Edge multiplicity 4 (each adjacent block pair shares 4 edges). This is the **sparsest quotient with full G₂(2) symmetry**.

## Full classification

All 85 quotients listed by block size (density = deg/(V-1)):

| |K| | K structure | V | Deg | Girth | Dens% | |Aut| |
|-----|-------------|------|-----|-------|-------|-------|
| 6 | C₆ | 2016 | 6 | 8 | 0.3 | ? |
| 6 | C₆ | 2016 | 6 | 7 | 0.3 | ? |
| 6 | C₆ | 2016 | 6 | 7 | 0.3 | ? |
| 6 | S₃ | 2016 | 6 | 3 | 0.3 | ? |
| 6 | S₃ | 2016 | 6 | 4 | 0.3 | ? |
| 6 | S₃ | 2016 | 6 | 4 | 0.3 | ? |
| 6 | S₃ | 2016 | 6 | 6 | 0.3 | ? |
| 9 | C₃×C₃ | 1344 | 9 | 6 | 0.7 | ? |
| 12 | D₁₂ | 1008 | 12 | 3–4 | 1.2 | ? |
| 12 | A₄ | 1008 | 12 | 4–6 | 1.2 | ? |
| 18 | C₃×S₃ | 672 | 9–18 | 3–4 | 1.3–2.7 | G₂(2) |
| 21 | C₇:C₃ | 576 | 7–21 | 4–6 | 1.2–3.7 | G₂(2) |
| 24 | S₄ / C₂×A₄ | 504 | 6–24 | 3–6 | 1.2–4.8 | G₂(2) |
| 27 | (C₃×C₃):C₃ | 448 | 27 | 4 | 6.0 | 193536 |
| 36 | S₃×S₃ | 336 | 9–36 | 3 | 2.7–10.7 | G₂(2) |
| 42 | C₇:C₆ | 288 | 14–42 | 3–4 | 4.9–14.6 | G₂(2) |
| 48 | C₂×S₄ | 252 | 12–48 | 3–4 | 4.8–19.1 | G₂(2) |
| 48 | (C₄×C₄):C₃ | 252 | 16–48 | 4 | 6.4–19.1 | 2×G₂(2) |
| 54 | various | 224 | 27–54 | 3–4 | 12.1–24.2 | varies |
| 96 | various | 126 | 16–48 | 3–4 | 12.8–38.4 | varies |
| 108 | various | 112 | 27–54 | 3–4 | 24.3–48.6 | varies |
| 168 | PSL(3,2) | 72 | 7–21 | 4 | 9.9–29.6 | G₂(2) |
| 192 | various | 63 | 24–32 | 3 | 38.7–51.6 | varies |
| 216 | various | 56 | 27 | 3–4 | 49.1 | varies |
| 336 | PSL(3,2):C₂ | 36 | 14–21 | 3 | 40–60 | G₂(2) |
| 432 | ((C₃×C₃):C₃):QD₁₆ | 28 | 27 | — | 100 | K₂₈ |
| 6048 | PSU(3,3) | 2 | 1 | — | 100 | K₂ |

## Orbital dependence

Different self-paired connected orbitals give non-isomorphic Clayworth graphs with different quotient structures. The block system depends only on K (not the orbital), but the quotient graph's edge set — and hence its degree and girth — depends on which orbital defines the Clayworth edges.

There are 1340 self-paired cubic suborbits of G₂(2) on G₂(2)/C₃. Each connected one gives a cubic arc-transitive graph on 4032 vertices. The specific surface-code Clayworth graph (from Breuckmann-Terhal) corresponds to a particular orbital, not yet identified in this analysis.

## Verified by

`prove_clayworth_quotients.g` — GAP computation using GRAPE package, EdgeOrbitsGraph for correct orbital construction, ConnectedComponents for quotient analysis.
