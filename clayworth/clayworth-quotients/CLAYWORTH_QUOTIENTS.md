# Clayworth Graph: Block System Quotients

The Clayworth graph is Sab(G₂(2), C₃) — 4032 vertices, cubic, on a genus-505 surface. The C₃ subgroup belongs to conjugacy class 1 (28 conjugates, fixes 0 points of GH(2,2)).

## The Langer graph is NOT a quotient

The Langer graph (63 vertices, 6-regular) lives on the 63 points of GH(2,2). It is Sab(G₂(2), H₁₉₂) where H₁₉₂ is the point stabiliser (order 192). For the coset quotient map G₂(2)/C₃ → G₂(2)/H₁₉₂ to be well-defined, we need C₃ ≤ H₁₉₂. But the Clayworth C₃ fixes **zero** points of GH(2,2), so it is not contained in any point stabiliser. The two graphs are algebraic siblings (same parent group G₂(2)) but not related by quotient.

There are 9 subgroups K₁₉₂ ≅ (SL(2,3) ⋊ C₄) ⋊ C₂ of order 192 containing C₃. These give 63-block quotients, but ALL 6048 edges are internal to the 64-element blocks — the quotient graph has degree 0. No 63-vertex quotient of Clayworth has any edges.

## Conjugacy classes of C₃ in G₂(2)

| Class | Conjugates | Fixed points on GH(2,2) | Clayworth connected? |
|-------|------------|-------------------------|---------------------|
| 1 | 28 | 0 | **Yes** (4032 vertices) |
| 2 | 336 | 3 | No (56-vertex component = Klein graph) |

Only class 1 gives a connected Clayworth graph.

## Intermediate subgroup lattice

G₂(2) has 252 subgroups between C₃ and G₂(2). These give block systems on the 4032 Clayworth vertices. Below are the distinct regular quotients (one representative per distinct quotient type):

| |K| | K structure | Blocks | Block size | Degree | Edges | Internal edges |
|-----|-------------|--------|------------|--------|-------|----------------|
| 6 | C₆ | 2016 | 2 | 3 | 3024 | 0 |
| 6 | C₆ | 2016 | 2 | 6 | 6048 | 0 |
| 9 | C₃ × C₃ | 1344 | 3 | 9 | 6048 | 0 |
| 12 | C₁₂ | 1008 | 4 | 3 | 1512 | 0 |
| 12 | C₁₂ | 1008 | 4 | 6 | 3024 | 0 |
| 12 | C₁₂ | 1008 | 4 | 12 | 6048 | 0 |
| 18 | C₃ × C₆ | 672 | 6 | 9 | 3024 | 0 |
| 18 | C₃ × C₆ | 672 | 6 | 18 | 6048 | 0 |
| 24 | SL(2,3) | 504 | 8 | 1 | 252 | 0 |
| 24 | SL(2,3) | 504 | 8 | 3 | 756 | 0 |
| 24 | SL(2,3) | 504 | 8 | 8 | 2016 | 0 |
| 24 | SL(2,3) | 504 | 8 | 24 | 6048 | 0 |
| 27 | C₃ × C₃ × C₃ | 448 | 9 | 27 | 6048 | 0 |
| 36 | — | 336 | 12 | 18 | 3024 | 0 |
| 36 | — | 336 | 12 | 36 | 6048 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 0 | 0 | 6048 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 1 | 126 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 2 | 252 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 3 | 378 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 8 | 1008 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 16 | 2016 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 24 | 3024 | 0 |
| 48 | SL(2,3) ⋊ C₂ | 252 | 16 | 48 | 6048 | 0 |
| 54 | (C₃ × C₃ × C₃) ⋊ C₂ | 224 | 18 | 27 | 3024 | 0 |
| 54 | (C₃ × C₃ × C₃) ⋊ C₂ | 224 | 18 | 54 | 6048 | 0 |
| 96 | SL(2,3) ⋊ C₄ | 126 | 32 | 32 | 2016 | 0 |
| 108 | ((C₃×C₃) ⋊ C₃) ⋊ C₄ | 112 | 36 | 27 | 1512 | 0 |
| 192 | (SL(2,3) ⋊ C₄) ⋊ C₂ | 63 | 64 | 0 | 0 | 6048 |
| 216 | ((C₃×C₃) ⋊ C₃) ⋊ C₈ | 56 | 72 | 27 | 756 | 0 |
| 432 | ((C₃×C₃) ⋊ C₃) ⋊ QD₁₆ | 28 | 144 | 27 | 378 | 0 |
| 6048 | PSU(3,3) = G₂(2)' | 2 | 2016 | 1 | 1 | 0 |

## Small quotients are trivial

The quotients with ≤ 56 vertices are all degenerate:

- **28v, deg 27**: K₂₈ — the complete graph. |Aut| = 28!. Every block connects to every other.
- **56v, deg 27**: 2 × K₂₈ — disjoint union of two complete graphs. |Aut| = 2 × (28!)². The K₂₁₆ subgroup splits each K₂₈-block in half, but inter-block edges remain complete.
- **63v, deg 0**: edgeless — all 6048 edges lie within the 64-element blocks.
- **2v, deg 1**: K₂ via PSU(3,3) = G₂(2)' — the Clayworth graph is bipartite under this block system.

This is a consequence of the C₃ fixing zero points of GH(2,2): the quotient lattice jumps from trivially complete to large sparse structures with nothing small and sparse in between.

## The interesting quotients are large

The first non-trivial sparse quotients appear at **252 vertices** (|K| = 48, K ≅ SL(2,3) ⋊ C₂):

| Blocks | Degree | Edges | Notes |
|--------|--------|-------|-------|
| 252 | 1 | 126 | perfect matching |
| 252 | 2 | 252 | 2-regular (union of cycles) |
| 252 | 3 | 378 | cubic — the first interesting quotient |
| 252 | 8 | 1008 | |
| 252 | 16 | 2016 | |
| 252 | 24 | 3024 | |
| 252 | 48 | 6048 | |

However, the **252-vertex cubic** quotient turns out to be **63 disjoint copies of K₄** (disconnected, girth 3, |Aut| ≈ 63! × 24⁶³). Similarly, the degree-1 quotient is 126 × K₂ (perfect matching) and the degree-2 quotient is a union of cycles.

**The Clayworth graph has no interesting small quotients.** Every quotient is either disconnected (disjoint complete graphs) or trivially dense (complete graphs). This is because the C₃ from conjugacy class 1 sits in a "generic" position in G₂(2) — it doesn't align with any of the natural geometric structures (points, lines, flags of GH(2,2)). The subgroup lattice above C₃ consists of groups that either absorb all edges internally (edgeless quotients) or spread edges uniformly (complete graphs).

| Quotient | Actual structure |
|----------|-----------------|
| 28v, deg 27 | K₂₈ |
| 56v, deg 27 | 2 × K₂₈ |
| 63v, deg 0 | 63 isolated vertices |
| 252v, deg 3 | 63 × K₄ |
| 252v, deg 1 | 126 × K₂ |
| 2v, deg 1 | K₂ |

### The Tutte 12-cage is NOT a quotient
The 126-block quotient (|K| = 96) has degree 32, not 3. The Tutte 12-cage does not appear as any quotient of the Clayworth graph.

### K₂ — the bipartition (|K| = 6048)
The simple subgroup G₂(2)' = PSU(3,3) (order 6048, index 2) gives a 2-block partition. The quotient is K₂ — the Clayworth graph is bipartite, with 2016 vertices in each part and all 6048 edges crossing the cut.
