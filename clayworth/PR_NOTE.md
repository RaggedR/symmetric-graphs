# PR Plan: Three PRs

## PR 1: CSS Surface Codes

CellularSurface formalism, ∂²=0, rank bounds, k=2g theorem, and four surface instances.

| File | What it proves |
|------|----------------|
| `Basic.lean` | `BinaryLinearCode` over F₂ |
| `CSSCode.lean` | `SurfaceTiling`, `CSSFromTiling`, `css_k_eq_2g` |
| `CellularSurface.lean` | `CellularSurface` structure, ∂₁, ∂₂, **∂₁∘∂₂ = 0** |
| `IncMatrixRank.lean` | rank(d₁) = V−1, rank(d₂) = F−1 |
| `SurfaceCSS.lean` | Assembly: CellularSurface → CSSFromTiling → k = 2g |
| `Heawood.lean` | Heawood surface (14V, 21E, 7F, genus 1) |
| `Bolza.lean` | Bolza / Möbius-Kantor surface (16V, 24E, 6F, genus 2) |
| `KleinSurface.lean` | Klein quartic surface (56V, 84E, 24F, genus 3) |
| `ClayworthSurface.lean` | Clayworth surface (4032V, 6048E, 1008F, genus 505) |
| `Q62.lean` | Q(6,2) self-orthogonality, Langer graph via Zorn product |

## PR 2: Sabidussi Quotients and Primitivity (this PR)

Quotient graph constructions for arc-transitive graphs, proved regular quotients from block systems, and primitivity classification of named graphs.

### How quotients drive the drawing algorithm

The layout engine (`draw_graph.py`) uses quotients as a **layout skeleton**:

1. Compute Aut(G) via GAP, find a subgroup H whose orbits form a block system
2. The quotient G/H is a small graph — this becomes the **radial structure**: each orbit = one concentric ring
3. Cyclic symmetry within each orbit gives rotational symmetry
4. Crossing minimisation tweaks angular offsets between rings

A voltage graph on K₂ with Z₇ (Heawood) gives 2 rings of 7 vertices each. The quotient determines the topology of the layout; the fibres fill in the detail.

### Named Sabidussi graphs

All graphs in this project are **arc-transitive** (symmetric) — Aut(G) acts transitively on directed edges.

| Graph | Sabidussi | Primitive? | Quotient proved | Quotient by |
|-------|-----------|------------|-----------------|-------------|
| Dodecahedron (20v) | Sab(A₅, C₃) | No | → Petersen (10v, 3-reg) | Z₂ antipodal |
| Cube (8v) | Sab(S₄, C₃) | No | → K₄ (4v, 3-reg) | Z₂ antipodal |
| Heawood (14v) | Sab(G₄₂, C₃) | No | → K₂ (2v, 1-reg) | Z₇ fibre translation |
| Möbius-Kantor (16v) | Sab(GL(2,3), C₃) | No | → K₂ (2v, 1-reg) | Z₈ fibre translation |
| Petersen (10v) | Sab(S₅, S₂×S₃) | ? | — | — |
| Klein (56v) | Sab(PSL(2,7), C₃) | **Yes** | — | — |
| Zhou (182v) | Sab(PSL(2,13), S₃) | **Yes** | — | — |
| Langer (63v, 6-reg) | Sab(G₂(2), H₁₉₂) | **Yes** | — | — |
| Clayworth (4032v) | Sab(G₂(2), C₃) | No* | no interesting quotients | — |

\* The Clayworth graph has block systems but all quotients are degenerate (complete graphs, disjoint unions of complete graphs, or edgeless). See `clayworth-quotients/CLAYWORTH_QUOTIENTS.md`.

### Primitivity

A symmetric graph is **primitive** when Aut(G) acts on its vertices with no non-trivial block system — there is no way to partition vertices into equal blocks that every automorphism permutes. Primitive graphs cannot be expressed as non-trivial coverings of smaller graphs.

- **Klein** (F056A): Aut ≅ PGL(2,7), order 336. Primitive. GAP's `AllBlocks` confirms.
- **Zhou** (F182A): Aut ≅ PSL(2,13), order 1092. Primitive.
- **Langer**: Aut ≅ G₂(2), order 12096. Primitive. Arc-transitive (arc stabiliser order 32, vertex stabiliser order 192).

### Generalized polygon hierarchy over F₂

| Geometry | Type | Incidence graph | Girth | Aut | Polarity? | Code |
|----------|------|-----------------|-------|-----|-----------|------|
| PG(2,2) (Fano) | gen. 3-gon | Heawood (14v) | 6 | PGL(2,7) | **Yes** → arc-transitive | [7,4,3] Hamming |
| GH(2,2) | gen. 6-gon | Tutte 12-cage (126v) | 12 | G₂(2) | **No** → semisymmetric | [63,56,3] |

PSL(2,7) ≅ GL(3,2) (Aut of Fano) is a subgroup of G₂(2). The Fano plane is the "seed" geometry; the split Cayley hexagon is the exceptional geometry built from G₂.

### The Tutte 12-cage: semisymmetric, not Sabidussi

The Tutte 12-cage is **edge-transitive but not vertex-transitive** (semisymmetric). Aut ≅ G₂(2) (order 12096) preserves the bipartition into 63 points + 63 lines. No polarity exists for GH(2,2) over F₂ (verified by nauty and bliss).

Not a Sabidussi graph — it is the incidence graph of a rank-2 coset geometry (Tits building of type G₂), with vertices from two coset spaces G₂(2)/P₁ ∪ G₂(2)/P₂ joined by incidence.

**Proved in Lean** (`Tutte12Cage.lean`):
- 3-regular, 189 edges
- Bipartite (every edge connects a point to a line)
- `langer_eq_tutte12_distance2`: the Langer graph equals the distance-2 graph on the point orbit — connecting the algebraic definition (Zorn product on Q(6,2)) to the geometric definition (collinearity in GH(2,2))

### Lean files

| File | What it proves |
|------|----------------|
| `QuotientGraph.lean` | `SimpleGraph.quotientGraph` by surjection π : V → W |
| `CosetGraph.lean` | Sabidussi construction: `IsConnectionSet`, `cosetGraph`, well-definedness |
| `PetersenQuotient.lean` | Dodecahedron → Petersen (Z₂ antipodal), both 3-regular, edge counts |
| `CubeQuotient.lean` | Cube → K₄ (Z₂ antipodal), both 3-regular, edge counts |
| `NamedCosetGraphs.lean` | Heawood/MK as voltage graphs on K₂, fibre quotients to K₂ |
| `NamedGraphs.lean` | Klein, Zhou, Langer parameter stubs with primitivity notes |
| `Tutte12Cage.lean` | Tutte 12-cage: 3-regular, bipartite, Langer = distance-2 on points |

## The Langer graph is NOT a quotient of Clayworth

Previously assumed but **disproved** in this work. The Clayworth C₃ belongs to conjugacy class 1 in G₂(2) (28 conjugates, fixes **zero** points of GH(2,2)). The Langer graph lives on the 63 points of GH(2,2) as Sab(G₂(2), H₁₉₂). For the coset map gC₃ ↦ gH₁₉₂ to be well-defined, C₃ must be contained in H₁₉₂ — but it isn't, because C₃ fixes no points.

Furthermore, all 9 subgroups of order 192 containing C₃ give 63-block quotients where ALL 6048 edges are internal (degree-0 quotient). And every small quotient of Clayworth is degenerate: K₂₈, 2×K₂₈, 63×K₄, etc. Full analysis in `clayworth-quotients/CLAYWORTH_QUOTIENTS.md`.

The Clayworth graph and Langer graph are **algebraic siblings** (both arise from G₂(2)) but are not related by any graph quotient. The two conjugacy classes of C₃ in G₂(2) give qualitatively different behaviour:

| C₃ class | Conjugates | Fixed points | Clayworth | Quotient to Langer? |
|----------|------------|--------------|-----------|-------------------|
| Class 1 | 28 | 0 | connected (4032v) | **No** |
| Class 2 | 336 | 3 | disconnected (56v component = Klein) | N/A |

## References

- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
- Breuckmann & Terhal, *Constructions and noise threshold of hyperbolic surface codes*, arXiv:1506.04029
- Sabidussi, *Vertex-transitive graphs*, Monatsh. Math. 68 (1964)
