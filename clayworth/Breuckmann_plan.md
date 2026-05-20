# Formalising ∂² = 0: Surface Tilings → CSS Quantum Codes in Lean 4

## Goal

Prove in Lean 4 that any 2-cell embedding of a graph on a closed surface gives a valid CSS quantum code, with k = 2g logical qubits. Currently we assume the CSS orthogonality condition (H_X · H_Zᵀ = 0) and the rank hypotheses as parameters. This project removes both assumptions.

## The mathematics

A 2-cell embedding of a graph on a closed orientable surface gives a chain complex:

```
C₂ --∂₂--> C₁ --∂₁--> C₀
```

where C₀ = F₂^V (vertices), C₁ = F₂^E (edges), C₂ = F₂^F (faces). The CSS code has:
- H_X = ∂₁ᵀ (vertex-edge incidence, V × E)
- H_Z = ∂₂ᵀ (face-edge incidence, F × E, transposed)
- CSS condition: H_X · H_Zᵀ = ∂₁ᵀ · ∂₂ = (∂₁ ∘ ∂₂)ᵀ = 0ᵀ = 0

The proof that ∂₁ ∘ ∂₂ = 0 over F₂: each column of ∂₂ represents a face boundary. Each entry (∂₁ ∘ ∂₂)_{v,f} counts the number of edges shared between vertex v and face f. In a 2-cell embedding, a vertex appears on a face boundary 0 or 2 times (it enters and exits). Over F₂, 2 = 0. Done.

The rank theorem: for a connected graph, rank(∂₁ᵀ) = V - 1. The all-ones vector is in the kernel (every vertex is incident to an even number of edges in a cycle — wait, that's not right for ∂₁). Actually: the kernel of ∂₁ᵀ over F₂ is spanned by the all-ones vector when the graph is connected, because ∂₁ᵀ · 1 = (degree of each vertex mod 2). For a surface tiling where every vertex has the same degree q, this is either 0 (q even) or the all-ones vector (q odd). The general argument: rank(B) over F₂ where B is the incidence matrix of a connected graph is V - 1. This holds over ANY field (the kernel is 1-dimensional, spanned by the all-ones vector over F₂ when the graph is bipartite, or trivial otherwise — actually this needs care).

**Correction**: Over F₂, the rank of the incidence matrix B of a connected graph is:
- V - 1 if the graph has at least one odd cycle (i.e., is not bipartite)
- V - 1 if bipartite (the all-ones vector on each part gives the kernel — wait, let me think again)

Actually: Over any field, the incidence matrix of a connected graph has rank V - 1. The kernel over ℝ is spanned by the all-ones vector. Over F₂, the situation is different because B over F₂ is the same as the adjacency-to-edge map mod 2. The standard proof (over ℝ) uses the fact that deleting one row gives a full-rank matrix, which works over any field. So rank = V - 1 over F₂ as well. The same argument applies to F - 1 for the dual.

## What Mathlib provides

| Component | Mathlib status |
|-----------|---------------|
| `ZMod 2` as a field | ✅ Ready |
| `Matrix.rank` over fields | ✅ Ready |
| Rank-nullity theorem | ✅ Ready |
| `SimpleGraph` with `incMatrix` | ✅ Ready |
| `SimpleGraph.Connected` | ✅ Ready |
| Chain complexes (`HomologicalComplex`) | ✅ Ready (but categorical, not matrix-based) |
| CW-complex / 2-cell embedding | ❌ Missing |
| Incidence matrix rank = V-1 over F₂ | ❌ Missing |
| Face-edge incidence from embedding | ❌ Missing |

## Architecture

### Phase 1: CellularSurface (the scaffolding)

Define a purely combinatorial structure — no topology, just the data that a topologist would extract from a 2-cell embedding.

```lean
structure CellularSurface where
  V : ℕ                -- number of vertices (0-cells)
  E : ℕ                -- number of edges (1-cells)
  F : ℕ                -- number of faces (2-cells)
  -- Incidence data
  edge_src : Fin E → Fin V         -- source vertex of each edge
  edge_tgt : Fin E → Fin V         -- target vertex of each edge
  edge_irrefl : ∀ e, edge_src e ≠ edge_tgt e   -- no loops
  -- Face boundary: each face is a cycle of edges
  face_len : Fin F → ℕ             -- number of edges bounding each face
  face_edge : (f : Fin F) → Fin (face_len f) → Fin E   -- edges in order
  -- The key axiom: every edge borders exactly two faces
  two_faces : ∀ e : Fin E, (Finset.univ.filter (fun f => 
    ∃ i, face_edge f i = e)).card = 2
  -- Boundary consistency: face edges form a closed walk
  face_closed : ∀ (f : Fin F) (i : Fin (face_len f)),
    edge_tgt (face_edge f i) = edge_src (face_edge f ⟨(i + 1) % face_len f, ...⟩)
```

**Design decision**: faces are ordered sequences of edges forming closed walks, not just sets. This makes the "vertex appears 0 or 2 times" argument concrete.

### Phase 2: Boundary operators as matrices

```lean
-- ∂₁: the vertex-edge incidence matrix (V × E over F₂)
def CellularSurface.d1 (S : CellularSurface) : Matrix (Fin S.V) (Fin S.E) (ZMod 2) :=
  fun v e => if edge_src e = v then 1 else 0 + if edge_tgt e = v then 1 else 0

-- ∂₂: the edge-face boundary matrix (E × F over F₂)
def CellularSurface.d2 (S : CellularSurface) : Matrix (Fin S.E) (Fin S.F) (ZMod 2) :=
  fun e f => if (∃ i, face_edge f i = e) then 1 else 0
  -- Over F₂ orientation doesn't matter: we just count incidences mod 2
```

**Key simplification over F₂**: we don't need to orient edges or faces. The incidence matrices are just 0/1 matrices counting membership mod 2. This is much simpler than the ℤ chain complex in `bolza_cohomology.py`.

### Phase 3: ∂² = 0 (the satisfying proof)

```lean
theorem d1_comp_d2_eq_zero (S : CellularSurface) : S.d1 * S.d2 = 0
```

Proof sketch:
- (d1 * d2)_{v,f} = Σ_e d1_{v,e} · d2_{e,f}
- This counts how many boundary edges of face f are incident to vertex v, mod 2
- The face boundary is a closed walk. Vertex v can appear on this walk:
  - 0 times → contributes 0
  - k times → each appearance contributes 2 incident edges (entering + leaving), so total = 2k ≡ 0 mod 2
- Therefore (d1 * d2)_{v,f} = 0 for all v, f

The formal proof needs to argue about the closed walk structure of face boundaries. The `face_closed` axiom gives us the tools.

### Phase 4: rank(incidence matrix) = V - 1 for connected graphs

```lean
theorem incMatrix_rank_connected (S : CellularSurface) 
    (hconn : connected S) : S.d1.rank = S.V - 1
```

Proof strategy (works over any field, including F₂):
1. Show that the rows of d1 span a space of dimension ≤ V (trivially)
2. Show that any V-1 rows are linearly independent: take a spanning tree T of the graph (exists by `hconn`). The submatrix of d1 restricted to the V-1 tree edges is upper-triangular (after reordering by BFS from root) with 1s on the diagonal. Upper-triangular with nonzero diagonal → full rank.
3. Show that rank < V: the sum of all rows (over F₂) is zero iff all vertices have even degree, or equals something nonzero otherwise. Actually: rank ≤ V - 1 because the image of d1ᵀ (as a linear map from F₂^E → F₂^V) lies in the hyperplane {x : Σ xᵢ · deg(i) mod 2 = 0}... hmm, this needs more care.

**Alternative approach**: Show rank(d1) = V - 1 directly via spanning trees. A spanning tree of a connected graph on V vertices has V - 1 edges. The submatrix of the incidence matrix restricted to tree edges has rank V - 1 (each tree edge connects a new vertex to the existing tree — this is a triangular argument). Adding non-tree edges doesn't increase rank beyond V - 1 because... actually it could. Over F₂ we need: the all-ones vector is NOT in the row space (it has weight V, and each row has weight 2, so over F₂ the all-ones vector is in the row space iff V is even... no, that's not right either).

**Correct argument**: The incidence matrix B has B · 1 = (deg(v₁), ..., deg(vₙ))ᵀ. Over ℝ, each row has entries summing to 0 (one +1 and one -1), so rank ≤ V-1. Over F₂, each row sums to 0 (one 1 and one 1, sum = 0 mod 2), so the all-ones vector 1_E is in the kernel of Bᵀ... no, that's the wrong direction.

Let me restart. B is V × E. Each column has exactly two 1s (the two endpoints). So Bᵀ · 1_V = (2, 2, ..., 2)ᵀ = 0 over F₂. Therefore 1_V ∈ ker(Bᵀ), so rank(B) = rank(Bᵀ) ≤ V - 1. Equality follows from the spanning tree argument.

**This is clean and provable.** The key lemma: each column of the incidence matrix has exactly two 1s → the all-ones vector is in the left kernel → rank ≤ V - 1. Then spanning tree → rank ≥ V - 1.

### Phase 5: Assembling the CSS code

With Phases 3 and 4, we can construct a `CSSFromTiling` from a `CellularSurface` with both the orthogonality condition and rank hypotheses proved. The existing `css_k_eq_2g` theorem then gives k = 2g with no sorry anywhere in the chain.

## File structure

```
clayworth/CSSCode/
  Basic.lean              -- BinaryLinearCode (done)
  CSSCode.lean            -- SurfaceTiling, CSSFromTiling, css_k_eq_2g (done)
  Clayworth.lean          -- Concrete instances (done)
  CellularSurface.lean    -- Phase 1: CellularSurface structure
  BoundaryOperators.lean  -- Phase 2: d1, d2 as matrices over F₂
  BoundarySquaredZero.lean -- Phase 3: d1 * d2 = 0
  IncMatrixRank.lean      -- Phase 4: rank = V - 1 for connected graphs
  SurfaceCSS.lean         -- Phase 5: CellularSurface → CSSFromTiling
```

## Difficulty assessment

| Phase | Difficulty | Estimated effort | Sorry-free? |
|-------|-----------|-----------------|-------------|
| 1. CellularSurface | Medium | Getting the axioms right is the design challenge. Multiple iterations likely. | N/A (definitions) |
| 2. Boundary operators | Easy | Straightforward matrix construction from incidence data. | Yes |
| 3. ∂² = 0 | Medium | The closed-walk argument needs careful manipulation of Fin arithmetic and Finset sums. The *mathematics* is trivial but the *formalisation* involves modular index arithmetic in face boundaries. | Yes — this is the target |
| 4. rank = V - 1 | Hard | Spanning trees in Lean, triangular matrix argument, or the "each column has two 1s" + "connected → rank ≥ V-1" approach. Spanning trees may not be in Mathlib. | Likely needs sorry initially |
| 5. Assembly | Easy | Just plugging phases 3-4 into existing theorems. | Depends on 3 and 4 |

## Critical risks

1. **CellularSurface axioms too strong or too weak.** If the structure is too rigid, concrete instances (Bolza, Klein, Clayworth) won't fit. If too loose, ∂² = 0 won't follow. Need to iterate.

2. **Fin arithmetic hell.** Face boundaries involve cyclic indexing (`(i+1) % n`), which is notoriously painful in Lean. May need helper lemmas about `Fin` and `ZMod`.

3. **Spanning trees.** Mathlib has `SimpleGraph.IsTree` and `SimpleGraph.IsSpanning` but may not have "every connected graph has a spanning tree" as a constructive term. May need to build this.

## What this achieves

If completed, this would be:
- The **first formalisation** of the surface code construction (∂² = 0) from combinatorial topology in any proof assistant
- A **Mathlib-quality contribution**: the incidence matrix rank theorem (Phase 4) is independently useful
- A verified proof that the Clayworth, Klein quartic, and Bolza CSS codes are valid — not just that their parameters are correct, but that the orthogonality condition actually holds

## References

- Breuckmann & Terhal, "Constructions and Noise Threshold of Hyperbolic Surface Codes", IEEE TIT 2016, arXiv:1506.04029
- QECLean (Stavan Jain), github.com/Stavan-Jain/QECLean — existing CSS formalisation, toric code only
- Mathlib: `Combinatorics.SimpleGraph.IncMatrix`, `LinearAlgebra.Matrix.Rank`
