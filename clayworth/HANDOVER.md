# Handover: CSS Surface Codes in Lean 4

## What this is

A Lean 4 formalisation of CSS quantum error-correcting codes from surface tilings, following Breuckmann & Terhal (arXiv:1506.04029). The pipeline: combinatorial surface data → chain complex ∂² = 0 → rank theorems → k = 2g logical qubits.

## Current state

**Builds clean. 1 sorry (in a concrete instance, not a theorem).**

```bash
cd /Users/robin/git/graphviz/clayworth
lake build  # ~30s, 1 sorry warning from KleinQuartic.lean
```

### Files

| File | Purpose | Sorry |
|------|---------|-------|
| `CSSCode/Basic.lean` | `BinaryLinearCode` over F₂ | 0 |
| `CSSCode/CSSCode.lean` | `SurfaceTiling`, `CSSFromTiling`, `css_k_eq_2g` | 0 |
| `CSSCode/CellularSurface.lean` | `CellularSurface` structure, ∂₁, ∂₂, **∂₁∘∂₂ = 0** | 0 |
| `CSSCode/IncMatrixRank.lean` | rank(d1)=V-1, rank(d2)=F-1, SimpleGraph bridge | 0 |
| `CSSCode/SurfaceCSS.lean` | Assembly: `CellularSurface` → `CSSFromTiling` → k=2g | 0 |
| `CSSCode/Clayworth.lean` | Parameter-only instances (Clayworth, Klein, Heawood) | 0 |
| `CSSCode/Heawood.lean` | Concrete `CellularSurface` (14V, 21E, 7F, genus 1) | 0 |
| `CSSCode/Bolza.lean` | Concrete `CellularSurface` (16V, 24E, 6F, genus 2) | 0 |
| `CSSCode/KleinQuartic.lean` | Concrete `CellularSurface` (56V, 84E, 24F, genus 3) | **3** |
| `cellular_surfaces_data.json` | All embeddings (Bolza, Klein, Clayworth) as JSON | — |
| `Breuckmann_plan.md` | Original 5-phase plan (historical) | — |

### The 3 sorry's in KleinQuartic.lean

All in the `kleinQuarticSurface` definition:
- `edge_ne` — each edge has distinct endpoints
- `face_inj` — face boundaries are trails (no repeated edges)
- `face_closed` — face boundaries are closed walks

These are **data verification sorry's**, not theorem sorry's. The data was computed by GAP using the (2,3,7) triangle group and verified by Python. The sorry's exist because `native_decide` on 84-element `![...]` vectors takes >20 minutes and >1GB RAM. Solutions:
1. Wait for `native_decide` (works, just slow)
2. Use `Array`-backed lookup instead of `![...]`
3. Write a custom `Decidable` instance

## Key proof ideas

### ∂² = 0 (CellularSurface.lean:162)
Each face boundary is a closed walk. For vertex v and face f, the number of incident edges is 2 × (visits to v). Over F₂, 2k = 0. Formal proof: define `stepStart`/`stepEnd`, show `[src=v]+[tgt=v] = [start=v]+[end=v]` by `add_comm`, reindex by cyclic bijection `nextIdx`, conclude `x + x = 0` by `decide`.

### rank(d1) = V-1 (IncMatrixRank.lean:190)
**Upper bound** (d1_rank_le): Each column of d1 has two 1s → 1_V ∈ ker(d1ᵀ) → rank-nullity gives rank ≤ V-1.

**Lower bound** (d1_rank_ge): ker(d1ᵀ) = span{1_V} for connected graphs. Proof: x ∈ ker means x(src e) = x(tgt e) for all edges (from `d1T_mulVec_entry` + char 2). Walks propagate this (`walk_preserves_ker`). Connectivity makes x globally constant → x ∈ span{1_V}. So dim(ker) = 1, rank = V-1.

### rank(d2) = F-1 (IncMatrixRank.lean:310)
Same structure as d1 but for the dual graph. Key helper: `d2_entry` shows d2(e,f) is an indicator function (from `face_inj`). The `two_sides` hypothesis (exactly 2 faces per edge, as ℕ) isolates y(f₁)+y(f₂)=0 for dual-adjacent faces. Dual connectivity → y constant → ker = span{1_F}.

**Important**: the hypothesis was strengthened from `d2 *ᵥ 1 = 0` (mod 2) to `∀ e, card{f | d2 e f ≠ 0} = 2` (exactly 2 as ℕ). The weaker form is derived as `d2_mulVec_one_of_two_sides`.

## What's left to do

### 1. Create the PR
Reference arXiv:1506.04029 (Breuckmann & Terhal) and arXiv:1306.4798 (Robin's thesis). The PR description should note that all general theorems are sorry-free. **Include the images from `named_graphs/`** in the PR body — they show the symmetry-aware drawings of all the named graphs (Heawood, Möbius-Kantor, Klein, Zhou, Tutte 12-cage). Upload them inline in the GitHub PR description.

**PR body should include this paragraph about the Zhou and Langer graphs:**

> The named graphs folder also includes the Zhou graph (F182A, Sab(PSL(2,13), S₃), 182 vertices) and the Langer graph (63 vertices, 6-regular, quotient of the Clayworth graph by GH(2,2)). The Zhou graph illustrates the boundary of the construction: its vertex stabiliser S₃ is not cyclic, so it does not arise from a triangle group and has no known surface tiling. The Langer graph raises the question of how CSS code structure behaves under graph quotients.

### 2. Remove Klein quartic sorry's
Options:
- **Patience**: `native_decide` works but takes 20+ minutes. Try `set_option maxHeartbeats 4000000`.
- **Array representation**: Replace `![...]` with `fun i => Array.get! data i` where `data` is a compiled array.
- **Split the decision**: Verify `edge_ne`, `face_inj`, `face_closed` separately with smaller lemmas.

### 3. Clayworth CellularSurface instance
The 4032-vertex data exists in `cellular_surfaces_data.json` (6048 edges, 1008 dodecagonal faces). The `![...]` notation crashes Lean (6048-deep term tree exceeds memory). Solutions:
- `Array`-backed `Fin n → α` using `Array.get` with a proof of size
- Code-generate a Lean file using `DecidableEq` on a serialized format
- Use `native_decide` on `ByteArray`-backed data

### 4. Sabidussi bridge (ambitious)
Robin's mathlib4 PRs:
- #39548: `CosetGraph` construction `Sab(G, H, D) : SimpleGraph (G ⧸ H)`
- #39550: Sabidussi representation theorem (every vertex-transitive graph is a coset graph)
- #39551: Lorimer's theorem + quotient theorem

The chain G₂(2) + C₃ → `SimpleGraph` is formalised by #39548. The gap is the rotation system: going from `SimpleGraph` to `CellularSurface` requires the embedding data (currently in JSON). This is the "Lie group → surface → code" vision from `project_lie_group_programme.md` in memory.

### 5. Mathlib contributions
Independently useful results:
- `d1_col_sum_eq_zero`: incidence matrix columns sum to 0 over F₂
- `ker_d1T_le_span_one`: kernel characterisation of connected graph incidence matrices
- `d2_entry`: face-edge boundary matrix is an indicator function (from trail condition)

## Architecture notes

- `ZMod 2` needs `attribute [local instance] ZMod.instField` (from `Mathlib.Algebra.Field.ZMod`) for `Matrix.rank_transpose` and `finrank_range_add_finrank_ker`.
- `SimpleGraph` loopless field uses `Std.Irrefl` (class wrapping), not a bare proposition. Construct with `⟨fun v h => ...⟩`.
- `Submodule.one_le_finrank_iff` requires `[StrongRankCondition R]` and `[Module.Finite R S]`.
- `Finset.sum_ite_eq'` doesn't match constant functions in `![...]` vectors. Use `Finset.sum_eq_single_of_mem` instead.
- `Equiv.sum_comp` needs the function argument explicit for unification when the goal has compound expressions.
- `omega` does NOT handle `%` (modular arithmetic). Eliminate via `Nat.mod_eq_of_lt` / `Nat.mod_self` + case split, then `omega` handles linear arithmetic.

## Named graphs

| Name | Foster | Vertices | Sabidussi |
|------|--------|----------|-----------|
| Heawood | F014A | 14 | Sab(G₄₂, C₃) |
| Möbius-Kantor | F016A | 16 | Sab(GL(2,3), C₃) |
| Klein | F056A | 56 | Sab(PSL(2,7), C₃) |
| Zhou | F182A | 182 | Sab(PSL(2,13), S₃) |
| Langer | — | 63, 6-reg | Quotient of Clayworth by GH(2,2) |
| Clayworth | — | 4032 | Sab(G₂(2), C₃) |

## Commits (oldest first)

1. `20c647c` — Initial: CellularSurface + ∂² = 0 + rank bounds + assembly
2. `950012e` — Prove rank(d1) = V-1 via kernel characterisation
3. `d93ddee` — Heawood instance (Fano plane on torus)
4. `702d832` — Prove rank(d2) = F-1 (project becomes sorry-free)
5. `d4149e2` — Bolza instance (Möbius-Kantor on genus 2)
6. `b33733d` — Klein quartic instance (PSL(2,7) on genus 3)
7. `39e9277` — Computed embedding data for all three surfaces (JSON)
8. `290ea66` — Foster census names and Sabidussi descriptions
9. `9fd635c` — HANDOVER.md
10. `d91df12` — Named graphs: Zhou (F182A), Langer (63v quotient), Klein (F056A)
