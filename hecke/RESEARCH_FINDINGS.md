# Research Findings: Symplectic Groups, Hecke Algebras, and the Schläfli Graph

## Date: 2026-06-05

## The Three Pieces

### 1. The Framework — `hecke.tex`

Robin's paper develops the Hecke algebra H(G,H) = End_G(C[G/H]) from scratch:
- Coset space G/H, linearized to C[G/H]
- Double cosets H\G/H partition G; each HgH defines a Sabidussi orbital graph
- Adjacency matrices {T_{HgH}} form a basis for End_G(C[G/H])
- Proof: equivariance cuts n² → n parameters; H-invariance cuts n → r = |H\G/H|
- Structure constants c_d count two-step paths through orbital graphs
- The algebra End_G(C[G/H]) with this multiplication is the Hecke algebra H(G,H)

### 2. The Encyclopedia — Wan Zhexian's "Geometry of Classical Groups"

**Book**: `geometry-finite-groups.pdf` (Science Press, 2002, 2nd ed.)
**Survey**: `chinese-book.pdf` (Discrete Math 174, 1997, pp. 365–381)

Wan's four central questions for a classical group G_n acting on subspaces of F_q^(n):
- **(i)** What are the orbits? → **Double cosets H\G/H** (H = subspace stabilizer)
- **(ii)** How many orbits? → **r = dim H(G,H)**, rank of the Hecke algebra
- **(iii)** Orbit lengths? → **Valencies** of orbital adjacency matrices (Anzahl formulas)
- **(iv)** Subspaces in one orbit contained in a given subspace? → **Structure constants c_d**

Key Anzahl results:
- GL_n(F_q): orbits classified by dimension, N(m,n) = Gaussian binomial [n choose m]_q
- Sp_{2v}(F_q): orbits classified by type (m,s), number of orbits = ½(v+1)(v+2), lengths by formula (10)
- U_n(F_q), O_n(F_q,S): analogous classification and formulas

Applications: association schemes (Section 5), authentication codes (Section 6), projective codes (Section 7), lattices (Section 8), representation of forms (Section 9).

### 3. The Graphs — `~/git/graphviz` + new results

Robin's graphviz project builds Sabidussi coset graphs Sab(G, H, HgH) and studies arc-transitivity, surface embeddings, CSS quantum codes.

---

## Key Connection: Wan's Anzahl = Hecke Structure Constants

The Anzahl theorems throughout Wan's Chapters 3–7 literally compute the multiplication table of H(G,H) for each classical group. Robin's `hecke.tex` gives the abstract framework; Wan fills in the explicit numbers.

The association schemes in Wan's Section 5 are the Bose-Mesner algebras = commutative Hecke algebras (for Gelfand pairs).

---

## PSL(2,q) ≅ SO(3,q) — The Orthogonal Geometry Connection

PSL(2,q) ≅ SO(3,q) via the adjoint representation on sl(2) (documented in `~/git/graphviz/CLASSIFICATION.md`). This means Robin's three families of cubic arc-transitive graphs:

| Family | H | |V| | Notes |
|--------|---|-----|-------|
| 1 | C₃ | q(q²−1)/6 | All primes q ≥ 5 |
| 2 | S₃ | q(q²−1)/12 | 2-arc-transitive, Aut = PGL(2,q) |
| 3 | S₄ | q(q²−1)/48 | q ≡ ±1 mod 8 only |

are **orbital graphs of the 3-dimensional orthogonal geometry over F_q**. The chain C₃ ⊂ S₃ ⊂ S₄ exhausts the possibilities because A₅ and larger subgroups fail the valency-3 condition.

---

## Computational Census: Sp(2n, q)

### Setup

Script: `sp_cubic_census.g`

Searched PSp(4,q) = Sp(4,q)/Center for q = 2, 3 with valencies 3–10. Used symmetric double cosets (a⁻¹ ∈ HaH) for undirected graphs. Checked all conjugacy classes of subgroups and all conjugacy class representatives as double coset generators.

### Key observation: must use projective groups

Sp(4,2) ≅ S₆ (not simple) gives only disconnected cubic graphs. PSp(4,2)' = A₆ ≅ PSL(2,9) (simple, order 360) is the correct object, analogous to using PSL not SL.

### Results: PSp(4,2)' = A₆

| H | |V| | Val | Girth | |Aut| | Aut |
|---|-----|-----|-------|-------|-----|
| S₄ | 15 | 6 | 3 | 720 | S₆ |
| A₅ | 6 | 5 | 3 | 720 | S₆ |

No connected cubic graphs. (A₆ ≅ PSL(2,9), so cubic graphs from PSL(2,9) use the PSL construction instead.)

### Results: Sp(4,2) = S₆ (selected)

| H | |V| | Val | Girth | |Aut| | Notes |
|---|-----|-----|-------|-------|-------|
| C₆ | 120 | 6 | 5 | 720 | S₆ |
| D₁₀ | 72 | 5 | 6 | 2880 | C₂ × ((A₆.C₂):C₂) |
| S₄ | 30 | 6 | 4 | 1440 | C₂ × S₆ |
| A₅ | 12 | 5 | 4 | 1440 | C₂ × S₆ |
| S₅ | 6 | 5 | 3 | 720 | S₆ |

Again no connected cubic graphs from S₆.

### Results: PSp(4,3) ≅ PΩ(5,3) (order 25920, simple)

**NO connected cubic graphs.** All cubic double cosets give disconnected graphs with components of 4 or 8 vertices.

**Connected graphs found:**

| H | |H| | |V| | Val | Girth | Diam | |Aut| | Aut | SRG? |
|---|-----|-----|-----|-------|------|-------|-----|------|
| **(C₂⁴):A₅** | **960** | **27** | **10** | **3** | **2** | **51840** | **O(5,3):C₂** | **(27,10,1,5) — SCHLÄFLI** |
| S₅ | 120 | 216 | 10 | 4 | 4 | 51840 | O(5,3):C₂ | No (λ=0, μ varies) |
| C₂ × S₄ | 48 | 540 | 6 | ? | ? | ? | ? | (needs further analysis) |

---

## The Schläfli Graph — Headline Result

**Sab(PSp(4,3), (C₂⁴):A₅, HaH) = Schläfli graph**, where a is any involution outside H.

### Properties
- **SRG(27, 10, 1, 5)**: strongly regular, 27 vertices, 10-regular, λ=1, μ=5
- **Girth 3, diameter 2**
- **Aut = O(5,3):C₂**, order 51840 = |W(E₆)|, the Weyl group of E₆
- **Geometric meaning**: The 27 vertices are the **27 lines on a cubic surface**. Two vertices are adjacent iff the corresponding lines are **skew** (disjoint).

### Hecke algebra interpretation
- H(PSp(4,3), (C₂⁴):A₅) has rank = |H\G/H| (number of double cosets)
- The Schläfli adjacency matrix T_{HaH} is one basis element
- The structure constants count paths through this orbital graph
- The full Hecke algebra encodes the complete incidence geometry of the 27 lines

### Connection to Wan
- PSp(4,3) ≅ PΩ(5,3): this is an **orthogonal group** in Wan's classification
- The 27 cosets PSp(4,3)/(C₂⁴):A₅ are an orbit of subspaces in the 5-dimensional orthogonal geometry over F₃
- The Anzahl theorems for O(5,F₃,S) govern the structure of this Hecke algebra

### Connection to E₆
- |Aut(Schläfli)| = 51840 = |W(E₆)| = |Sp(4,3)| = 2|PSp(4,3)|
- The 27 lines on a cubic surface correspond to the 27 weights of the fundamental representation of E₆
- The E₆ root system has 72 roots; 72 = 27 × 16/6 (combinatorial relation)
- This connects symplectic geometry over F₃ to exceptional Lie theory

---

## The 216-Vertex Graph

**Sab(PSp(4,3), S₅, HaH)**: 216 vertices, 10-regular, girth 4, diameter 4.

### Remarkable properties
- **Triangle-free** (λ = 0) despite being 10-regular — very rare
- **Not strongly regular** (μ varies between non-adjacent pairs)
- Aut = O(5,3):C₂ = W(E₆), same as Schläfli
- 216 = 6³ = |PSp(4,3)|/|S₅| = 25920/120

### Geometric meaning (tentative)
- S₅ ≤ PSp(4,3) is a subgroup of index 216
- May correspond to the 216 roots of E₆ (72 positive + 72 negative + 72 ???) — needs verification
- Or to tritangent planes / Steiner trihedra of the cubic surface

---

## Structural Observations

### The Rank Barrier Conjecture

| Group | Rank | Min connected valency found |
|-------|------|---------------------------|
| PSL(2,q) ≅ SO(3,q) | 1 | 3 (cubic) |
| G₂(2) | exceptional | 3 (Clayworth, 4032v) |
| PSp(4,3) ≅ PΩ(5,3) | 2 | 5 or 6 |

**Conjecture**: For G a classical group of BN-rank r over F_q, the minimum valency of a connected Sabidussi graph grows with r. Rank-1 groups (PSL) achieve cubic; rank-2 and above cannot.

**Why?** For a connected Sabidussi graph Sab(G, H, HaH), we need ⟨H, a⟩ = G. In rank-1 groups, any two non-commuting elements tend to generate the whole group (simplicity + rank 1 = "thin" structure). In higher-rank groups, H and a can be trapped in a proper parabolic subgroup without generating G.

### Why PSp(4,3) cubic graphs disconnect

All cubic double cosets give components of size 4 or 8. The component size divides |G| but is much smaller than |G/H|. This means ⟨H, a⟩ is a proper subgroup — the cubic double coset doesn't "spread" enough to reach all cosets.

---

## Open Questions

1. **Is the 216-vertex graph a known graph?** Check Brouwer's database of distance-regular and strongly regular graphs. The parameters (216, 10, 0, μ-variable) may match a known object.

2. **What is the 540-vertex, 6-regular graph?** From PSp(4,3)/(C₂ × S₄). Need to identify which conjugacy class of C₂ × S₄ gives connectivity and compute full invariants.

3. **Unitary groups**: Does PSU(3,q) produce cubic graphs? PSU(3,4) has order 62400 — tractable.

4. **Can the rank barrier be made precise?** What is the exact minimum valency for Sp(2n,q) as a function of n and q?

5. **Wan's Anzahl for O(5,3)**: Compute the full Hecke algebra H(PSp(4,3), (C₂⁴):A₅) explicitly — all basis elements, all structure constants. This gives the complete incidence calculus of the 27 lines.

6. **CSS codes from PSp(4,3)?** The Schläfli graph has girth 3, which is too small for good CSS codes. But the 216-vertex graph has girth 4 — does it tile a surface? If so, what genus, and what CSS parameters?

---

## Files

| File | Description |
|------|-------------|
| `hecke.tex` | Main paper: Hecke algebra of (G,H) |
| `hecke.pdf` | Compiled output |
| `geometry-finite-groups.pdf` | Wan Zhexian, 2nd ed. (2002) |
| `chinese-book.pdf` | Wan survey, Discrete Math 174 (1997) |
| `sp_cubic_census.g` | GAP script: census of graphs from Sp(2n,q) |
| `RESEARCH_FINDINGS.md` | This document |
