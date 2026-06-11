# Hecke Algebras of the Catalogue: H(G,H) for every graph

**Companion to:** `CATALOGUE.md` (the graphs), `RESEARCH_PROGRAM.md` (the Sabidussi programme),
`hecke/why13.pdf` (the worked G₂(2) example).
**Computed by:** `hecke/hecke_census.g` → `hecke_census.out` (GAP 4.15.1 + GRAPE, 2026-06-10).

---

## 0. Methodology

For a transitive pair (G, H) the **Hecke algebra** (centralizer algebra) is

```
H(G,H) = End_G(ℂ[G/H]) ≅ ⊕_χ M_{m_χ}(ℂ),    m_χ = ⟨ind_H^G 1, χ⟩.
```

Everything follows from the permutation character π = ind_H^G 1:

- **rank** dim H(G,H) = Σ m_χ² = number of H-orbits on G/H (suborbits);
- **irreps of H(G,H)** ↔ constituents χ of π, with dim = m_χ;
- **Gelfand pair** (H commutative) ⟺ π is multiplicity-free ⟺ all m_χ ≤ 1.

So no structure constants are needed: `PermutationCharacter(G,H)` + `Irr(G)` give the
full Wedderburn decomposition. For graphs without a known group model we build the graph
in GRAPE, take G = Aut(Γ), H = vertex stabilizer. Suborbit lengths are the valencies of
the orbital graphs (the standard basis T₁,…,T_r of H).

**Reading the tables.** "Wedderburn" lists blocks as `d^(m)` meaning a constituent of
degree d with multiplicity m, contributing a block M_m(ℂ); exponent ×k means k distinct
constituents of that shape. Commutative blocks (m = 1) are written ℂ.

---

## 1. Section A — coset pairs (G, H)

### A₁ family: PSL(2,q) / PGL(2,q)  (the C₃ ⊂ S₃ ⊂ S₄ tower)

| Pair | n | suborbits | rank | Wedderburn | Gelfand? |
|---|---|---|---|---|---|
| PSL(2,5)/S₃ (**Petersen**) | 10 | 1,3,6 | 3 | ℂ³ (χ deg 1,4,5) | ✔ |
| PSL(2,5)/C₃ (double cover of Petersen) | 20 | 1,1,3⁶ | 8 | ℂ⁴ ⊕ **M₂(ℂ)** (1, 3, 3, 4^(2), 5) | ✘ |
| PSL(2,7)/C₃ | 56 | 1,1,3¹⁸ | 20 | ℂ³ ⊕ M₂ ⊕ M₂ ⊕ **M₃** (1, 3, 3, 6^(2), 7^(3), 8^(2)) | ✘ |
| PGL(2,7)/S₄ (**Heawood**) | 14 | 1,3,4,6 | 4 | ℂ⁴ (1, 1, 6, 6) | ✔ |
| PGL(2,7)/D₁₂ (**Coxeter**) | 28 | 1,3,6,6,12 | 5 | ℂ⁵ (1, 6, 6, 7, 8) | ✔ |
| PGL(2,7)/A₄ (other order-12 subgroup) | 28 | 1,1,4,4,6,12 | 6 | ℂ⁶ (1, 1, 6, 6, 7, 7) | ✔ |
| PSL(2,13)/S₃ (**Foster 182**) | 182 | 1,1,3⁴,6²⁸ | 34 | ℂ⁴ ⊕ M₂×3 ⊕ M₃×2 (1, 7, 7, 12^(2)×3, 13^(3), 14, 14^(3)) | ✘ |
| PSL(2,17)/S₄ (**Biggs–Smith**) | 102 | 1,3,6,8,12,24³ | 8 | ℂ⁸ (1, 9, 9, 16, 16, 16, 17, 18) | ✔ |
| PSL(2,19)/S₃ (cr19, 570) | 570 | 1,2,3⁹,6⁹⁰ | 101 | ℂ ⊕ M₂×4 ⊕ M₃×4 ⊕ M₄×3 (1, 9^(2)×2, 18^(2)×2, 18^(4)×2, 19^(3), 20^(3)×3, 20^(4)) | ✘ |
| PSL(2,31)/S₄ (cr31, 620) | 620 | 1,3,6²,8²,12⁵,24²² | 33 | ℂ⁹ ⊕ M₂×6 (1, 30×2, 30^(2)×4, 31, 32×5, 32^(2)×2) | ✘ |

Pattern: the **S₄ stabilizers are Gelfand** for small q (Heawood, Biggs–Smith) but fail at
q = 31; the **C₃ stabilizers are never Gelfand** (the order-2 normalizer quotient
N(C₃)/C₃ acting on the two fixed cosets forces a 1,1 start of the suborbit list, and the
2-dimensional sign-ambiguity propagates into 2×2 blocks). Note PSL(2,19)/S₃ has a
suborbit of length **2** — H meets a conjugate in C₃, so the coset graph there is not simple.

### A₂ family: projective planes (generalized 3-gons)

| Pair | n | suborbits | rank | Wedderburn | Gelfand? |
|---|---|---|---|---|---|
| PSL(3,3)/P₁ (points of PG(2,3)) | 13 | 1,12 | 2 | ℂ² (1, 12) | ✔ (2-transitive) |
| PSL(3,4)/P₁ (points of PG(2,4)) | 21 | 1,20 | 2 | ℂ² (1, 20) | ✔ (2-transitive) |

Doubly transitive ⇒ rank 2 ⇒ H = ℂ², always. The interesting algebra for planes lives on
**flags/Levi graphs** (see Heawood above: rank 4, two linear characters = bipartition).

### B₂ family: generalized quadrangles

| Pair | n | suborbits | rank | Wedderburn | Gelfand? |
|---|---|---|---|---|---|
| PSp(4,3)/P₁ (points of GQ(3,3)) | 40 | 1,12,27 | 3 | ℂ³ (1, 15, 24) | ✔ |

This is the strongly regular graph SRG(40,12,2,4) — rank 3 = distance-regular of diameter 2.

### G₂ family: generalized hexagons — the why-13 home

| Pair | n | suborbits | rank | Wedderburn | Gelfand? |
|---|---|---|---|---|---|
| G₂(2)/P₁ (points of GH(2,2)) | 63 | 1,6,24,32 | 4 | ℂ⁴ (1, 14, 21, 27) | ✔ |
| **G₂(2)/PSL(3,2) (Meinhold-72)** | 72 | 1,1,7⁴,21² | 8 | **ℂ² ⊕ M₂(ℂ) ⊕ ℂ²** (1, 1, 14^(2), 21, 21) | ✘ |

The Meinhold-72 M₂(ℂ) block (χ = 14a with multiplicity 2) is where **B² = 13I** lives —
fully worked out in `hecke/why13.pdf` (13 = 7 + 2·3 via NNᵀ = 7I + 2·A_SRG(36,21,12,12)).

### Future-direction groups (CATALOGUE §9.2)

| Pair | n | suborbits | rank | Wedderburn | Gelfand? |
|---|---|---|---|---|---|
| PSU(3,3)/H₂₈ (Hermitian unital) | 28 | 1,27 | 2 | ℂ² (1, 27) | ✔ (2-transitive) |
| Sz(8)/H₆₅ (Suzuki ovoid) | 65 | 1,64 | 2 | ℂ² (1, 64) | ✔ (2-transitive) |
| J₂/U₃(3) (SRG(100,36,14,12)) | 100 | 1,36,63 | 3 | ℂ³ (1, 36, 63)* | ✔ (rank 3) |

\* J₂ exceeded the Irr size limit in the census; rank 3 from suborbits. The decomposition
π = 1 + 36a + 63a is classical (rank-3 graph; the 36 here is the *same* SRG(36,…) family
echo as in why13 — J₂ ⊃ U₃(3) extends the G₂(2) ⊃ PSU(3,3) story one level up).

---

## 2. Section B — Foster-census graphs (G = Aut Γ, H = vertex stabilizer)

| Graph (n) | |Aut| | suborbits | rank | Wedderburn | Gelfand? |
|---|---|---|---|---|---|
| K₄ (4) | 24 | 1,3 | 2 | ℂ² (1, 3) | ✔ |
| K₃,₃ (6) | 72 | 1,2,3 | 3 | ℂ³ (1, 1, 4) | ✔ |
| Cube (8) | 48 | 1,1,3,3 | 4 | ℂ⁴ (1, 1, 3, 3) | ✔ |
| Petersen (10) | 120 | 1,3,6 | 3 | ℂ³ (1, 4, 5) | ✔ |
| Heawood (14) | 336 | 1,3,4,6 | 4 | ℂ⁴ (1, 1, 6, 6) | ✔ |
| Möbius–Kantor (16) | 96 | 1,1,2,3,3,6 | 6 | ℂ⁶ (1, 1, 3, 3, 4, 4) | ✔ |
| Pappus (18) | 216 | 1,2,3,6,6 | 5 | ℂ⁵ (1, 1, 4, 6, 6) | ✔ |
| Desargues (20) | 240 | 1,1,3,3,6,6 | 6 | ℂ⁶ (1, 1, 4, 4, 5, 5) | ✔ |
| Dodecahedron (20) | 120 | 1,1,3,3,6,6 | 6 | ℂ⁶ (1, 3, 3, 4, 4, 5) | ✔ |
| Nauru (24) | 144 | 1,2,3,3,3,6,6 | 7 | ℂ⁷ (1, 1, 3, 3, 4, 6, 6) | ✔ |
| **F26A (26)** | 78 | 1,1,3⁸ | 10 | **ℂ² ⊕ M₂(ℂ) ⊕ M₂(ℂ)** (1, 1, 6^(2), 6^(2)) | ✘ |
| Tutte–Coxeter (30) | 1440 | 1,3,6,8,12 | 5 | ℂ⁵ (1, 1, 9, 9, 10) | ✔ |
| Tutte 12-cage (126) | **12096** | — | — | **not vertex-transitive**: orbits 63 + 63 | — |
| ↳ each 63-orbit | | 1,6,24,32 | 4 | ℂ⁴ (1, 14, 21, 27) | ✔ |

Sanity checks all passed: girths 3,4,4,5,6,6,6,6,5,6,6,8,12 and the order-only Aut
attributions in CATALOGUE §6 (MK 96, Pappus 216, Nauru 144, F26A 78) confirmed by GRAPE.

---

## 3. Highlights and corrections

### 3.1 The Tutte 12-cage is not vertex-transitive (CATALOGUE §6 corrected)

`AutGroupGraph` gives |Aut| = 12096 = |G₂(2)| with **two orbits of 63** (points and lines
of GH(2,2)). There is no graph automorphism swapping them: a point–line duality of H(q)
exists only for q = 3^h (the exceptional triality-flavoured duality), and q = 2 is not such.
The old "G₂(2)·2 / F126A vertex-transitive" entry was wrong; the Foster census lists F126A
as arc-transitive *on each bipartition class action*, but the full graph is merely
edge- but not vertex-transitive (a "semisymmetric-style" bipartite situation). Each orbit
carries the same Gelfand ℂ⁴ algebra as G₂(2)/P₁ — points and lines are abstractly
isomorphic G₂(2)-sets without being graph-equivalent.

### 3.2 PSL(2,5)/C₃ — a 20-point miniature of Meinhold-72

Exactly the same shape as the why-13 algebra one level down the tower:

| | Meinhold-72 | PSL(2,5)/C₃ |
|---|---|---|
| n | 72 | 20 |
| rank | 8 | 8 |
| suborbits | 1,1,7,7,7,7,21,21 | 1,1,3,3,3,3,3,3 |
| Wedderburn | ℂ²⊕M₂⊕ℂ² | ℂ⁴⊕M₂ |
| 2×2 block from | χ₁₄ (m=2) | χ₄ (m=2) |

**Answered (why-N junior):** explicit matrices in `hecke/HECKE-IRREPS.md`
(computed by `hecke/hecke_irreps.g`, 2026-06-10). The role of 13 is played by 4:
the symmetric generator satisfies (ρ(T₆) − I)² = 4I. There is **no universal formula**
N = k + 2√(k−µ) — the SRG mechanism 13 = 7 + 2·3 is special to the Meinhold geometry
(PSL(2,7)/C₃ even has N = 2 < k = 3).

### 3.3 The non-Gelfand census

Noncommutative Hecke algebras found (each M_m block, m ≥ 2, is a "hidden matrix unit"
worth a why-N question):

| Pair | rank | largest block | source of multiplicity |
|---|---|---|---|
| PSL(2,5)/C₃ | 8 | M₂ | χ₄ twice |
| PSL(2,7)/C₃ | 20 | M₃ | χ₇ three times |
| PSL(2,13)/S₃ | 34 | M₃ ×2 | χ₁₃, χ₁₄ |
| PSL(2,19)/S₃ | 101 | M₄ ×3 | three χ₂₀'s |
| PSL(2,31)/S₄ | 33 | M₂ ×6 | — first S₄ failure |
| G₂(2)/PSL(3,2) | 8 | M₂ | χ₁₄ — **B²=13I** |
| F26A (Aut, 26) | 10 | M₂ ×2 | two χ₆'s of C₁₃⋊C₆ |

Pattern: Gelfand holds when H is "large enough relative to its normalizer structure"
(maximal-ish S₄, point stabilizers, vertex stabilizers of small Foster graphs); it fails
for the C₃ stabilizers (free normalizer action on fixed points), for big coset spaces in
the tower, and for F26A whose automorphism group C₁₃⋊C₆ is sharply flag-transitive with
no room for multiplicity-freeness (78 = 3·26 exactly arc-transitive).

### 3.4 Gelfand pairs = distance-transitive backbone

Every distance-transitive graph in the catalogue gives a Gelfand pair whose rank is
diameter + 1, and the constituent degrees are the multiplicities of the graph eigenvalues:
Petersen ℂ³ (1,4,5), Heawood ℂ⁴, Tutte–Coxeter ℂ⁵, Biggs–Smith ℂ⁸ (diameter 7!),
GH(2,2) points ℂ⁴ (1,14,21,27 ↔ eigenvalues 6, 3²¹, −1²⁷... matched in why13.tex §3),
GQ(3,3) points ℂ³, Foster baseline ≤ 30 all commutative except F26A.

---

## 4. Where to go next

1. ~~**why-N junior**~~ DONE — see §3.2 and `hecke/HECKE-IRREPS.md` (all small irreps
   of the non-Gelfand cases: exact matrices, relations, charpolys).
2. **F26A** (§3.3): two 2×2 blocks over a Frobenius group C₁₃⋊C₆ — the smallest "wild"
   Foster graph. DONE in HECKE-IRREPS.md §3: the blocks are Galois conjugates over ℚ(√13),
   with T₃² = ((5−√13)/2)·I and N₊N₋ = 3 = valency. The 13 here comes from the Gauss
   periods of ℚ(ζ₁₃) — a genuinely *different* 13 from Meinhold's 7+2·3.
3. **J₂/U₃(3)**: rank-3 Gelfand here, but the *flag* version (J₂ on cosets of L₂(7)-type
   subgroups, mirroring Meinhold) is unexplored — does the 36→72 chirality doubling lift?
4. **³D₄(2) hexagon GH(2,8)** (819 points, CATALOGUE §8): the point-graph Hecke algebra is
   ℂ⁴ by Gelfand; the analogue of the Meinhold construction needs a subgroup playing PSL(3,2)'s role.
5. Re-run with `IRR_SIZE_LIMIT` raised to get J₂ multiplicities from GAP's character
   table library (`CharacterTable("J2")` + `PermChars`) rather than by citation.

---

*Generated from `hecke/hecke_census.g`; raw output `hecke/hecke_census.out`.
See `hecke/why13.pdf` for the full G₂(2) story and `hecke/MEINHOLD-HECKE.md`
for the rank-8 multiplication table.*
