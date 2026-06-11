# Explicit irreps of the non-Gelfand Hecke algebras: generators, relations, matrix entries

**Computed by:** `hecke_irreps.g` → `hecke_irreps.out` (GAP 4.15.1 + GRAPE, 2026-06-10).
**Companion to:** `~/git/graphviz/docs/HECKE_RESEARCH.md` (the census), `MEINHOLD-HECKE.md`,
`why13.pdf`.

## 0. Method and conventions

For each pair (G,H), the orbital basis T₁,…,T_r of H(G,H) acts on ℂ[G/H] by the
n×n orbital matrices A_i. For a constituent χ of π = ind_H^G 1 with multiplicity m,
the χ-isotypic component W ≅ (m-dim Hecke irrep) ⊗ ℂ^{χ(1)}, computed exactly via the
projector P_χ = (χ(1)/|G|) Σ_g χ(g⁻¹)π(g). Any eigenvector of an orbital matrix inside W,
for an eigenvalue that is simple in the irrep, is a pure tensor; its closure under
{A_i} is a single m-dimensional irreducible H-module. Restricting gives exact m×m
matrices over cyclotomics.

Conventions: **T₁ = identity orbital**; suborbits sorted by length (so the numbering
differs from MEINHOLD-HECKE.md / why13.tex, see §4). Matrices act on **row** vectors.
Basis = eigenbasis of the generator named in each section; entries are exact.
For m = 1 the irrep is the spherical function χ ↦ (ω(T₁),…,ω(T_r)).
Verified in every case: dim of the algebra generated = m², and ρ(Σ T_i) = 0
(image of the all-ones matrix on a nontrivial block).

Notation for quadratic irrationalities:

| symbol | value | satisfies |
|---|---|---|
| α | (−1+√−7)/2 = ζ₇+ζ₇²+ζ₇⁴ | α² + α + 2 = 0, αᾱ = 2 |
| η₀ | (−1+√13)/2 = Σ ζ₁₃^{QR} | η² + η − 3 = 0, η₀η₁ = −3, η₀−η₁ = √13 |
| η₁ | (−1−√13)/2 = Σ ζ₁₃^{NQR} | (Galois conjugate of η₀) |
| φ-sum | √5 = (ζ₅+ζ₅⁴)−(ζ₅²+ζ₅³) | Gauss sum mod 5 |

---

## 1. PSL(2,5)/C₃ (n = 20) — H ≅ ℂ⁴ ⊕ M₂(ℂ)

Suborbits (1,1,3,3,3,3,3,3). The drawn coset graph here is the **dodecahedron**
(as a graph): see the spectral reassembly below.

**1-dim irreps** (values on T₁..T₈):

| on χ | values |
|---|---|
| χ₁ (deg 1) | 1, 1, 3, 3, 3, 3, 3, 3 |
| χ₂ (deg 3) | 1, −1, **−√5, +√5**, 1, −1, 1, −1 |
| χ₃ (deg 3) | 1, −1, **+√5, −√5**, 1, −1, 1, −1 |
| χ₅ (deg 5) | 1, 1, 1, 1, −1, −1, −1, −1 |

**2-dim irrep** (on χ₄, deg 4; basis = eigenbasis of ρ(T₂)):

```
T₁ = I            T₂ = [ 1  0 ]        T₃ = T₄ = [ -2  0 ]
                       [ 0 -1 ]                  [  0  0 ]

T₅ = [  1/2    1 ]   T₆ = [  1/2   -1 ]   T₇ = [ 1/2   -1 ]   T₈ = [ 1/2    1 ]
     [ -15/4 -3/2]        [ -15/4  3/2]        [ 15/4 -3/2]        [ 15/4  3/2]
```

Characteristic polynomials: T₂: x²−1; T₃,T₄: x(x+2); T₅,T₇: **x²+x+3** (complex,
det = 3 = valency); T₆,T₈: (x−3)(x+1).

**Relations.** A := ρ(T₂), C := ρ(T₅): A² = I, C² + C + 3I = 0, and the deck
involution permutes the orbitals inside the block: AC = ρ(T₈), CA = ρ(T₆),
ACA = ρ(T₇). The traceless shift (ρ(T₆) − I)² = **4I** is the closest analogue of
B² = 13I here; note ρ(T₃) = ρ(T₄) (the block does not separate these two suborbits;
they are separated by χ₂, χ₃ via ±√5).

**Spectral reassembly (why √5):** the symmetric orbital T₃ has eigenvalue χ-data
{3 | ±√5 | 0, −2 | 1} with multiplicities {1 | 3,3 | 4,4 | 5} — exactly the
dodecahedron spectrum 3, √5³, 1⁵, 0⁴, (−2)⁴, (−√5)³. The dodecahedron's famous ±√5
lives in the **1-dimensional** irreps (the degree-3 characters of A₅); the 2×2 block
contributes only the integer eigenvalues.

---

## 2. PSL(2,7)/C₃ (n = 56) — H ≅ ℂ³ ⊕ M₂ ⊕ M₂ ⊕ M₃

Suborbits (1,1,3¹⁸); rank 20. Three small irreps (all matrices below exact).

### 2a. 2-dim irrep on χ₄ (deg 6) — the ℚ(√−7) block

ρ(T₂) = **I**: the free involution acts trivially here. Distinct images
(α = (−1+√−7)/2):

| matrix | which T_i | charpoly |
|---|---|---|
| diag(α, ᾱ) | T₃, T₅, T₁₂, T₁₇ | x²+x+2 |
| diag(ᾱ, α) | T₄, T₆, T₇, T₁₄ | x²+x+2 |
| [[1,2],[2,1]] | T₈, T₁₈ | (x−3)(x+1) |
| [[1, α−1],[ᾱ−1, 1]] | T₉, T₁₅ | (x−3)(x+1) |
| [[0, ᾱ],[α, 0]] | T₁₀, T₁₁, T₁₆, T₂₀ | **x² − 2** |
| [[−1, −ᾱ],[−α, −1]] | T₁₃, T₁₉ | x²+2x−1 |

**Relations.** D := diag(α, ᾱ), S := [[0,ᾱ],[α,0]]: D² + D + 2I = 0,
**S² = 2I** (2 = αᾱ = Nm(α), the prime 2 split in ℚ(√−7) — the binary Hamming/Paley
phenomenon at q=7), and SDS⁻¹ = D̄ (S swaps the eigenbasis, conjugating D to its
Galois conjugate). ρ(T₁₃) = −I − S. The six matrices in the table span M₂(ℂ), so
{D, S} generate the block. Eigenvalues −1±√2 appear (T₁₃): **√2 = the junior √13**.

### 2b. 3-dim irrep on χ₅ (deg 7) — all-integer matrices

Basis = eigenbasis of ρ(T₂) = diag(1,−1,−1). Distinct images:

```
T₃,(T₁₂ᵀ-pattern): [ 1 1 0 ]   T₄: [ 1 0 1 ]   T₅: [  1  1  0 ]   T₆: [  1  0  1 ]
                   [ 2 0 2 ]       [ 4 1 0 ]       [ -2  0 -2 ]       [ -4 -1  0 ]
                   [ 4 0 1 ]       [ 2 2 0 ]       [ -4  0 -1 ]       [ -2 -2  0 ]

T₇: [  1  0 -1 ]   T₈: [ -1  0  0 ]   T₉: [ -1  0  0 ]   T₁₀: [ -1 -1  1 ]
    [  4 -1  0 ]       [  0  1  0 ]       [  0  3 -2 ]        [  2  0  0 ]
    [  2 -2  0 ]       [  0  2 -3 ]       [  0  0 -1 ]        [ -2  0  0 ]

T₁₁: [ -1  1 -1 ]   T₁₃: [ -1  0  0 ]   (T₁₂,T₁₄..T₂₀ are the remaining sign/
     [  2  0  0 ]        [  0 -1  2 ]    transpose variants — full list in
     [ -2  0  0 ]        [  0 -2  1 ]    hecke_irreps.out lines 150–229)
```

Characteristic polynomials: x³−2x²−x−6 (T₃,T₄,T₁₂,T₁₄); x³+x−6 (T₅,T₆,T₇,T₁₇,
irreducible over ℚ — genuinely cubic eigenvalues); (x+3)(x−1)(x+1) (T₈,T₁₅);
(x−3)(x+1)² (T₉,T₁₈); x(x²+x+4) (T₁₀,T₁₆); **x(x²+x−4)** (T₁₁,T₂₀ — eigenvalues
0, (−1±**√17**)/2); (x+1)(x²+3) (T₁₃,T₁₉ — eigenvalues −1, ±i√3).

### 2c. 2-dim irrep on χ₆ (deg 8) — the rational block

Basis = eigenbasis of ρ(T₂) = diag(1,−1). Distinct images:

| matrix | which T_i | charpoly |
|---|---|---|
| [[−1/2, 1],[−7/4, −1/2]] | T₃, T₁₄ | x²+x+2 |
| [[−1/2, −1],[7/4, −1/2]] | T₄, T₁₂ | x²+x+2 |
| [[−1/2, 1],[7/4, 1/2]] | T₅, T₇ | **x² − 2** |
| [[−1/2, −1],[−7/4, 1/2]] | T₆, T₁₇ | **x² − 2** |
| diag(−1, 1) | T₈, T₁₅ | x²−1 |
| −I | T₉, T₁₈ | (x+1)² |
| [[1/2, 1],[−7/4, −3/2]] | T₁₀ (& T₁₆ variant) | **x²+x+1** (ρ³ = I!) |
| [[1/2, −1],[−7/4, 3/2]] | T₁₁ (& T₂₀ variant) | x²−2x−1 (eigenvalues 1±√2) |
| [[2, 0],[0, 0]] | T₁₃, T₁₉ | x(x−2) |

**Relations.** B := ρ(T₅): **B² = 2I**. The same x²+x+2 (norm-2) element as in χ₄,
now over ℚ. Remarkable: ρ(T₁₀) has **finite order 3** — a Hecke generator acting as a
rotation. And 1±√2 again: the ℤ[√2] thread runs through both 2-dim blocks.

---

## 3. F26A (n = 26, Aut = C₁₃⋊C₆) — H ≅ ℂ² ⊕ M₂ ⊕ M₂

Suborbits (1,1,3⁸). 1-dim irreps: χ₁: (1,1,3,…,3); χ₂: (1,−1,−3,3,−3,3,−3,3,−3,3).
The two 2-dim blocks (on the two degree-6 characters) are **Galois conjugate over
ℚ(√13)**: swap η₀ ↔ η₁ to pass between them.

**2-dim irrep on χ₇** (basis = eigenbasis of ρ(T₂) = diag(1,−1); γ := 2η₀ + (5/4)η₁):

```
T₃ = [  η₀/2   1  ]    T₄ = [  η₁/2    −(η₀+2η₁) ]    T₅ = [  η₀/2  −1  ]
     [  −γ   −η₀/2 ]        [ −√13/4     η₁/2    ]         [   γ   −η₀/2 ]

T₆ = [  η₀/2  −1  ]    T₇ = [  η₁/2     η₀+2η₁ ]      T₈ = [  η₁/2    η₀+2η₁ ]
     [  −γ    η₀/2 ]        [ −√13/4   −η₁/2   ]           [  √13/4   η₁/2   ]

T₉ = [  η₁/2   −(η₀+2η₁) ]   T₁₀ = [ η₀/2   1  ]
     [  √13/4   −η₁/2    ]         [  γ    η₀/2 ]
```

Characteristic polynomials (recall 2−η₀ = (5−√13)/2, 2−η₁ = (5+√13)/2):

| T_i | charpoly | relation |
|---|---|---|
| T₃, T₅ | x² + (η₀−2) | **T² = ((5−√13)/2)·I** |
| T₇, T₉ | x² + (η₁−2) | **T² = ((5+√13)/2)·I** |
| T₄, T₈ | x² − η₁x + (2−η₁) | T² = η₁T − (2−η₁)I |
| T₆, T₁₀ | x² − η₀x + (2−η₀) | T² = η₀T − (2−η₀)I |

**The F26A "why-N":** the traceless generators square to scalars
N± = (5∓√13)/2, with **N₊N₋ = 3** (the valency) and N₊+N₋ = 5. So √13 — the
*same* 13 as in Meinhold-72, but for the trivial reason 26 = 2·13 here — enters
through the quadratic Gauss periods of ℚ(ζ₁₃), and the two Galois-conjugate blocks
share the rational data (3, 5).

---

## 4. G₂(2)/PSL(3,2) — Meinhold-72 cross-check ✓

(Numbering here is sorted-by-valency; the dictionary to MEINHOLD-HECKE.md/why13.tex is
T₂ = the free involution "A", T₄/T₆ = the symmetric valency-7 pair "B/B′" — old T₃/T₅;
T₃/T₅ = the directed pair — old T₄/T₆.)

**1-dim irreps:** χ₁: (1,1,7,7,7,7,21,21); χ₂: (1,−1,7,−7,7,−7,−21,21);
χ₂₁: (1,−1,1,−1,1,−1,3,−3); χ₂₁′: (1,1,1,1,1,1,−3,−3).

**2-dim irrep on χ₁₄** (basis = eigenbasis of ρ(T₂)):

```
T₁ = I             T₂ = [ 1  0 ]      T₃ = [ -2  1 ]     T₄ = [ -2 -1 ]
                        [ 0 -1 ]           [ -9 -2 ]          [ -9  2 ]

T₅ = [ -2 -1 ]     T₆ = [ -2  1 ]     T₇ = [ 3  0 ]      T₈ = 3·I
     [  9 -2 ]          [  9  2 ]          [ 0 -3 ]
```

Charpolys: T₃,T₅: **x²+4x+13** (the −2±3i pair, |−2+3i|² = 13); T₄,T₆: **x²−13**;
T₇: x²−9; T₈ = 3I (scalar, as why13 §7 predicts).

**Relations reproduced exactly:** with A := ρ(T₂), B := ρ(T₄):
A² = I, **B² = 13I**, **AB + BA = −4I**, and BA = ρ(T₃), AB = ρ(T₅) (the directed
pair is literally B conjugated/translated by the chirality involution). This is an
independent re-derivation of why13.tex Theorem 8.1's algebra by a different method
(isotypic projection rather than orbital path-counting).

---

## 5. The emerging pattern: which N appears in B² = N·I?

| algebra | traceless generator | B² = N·I | where N lives |
|---|---|---|---|
| G₂(2)/L₂(7), χ₁₄ | valency-7 symmetric | N = **13** = 7 + 2·3 | ℤ; |−2±3i|² = 13 in ℤ[i] |
| PSL(2,7)/C₃, χ₄ & χ₆ | valency-3 (two blocks) | N = **2** | ℤ; 2 = αᾱ = Nm((−1+√−7)/2) |
| F26A, χ₇/χ₈ | valency-3 | N = **(5∓√13)/2** | ℤ[η], Nm(N) = 3 = valency |
| PSL(2,5)/C₃, χ₄ | (shifted) ρ(T₆)−I | N = **4** | ℤ; complex pair has det 3 = valency |

In every case the directed/complex orbitals have determinant or norm tied to the
valency, while the symmetric generator squares to an N that measures the
"interference" between the two chirality classes — for Meinhold-72 the precise
statement is NNᵀ = 7I + 2·A_SRG(36,21,12,12) ⇒ 13 = 7 + 2·3 (why13.tex §8).
Junior open question upgraded to data: there is **no universal formula N = k + 2√(k−µ)
yet** — PSL(2,7)/C₃'s N = 2 < k = 3 shows the SRG mechanism is special to the
Meinhold geometry. What is the general synthetic meaning of N? (→ why13.tex §9.)

Side findings worth keeping: a Hecke generator of **finite order 3** (PSL(2,7)/C₃, χ₆);
genuinely **cubic** eigenvalues (x³+x−6) in the 3-dim irrep; **√17** in PSL(2,7)/C₃'s
3-dim block; the dodecahedron's **√5** living entirely in the 1-dim irreps.

---

## 6. The PSL(2,q)/C₃ family: 2D blocks exist only for q = 5 and 7

Computed by `hecke_psl2_survey.g` (2026-06-11). The Wedderburn decomposition of
H(PSL(2,q), C₃) was computed for all primes 5 ≤ q ≤ 47. The minimum multiplicity
m > 1 (= dimension of the smallest non-commutative block) grows with q:

| q | n | rank | Wedderburn blocks | Min mult > 1 | Steinberg mult |
|---|---|------|-------------------|-------------|----------------|
| 5 | 20 | 8 | C⁴ ⊕ M₂ | **2** | 1 |
| 7 | 56 | 20 | C³ ⊕ M₂² ⊕ M₃ | **2** | 3 |
| 11 | 220 | 76 | C³ ⊕ M₃ ⊕ M₄⁴ | 3 | 3 |
| 13 | 364 | 124 | C ⊕ M₃² ⊕ M₄³ ⊕ M₅ ⊕ M₄² | 3 | 5 |
| 17 | 816 | 276 | C ⊕ M₃² ⊕ M₄ ⊕ M₅ ⊕ M₆⁶ | 3 | 5 |
| 19 | 1140 | 384 | C ⊕ M₃² ⊕ M₆⁴ ⊕ M₇ ⊕ M₈ ⊕ M₆³ | 3 | 7 |
| 23 | 2024 | 680 | C ⊕ M₃² ⊕ M₆ ⊕ M₇ ⊕ M₈⁹ | 3 | 7 |
| 29 | 4060 | 1360 | C ⊕ M₅² ⊕ M₈² ⊕ M₉ ⊕ M₁₀¹¹ | 5 | 9 |
| 31 | 4960 | 1660 | C ⊕ M₅² ⊕ M₁₀⁷ ⊕ M₁₁ ⊕ M₁₂² ⊕ M₁₀⁴ | 5 | 11 |

**The Steinberg multiplicity has an exact formula.** An element g of order 3 in
PSL(2,q) is semisimple; its eigenvalues live in F_q when 3|(q−1) and in F_{q²}
when 3|(q+1), giving Steinberg character value +1 or −1 respectively:

- q ≡ 1 (mod 3): m_St = (q + 2)/3
- q ≡ 2 (mod 3): m_St = (q − 2)/3

**Why mult 2 is impossible for q ≥ 11.** The multiplicity of χ in ind_{C₃}^G 1
is m_χ = dim(χ^{C₃}) = (χ(1) + χ(g) + χ(g²))/3. For m = 2 we need
χ(1) + χ(g) + χ(g²) = 6. As q grows, all character degrees grow (the smallest
non-trivial degrees are (q−1)/2 for q ≥ 11), and the character values on
order-3 elements are bounded, so m ≥ 3 is forced. **The B² = NI phenomenon
in the PSL(2,q)/C₃ family is confined to exactly two primes: q = 5 and q = 7.**

### PSL(2,11)/C₃ detail (computed by `hecke_psl2_11.g`)

Wedderburn: C(1) ⊕ C(5) ⊕ C(5) ⊕ M₄(10)² ⊕ M₃(11) ⊕ M₄(12)².

The 1-dim irreps on the degree-5 pair are Galois conjugates over **ℚ(√−11)**,
with spherical values involving Gauss periods η₀ = (−1+√−11)/2. The 3-dim
block (Steinberg, degree 11) has eigenvalues involving **√5** and the **golden
ratio** (1±√5)/2 — the same irrationality as the dodecahedron, now appearing
in the non-commutative block rather than the spherical functions.

Four orbitals act as −I on the 3-dim block (T₃₆ = T₆₀ = T₆₅ = T₆₆ = −I),
and elements of **finite order 6** appear (charpoly x² − x + 1, primitive
6th roots of unity).

---

## 7. Torsion in the non-commutative blocks

Computed by `hecke_torsion.g` (2026-06-11). For each non-commutative Hecke block,
what is the maximum order of ρ(T_i)?

| q | Block | χ (deg) | Max torsion | Cyclotomic |
|---|-------|---------|-------------|------------|
| 7 | M₂ | χ₄ (6) | 1 (no torsion) | — |
| 7 | M₃ | χ₅ (7) | 2 | Φ₁Φ₂ |
| 7 | M₂ | χ₆ (8) | **3** | **Φ₃** (x²+x+1) |
| 11 | M₃ | χ₆ (11) | **6** | x³+1 = Φ₂·**Φ₆** |
| 13 | M₃ | χ₂ (7) | 2 | Φ₁Φ₂ |
| 13 | M₃ | χ₃ (7) | 2 | Φ₁Φ₂ |

**Observations:**

1. **Torsion concentrates in specific blocks.** At q=7, the order-3 elements live
   in the rational M₂ block (χ₆, deg q+1), not the ℚ(√−7) block or the M₃.
   At q=11, order-6 elements live in the Steinberg M₃ block. At q=13, the
   M₃ blocks (on the deg-(q+1)/2 characters) have only involutions.

2. **The Steinberg block is special.** At q=11, eight orbitals act with order 6
   on the Steinberg block (charpoly x³+1 = (x+1)(x²−x+1)), and four act as −I.
   The 3×3 matrix has eigenvalues {−1, e^{πi/3}, e^{−πi/3}} — a 60° rotation
   in a 2D subspace. This connects to the Eisenstein integers ℤ[ω₃].

3. **Torsion order is NOT monotone in q.** q=7 gives order 3, q=11 gives order 6,
   but q=13 drops back to order 2 (in the blocks computed so far). The torsion
   depends on which block you're in, not just on q.
