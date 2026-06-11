# Hecke Algebras of the Cubic Arc-Transitive Catalogue

For every coset pair (G, H) producing a cubic arc-transitive graph, compute the
Hecke algebra H(G, H) = End_G(C[G/H]) — its Wedderburn decomposition, spherical
functions, and explicit matrix representations of every irreducible component.

## The story in three sentences

The Hecke algebra of a transitive permutation representation decomposes as
⊕ M_{m_χ}(C), with one block per constituent of ind_H^G 1. When every multiplicity
m_χ ≤ 1 (a Gelfand pair), the algebra is commutative and the theory reduces to
eigenvalue tables. When it isn't, you get 2×2 and 3×3 matrix blocks — and that's
where the interesting numbers live.

## Key results

- **B² = 13I** in the Meinhold-72 Hecke algebra (G₂(2)/PSL(3,2)), where
  13 = 7 + 2×3 via the strongly regular graph SRG(36, 21, 12, 12). Full analysis
  in [why13.pdf](why13.pdf).

- **The N-invariant**: traceless generators of 2D blocks satisfy B² = NI for
  varying N across different pairs — but no universal formula exists. See the
  table in [HECKE-IRREPS.md](HECKE-IRREPS.md) §5.

- **Spectral partition**: the Wedderburn decomposition sorts every graph eigenvalue
  into a specific Hecke block. The dodecahedron's ±√5 lives entirely in 1D irreps;
  the 2D block contributes only integers.

- **Finite-order Hecke elements**: PSL(2,7)/C₃ has a generator of order 3 in
  its rational 2D block — unexpected for Hecke algebras.

## Papers (LaTeX)

| File | Title |
|------|-------|
| [hecke.tex](hecke.tex) | The Hecke Algebra of a Pair (G, H) — foundations |
| [why13.tex](why13.tex) | Why 13? The Meinhold-72 Hecke algebra and G₂(2) |
| [hecke-irreps.tex](hecke-irreps.tex) | Explicit irreducible representations of the catalogue |
| [hecke-catalog.tex](hecke-catalog.tex) | Full Hecke census: Wedderburn types for every graph |

## GAP scripts

### Census and catalogue

| Script | What it computes |
|--------|-----------------|
| `hecke_census.g` | Wedderburn decomposition (ranks, multiplicities) for every pair |
| `hecke_catalog_full.g` | Full Hecke tables with spherical functions → LaTeX fragments in `catalog-data/` |
| `hecke_catalogue_gaps.g` | Identifies which catalogue entries are non-Gelfand |
| `hecke_spherical.g` | Spherical function tables for all Gelfand pairs |

### Irreducible representations

| Script | What it computes |
|--------|-----------------|
| `hecke_irreps.g` | Exact 2D/3D matrix irreps via isotypic projection (PSL(2,5)/C₃, PSL(2,7)/C₃, F26A, G₂(2)/PSL(3,2)) |
| `hecke_irreps2.g` | Extended irrep computations |
| `hecke_irreps3.g` | PSL(2,31)/S₄ and larger cases |

### Meinhold-72 deep dives

| Script | What it computes |
|--------|-----------------|
| `meinhold_hecke.g` | Full Hecke algebra of G₂(2)/PSL(3,2): structure constants, eigenvalues |
| `meinhold_hecke_tables.g` | Multiplication tables for various stabilizer pairs in G₂(2) |
| `meinhold72_2d_irrep.g` | The 2D irrep: A² = I, B² = 13I, AB + BA = −4I |
| `meinhold72_irreps.g` | All irreps of the rank-8 algebra |
| `meinhold72_center.g` | Center of the Hecke algebra |
| `meinhold72_codes.g` | Linear codes from Meinhold-72 orbital matrices |
| `meinhold72_codes2.g` | Extended code analysis |
| `meinhold72_generator.g` | Generator-relation presentation |

### Other

| Script | What it computes |
|--------|-----------------|
| `sp_cubic_census.g` | Cubic graph search in PSp(4, q) — found Schläfli graph, rank barrier |
| `why13_orbitals.g` | Orbital structure of G₂(2) on 72 vertices |
| `why13_classes.g` | Conjugacy class analysis for why13 |
| `why13_fano.g` | Fano plane structure inside the G₂(2) geometry |
| `why13_subhexagon.g` | Subhexagon analysis |
| `klein_adj.g` | Klein graph adjacency computations |
| `draw_schlafli.py` | Python script to draw the Schläfli graph |

## Documentation

| File | Contents |
|------|----------|
| [HECKE-IRREPS.md](HECKE-IRREPS.md) | Full explicit irrep data with relations and the N-invariant table |
| [MEINHOLD-HECKE.md](MEINHOLD-HECKE.md) | Detailed Meinhold-72 Hecke algebra analysis |
| [RESEARCH_FINDINGS.md](RESEARCH_FINDINGS.md) | PSp(4,3) → Schläfli graph, Wan's Anzahl connection |
| [no-more-hecke.md](no-more-hecke.md) | Blog post: "No more Hecke algebras, I've had enough..." |

## Dependencies

- [GAP](https://www.gap-system.org/) 4.15+ with the [GRAPE](https://gap-packages.github.io/grape/) package
- Python 3 with matplotlib (for `draw_schlafli.py`)
- LaTeX with amsmath, booktabs, hyperref (for the papers)
