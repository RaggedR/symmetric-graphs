# Spectral Rigidity of Tessellation Graphs

## Result

**Theorem (computational).** For both surfaces tested:

| Surface | Tessellation graph | |Aut(G)| | |Aut_full(X)| | Match? |
|---------|-------------------|---------|--------------|--------|
| Klein quartic (g=3, Hurwitz) | F056A (56v, cubic) | 336 | 336 | ✓ |
| Bolza (g=2, non-Hurwitz) | F016 Möbius-Kantor (16v, cubic) | 96 | 96 | ✓ |

In both cases: `|Aut(tessellation graph)| = |Aut_full(surface)|` including
anti-holomorphic/orientation-reversing isometries.

## The Spectral Mechanism

### Klein quartic (F056A)

The permutation representation of PSL(2,7) on 56 vertices decomposes as:
```
ℝ⁵⁶ = 1 ⊕ V₃ ⊕ V₃' ⊕ V₆ ⊕ V₆ ⊕ V₇ ⊕ V₇ ⊕ V₇ ⊕ V₈ ⊕ V₈
```

The 9 eigenvalues of the adjacency matrix and their irrep content:

| λ | dim | PSL(2,7)-content | PGL(2,7)-content |
|---|-----|------------------|------------------|
| 3 | 1 | trivial | trivial |
| 1+√2 | 8 | V₈ | V₈ (irred) |
| ≈1.56 | 7 | V₇ | V₇ (irred) |
| √2 | 6 | V₆ | V₆ (irred) |
| 0 | 7 | V₇ | V₇ (irred) |
| 1-√2 | 8 | V₈ | V₈ (irred) |
| -√2 | 6 | V₆ | V₆ (irred) |
| **-2** | **6** | **V₃ ⊕ V₃' (reducible!)** | **irred (fused)** |
| ≈-2.56 | 7 | V₇ | V₇ (irred) |

**The critical eigenspace**: At λ = -2, the 6-dimensional eigenspace
decomposes as V₃ ⊕ V₃' under PSL(2,7) — two conjugate 3-dim irreps.
They share an eigenvalue because the adjacency matrix is rational, and
χ₂, χ₃ are complex conjugates (evaluating to the same real value on the
connection double coset).

The outer automorphism (complex conjugation on the Klein quartic, which
extends PSL(2,7) to PGL(2,7) = Aut(F056A)) **swaps V₃ ↔ V₃'**, making
this eigenspace irreducible under the full automorphism group.

### Why this prevents hidden symmetries

Under PGL(2,7), every eigenspace is irreducible. By Schur's lemma,
any matrix commuting with all of PGL(2,7) must act as a scalar on each
eigenspace. But graph automorphisms commute with A and hence preserve
eigenspaces. If they also commute with PGL(2,7)'s action on each
eigenspace, they must be scalar — i.e., they're already in PGL(2,7).

Therefore: |Aut(G)| = |PGL(2,7)| = 336. No hidden symmetries.

### Bolza surface (F016)

```
Aut(F016) = GL(2,3) : C₂     (order 96)
Aut_full(Bolza) = GL(2,3) : C₂  (order 96)
```

The extra C₂ beyond GL(2,3) is the anti-conformal involution
(x,y) ↦ (x̄, ȳ) on the curve y² = x⁵ - x. Same mechanism:
some eigenspace of GL(2,3) is reducible but becomes irreducible
under the full group of order 96.

## Revised Conjecture

**Original (PRELIMINARY.md):** Hurwitz surfaces satisfy |Aut(G)| = |Aut_full(X)|;
non-Hurwitz surfaces can have |Aut(G)| > |Aut_full(X)|.

**Revised:** For ANY surface X = ℍ/Γ with triangle group tessellation,
if the tessellation graph G is a **simple graph** (not a multigraph):

> |Aut(G)| = |Aut_full(X)|

The K₃,₃ "failure" in PRELIMINARY.md was not about Hurwitz vs non-Hurwitz —
it was about the {3,8} dual tessellation degenerating to a multigraph
(|GL(2,3)|/8 = 6 vertices, degree 8 > 5, forcing multi-edges). When we
use the non-degenerate tessellation (the Möbius-Kantor {8,3}), equality holds.

**Condition for simplicity:** The tessellation graph for the {q,r} tiling of
X with Aut(X) = G is simple iff |G|/q > q (enough vertices for the valency).
For Hurwitz surfaces with G = PSL(2,p), |G| = p(p²-1)/2, and the (2,3,7)
tiling gives valency 7 with |G|/7 = p(p²-1)/14 vertices — always simple for p ≥ 7.

## The L-function Connection

The eigenspace V₃ ⊕ V₃' at λ = -2 on the Klein graph is directly related
to H¹(Klein quartic) = V₃ ⊕ V₃' (computed in `klein-quartic/`). This is
the cohomology of the surface — the space of holomorphic differentials
(dim g = 3) and its conjugate.

The Mellin transform of these holomorphic differentials gives the L-function
of the Klein quartic (as a modular curve Γ(7)\ℍ). The CM type
{ζ₇, ζ₇², ζ₇⁴} (quadratic residues mod 7) determines which eigenvalues
appear, and the fact that these are conjugate under Gal(ℚ(ζ₇)/ℚ) is
precisely why V₃ and V₃' share an adjacency eigenvalue.

**The chain:**
```
L-function (Mellin of automorphic forms)
  ↓ determines
Galois action on H¹ = V₃ ⊕ V₃' (conjugate irreps)
  ↓ forces
Shared eigenvalue on the graph (rationality of A)
  ↓ creates
Reducible eigenspace needing the outer automorphism to fuse
  ↓ produces
|Aut(G)| = 2 × |PSL(2,7)| = |Aut_full(X)|
```

The L-function doesn't just "constrain" the graph symmetries — it
**determines the precise mechanism** by which the anti-holomorphic
involution manifests in the graph's spectrum.

## Files

- `spectral_rigidity.py` — v1: eigenvalue dimensions vs irrep dimensions
- `spectral_rigidity_v2.py` — v2: distance-regularity and adjacency algebra
- `spectral_rigidity_v3.py` — v3: direct |Aut(G)| computation
- `klein-quartic/PRELIMINARY.md` — original conjecture (now superseded)
