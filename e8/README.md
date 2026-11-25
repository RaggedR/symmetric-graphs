# The exceptional ladder: E₆ ⊂ E₇ ⊂ E₈ as mandalas and Hecke algebras

Symmetry-aware drawings of the three largest exceptional Weyl-group coset
spaces, each rendered on its Coxeter plane, together with the orbital graphs
(= the basis of its Hecke algebra / association scheme). This extends the
graphviz programme — *cubic arc-transitive graphs as drawable shadows of
groups* — up the exceptional series, and closes it back onto the **Schläfli
graph** already studied in `../hecke/`.

Everything here is computed from scratch and verified before drawing
(`e8_mandala.py`, `e7_mandala.py`, `e6_one.py`, …).

## The three rungs

Each rung is a coset space G/H for a chain of Weyl groups
W(E₈) ⊃ W(E₇) ⊃ W(E₆) ⊃ W(D₅), realized geometrically as a *minuscule
descent*: the vertices of one rung are the vectors of the rung above lying at
60° (inner product +1) to a fixed vector.

| rung | G / H | \|G/H\| | Gosset polytope | Coxeter h | Coxeter-plane rings | scheme rank |
|------|-------|--------:|-----------------|----------:|---------------------|------------:|
| **E₈** | W(E₈)/W(E₇) | 240 | 4₂₁ | 30 | 8 × 30 | 5 |
| **E₇** | W(E₇)/W(E₆) | 56 | 3₂₁ | 18 | 2 + 3 × 18 | 4 |
| **E₆** | W(E₆)/W(D₅) | 27 | 2₂₁ | 12 | 3 + 2 × 12 | 3 |

|W(E₈)|=696,729,600, |W(E₇)|=2,903,040, |W(E₆)|=51,840, |W(D₅)|=1,920; the
indices 696729600/2903040 = **240**, 2903040/51840 = **56**, 51840/1920 =
**27** all check. The scheme rank drops 5 → 4 → 3 as we descend.

## The orbital graphs (the Hecke algebra basis)

For each rung, G acts on G/H and the orbitals (G-orbits on pairs = double
cosets H\G/H) are the inner-product classes. Each orbital is a symmetric graph
on the vertices; together they are a basis of the commutative Hecke algebra
H(G,H) = End_G ℂ[G/H] (a Bose–Mesner algebra, since all orbitals are
self-paired ⇒ Gelfand pair).

### E₈ — 240 roots, rank 5, suborbits {1, 1, 56, 56, 126}

| orbital (inner product) | graph | valency | file |
|---|---|---:|---|
| +2 | identity (the 240 points) | — | `e8_id.pdf` |
| +1 | **the mandala** (roots at 60°, 4₂₁ edge graph) | 56 | `e8_mandala.pdf` |
| 0 | the bloom (orthogonal roots) | 126 | `e8_ip0.pdf` |
| −1 | the weave (roots at 120°) | 56 | `e8_ipm1.pdf` |
| −2 | antipodal matching, 120·K₂ | 1 | `e8_ipm2.pdf` |

Spectrum of the +1 graph: 56¹, 28⁸, 8³⁵, (−2)¹¹², (−4)⁸⁴ (five eigenvalues ⇒
it generates the whole 5-dim algebra). Sample structure constant, verified:
**A₊₁² = 56·I + 27·A₊₁ + 12·A₀ + 1·A₋₁.**

Folding by the root antipode r ↦ −r (the 120 antipodal pairs) gives
**SRG(120, 56, 28, 24)** — the **E₈ root graph**, Aut = **O⁺(8,2):2** =
W(E₈)/{±1}. This is a celebrated strongly-regular graph (a quantum- but not
classically-isomorphic example; *Algebraic Combinatorics* alco.335 /
arXiv:2209.14906). Its complement is SRG(120, 63, 30, 36), and the 126-regular
"bloom" orbital is exactly 2× that complement.

### E₇ — 56 minuscule weights, rank 4, suborbits {1, 1, 27, 27}

| orbital | graph | valency | file |
|---|---|---:|---|
| +1 | **the Gosset graph** | 27 | `e7_mandala.pdf` |
| 0 | the partner 27-regular graph | 27 | `e7_partner.pdf` |
| −1 | **28·K₂ — the 28 bitangent pairs** | 1 | `e7_bitangents.pdf` |

Gosset-graph spectrum: 27¹, 9⁷, (−1)²⁷, (−3)²¹. The valency-1 orbital pairs
each weight with its unique opposite — the 28 pairs are the classical 28
bitangents / double-six structure.

### E₆ — 27 lines, rank 3 ⇒ STRONGLY REGULAR

| orbital | graph | valency | file |
|---|---|---:|---|
| 0 | **the Schläfli graph, SRG(27,10,1,5)** | 10 | `e6_schlafli.pdf` |
| +1 | its complement, SRG(27,16,10,8) | 16 | `e6_complement.pdf` |

This is the bottom of the ladder and the tie-back: the 27 = W(E₆)/W(D₅) are
the **27 lines on a cubic surface**, and the ip=0 orbital is exactly the
**Schläfli graph SRG(27,10,1,5)** built in `../hecke/` from PSp(4,3). Here it
is drawn in its *natural* E₆ Coxeter plane (12-fold, rings 3 + 12 + 12), and
its Hecke algebra is strongly regular because E₆ has only three orbital
classes.

## The tower of antipodal double covers

Each +1-graph is an **antipodal distance-regular graph of diameter 3** (for E₈
and E₇) or a primitive SRG (for E₆). "Folding" by the antipode — merging each
vertex with its unique distance-3 partner — gives the quotient:

| rung | +1-graph (degree) | vertices | antipodal quotient |
|------|-------------------|---------:|--------------------|
| **E₈** | 4₂₁ / Gosset (56) | 240 | **SRG(120,56,28,24)** — the E₈ root graph, O⁺(8,2):2 |
| **E₇** | Gosset graph (27) | 56 | **K₂₈** (complete — a trivial quotient) |
| **E₆** | Schläfli (10) | 27 | — (odd order ⇒ no antipodal pairs; terminal) |

So the ladder is a uniform tower of double covers, but the quotients are only
*sometimes* famous SRGs: E₈ folds to the O⁺(8,2) graph, E₇ folds merely to the
complete graph (the Gosset graph is the textbook antipodal 2-cover of K₂₈), and
E₆ cannot fold at all. **The once-conjectured SRG(28,12,6,4) does not appear** —
E₇'s even spectrum has no room for it (see below). Verify with `fold_analysis.py`.

### The fold, read off the spectrum

The antipode σ splits each cover's adjacency into a σ-even and a σ-odd part,
and the quotient's spectrum is *exactly* the σ-even half:

- **E₈:** {56, 28, 8, −2, −4} splits as even {56, 8, −4} (= SRG(120,…)) and
  odd {28⁸, (−2)¹¹²}.
- **E₇:** {27, 9, −1, −3} splits as even {27, −1} (= K₂₈, the only two-eigenvalue
  connected graph) and odd {9⁷, (−3)²¹}.

The **σ-odd part is the reflection representation** of the Weyl group: one
eigenvalue (28, 9, −5) with multiplicity **8 → 7 → 6**, tracking the rank
descent E₈ ⊃ E₇ ⊃ E₆. The geometry — the rank-n embedding itself — is precisely
the part that folding destroys. All three graphs are integral.

## Why the drawings work

`draw_graph.py`'s philosophy — vertices on concentric rings with rotational
symmetry — is exactly the Coxeter-plane projection: project onto the eigenplane
of a Coxeter element for the eigenvalue e^{2πi/h}, and the h-fold symmetry and
the ring structure appear automatically. The ring counts are the exponents of
the group (E₈: 1,7,11,13,17,19,23,29 → 8 rings).

**Caveat (honest):** the pictures use direct root/weight computation and the
Coxeter projection, not coset enumeration of the Weyl groups. The group/subgroup
pairs above are the Hecke-algebra interpretation of the same objects. All
cardinalities, valencies, spectra, Coxeter orders (30/18/12) and ring
structures were verified numerically.

## Files

**Build pipeline** (run in order; reproducible without the session transcript):
- `ladder.py` — shared machinery: the 240 E₈ roots, simple-root/Coxeter-plane finder, projection
- `e8_mandala.py` — E₈: 240 roots, inner-product signature, Coxeter projection, 8-ring verification → `/tmp/e8_{R,P}.npy`
- `e7_build.py` — E₇ minuscule descent (roots at ip +1 to r₀) + Coxeter plane → `/tmp/e7_{W,P,r0}.npy`, verifies rings 2 + 3×18
- `e6_build.py` — E₆ minuscule descent (E₇-weights at ip +1 to w₀) + Coxeter plane → `/tmp/e6_{W,P}.npy`, verifies rings 3 + 2×12

**Analysis & write-up:**
- `fold_analysis.py` — the antipodal-fold / spectral-split story above, fully verified
- `spectral-fold.tex` / `.pdf` — note (theorems + proofs): folding = discarding the
  reflection representation. Includes the full rank-5 E₈ Bose–Mesner multiplication
  table, the P-matrix/character table, and the ℂ⁵ decomposition (irreps ↔ the
  W(E₈)-constituents of dims 1, 8, 35, 112, 84)

**Drawers** (read the `/tmp/*.npy` caches):
- `e8_tikz.py`, `e8_one.py`, `e8_bloom.py` — E₈ orbital drawings
- `e7_mandala.py`, `e7_one.py` — E₇ (56 weights) drawings
- `e6_one.py` — E₆ (27 lines) drawings
- `e8_*.pdf`, `e7_*.pdf`, `e6_*.pdf` — the rendered mandalas
