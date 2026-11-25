# Handover: the E₆ ⊂ E₇ ⊂ E₈ Hecke-algebra / mandala work

*Written 2026-11-25, for the next instance. Companion to `README.md` (overview)
and `~/.claude/dreams/2026-11-25-birthday-e8-ladder.md` (the texture). The work
was done in conversation while Robin showed it to his friend **Miriam**, who is
interested in Hecke algebras — hence the name.*

## What this is

A trilogy of coset Hecke algebras / association schemes from the exceptional
Weyl-group chain **W(E₈) ⊃ W(E₇) ⊃ W(E₆) ⊃ W(D₅)**, each drawn as a
Coxeter-plane mandala (the `draw_graph.py` ring philosophy = projecting onto a
Coxeter element's eigenplane). It extends the graphviz programme up the
exceptional series and closes it onto the **Schläfli graph** already in
`../hecke/` (the E₆ rung *is* the 27 lines).

**Construction (minuscule descent):** the vertices of each rung are the vectors
of the rung above lying at inner product **+1** (60°) to a fixed vector.
Start from the 240 E₈ roots; the 56 E₇ weights are the roots at ip +1 to a
fixed root r₀; the 27 E₆ weights are the E₇ weights at ip +1 to a fixed weight
w₀. All of it is direct linear algebra in ℝ⁸.

## The three rungs — all numerically verified

| rung | G / H | \|G/H\| | Coxeter h | rings | scheme rank |
|---|---|--:|--:|---|--:|
| E₈ | W(E₈)/W(E₇) | 240 | 30 | 8 × 30 | 5 |
| E₇ | W(E₇)/W(E₆) | 56 | 18 | 2 + 3×18 | 4 |
| E₆ | W(E₆)/W(D₅) | 27 | 12 | 3 + 2×12 | 3 |

Indices check: 696729600/2903040 = 240, 2903040/51840 = 56, 51840/1920 = 27.
|W(E₈)|=696729600, |W(E₇)|=2903040, |W(E₆)|=51840, |W(D₅)|=1920.
The scheme rank **descends 5 → 4 → 3** down the ladder.

### E₈ (240 roots) — rank 5, suborbits {1,1,56,56,126}

Orbital graphs by inner product, valencies (1, 56, 126, 56, 1):
- **A = ip +1** (56-reg): the mandala, the 4₂₁/Gosset edge graph. Spectrum
  **{56¹, 28⁸, 8³⁵, (−2)¹¹², (−4)⁸⁴}** (5 eigenvalues ⇒ A generates the algebra).
- **B = ip 0** (126-reg): "the bloom" (orthogonal roots).
- **C = ip −1** (56-reg): "the weave" = A under the antipode.
- **D = ip −2** (1-reg): 120·K₂, the antipodal matching.
- Verified structure constant: **A² = 56·I + 27·A + 12·B + 1·C.**

**The P-matrix** (rows = eigenspaces by A-eigenvalue, cols = I A B C D):
```
 mult |  I    A    B    C    D
   1  |  1   56  126   56    1
   8  |  1   28    0  -28   -1
  35  |  1    8  -18    8    1
 112  |  1   -2    0    2   -1
  84  |  1   -4    6   -4    1
```
Patterns: **C = A ⊙ D** (the two 56-graphs are antipodal twins); **B is even**
(zero on the D=−1 eigenspaces); **D is a ℤ/2 grading** (240 = 120₊ ⊕ 120₋);
multiplicities 1,8,35,84,112 are the W(E₈)-irrep dims (8 = reflection rep).

**The headline result:** folding by the antipode (120 positive roots, adjacent
iff ip = ±1) gives **SRG(120, 56, 28, 24)** — verified: deg 56, λ=28, μ=24,
spectrum {56, 8³⁵, (−4)⁸⁴}. This is the **E₈ root graph**, Aut = **O⁺(8,2):2**
(= W(E₈)/{±1}). It is a celebrated graph (Q⁺(7,2) ovoids/spreads, S(3,4,8)
overlarge sets, Johnson-graph distance graphs; famous as a **quantum-isomorphism**
example, *Algebraic Combinatorics* alco.335 / arXiv:2209.14906; recent S₇
construction arXiv:2606.10652). Parameter set is NOT unique — 3 primitive SRGs
share it; the E₈ one is the O⁺(8,2):2 member. Its complement is SRG(120,63,30,36),
and **B = 2 × that complement** (eigenvalues 126, 6, −18 = 2·(63, 3, −9)).

### E₇ (56 minuscule weights) — rank 4, suborbits {1,1,27,27}

- **Gosset graph** (ip +1, 27-reg): spectrum **{27¹, 9⁷, (−1)²⁷, (−3)²¹}**.
- partner (ip 0, 27-reg).
- **28·K₂** (ip −1, 1-reg): each weight's unique opposite — **the 28 bitangents /
  double-six pairs.**

### E₆ (27 lines) — rank 3 ⇒ STRONGLY REGULAR

Inner-product classes are only {2, 1, 0} (no −1), so the scheme is rank 3:
- **Schläfli graph** (ip 0, 10-reg): **SRG(27, 10, 1, 5)** — the same graph built
  from PSp(4,3) in `../hecke/schlafli.*`, here drawn in its natural E₆ Coxeter
  plane.
- complement (ip 1, 16-reg): SRG(27, 16, 10, 8).

## Files (all in `~/git/graphviz/e8/`)

- `e8_mandala.py` — 240 roots, ip signature, Coxeter projection, 8-ring check. **Saves `/tmp/e8_R.npy`, `/tmp/e8_P.npy`.**
- `e8_one.py` — generic E₈ orbital drawer (args: ip fname op lw title).
- `e8_tikz.py` — the main mandala; `e8_bloom.py` — the 126-reg via matplotlib (too many edges for TeX).
- `e7_mandala.py`, `e7_one.py` — E₇; depend on `/tmp/e7_W.npy`, `/tmp/e7_P.npy`.
- `e6_one.py` — E₆; depends on `/tmp/e6_W.npy`, `/tmp/e6_P.npy`.
- PDFs: `e8_{id,mandala,ip0,ipm1,ipm2}.pdf`, `e7_{mandala,partner,bitangents}.pdf`, `e6_{schlafli,complement}.pdf`.
- `README.md` — the overview with the ladder table.

## ✅ Reproduction — now fully committed (was the old TODO)

The scripts cache numpy arrays in `/tmp`, which does NOT survive a reboot. The
build steps are now committed (the old TODO — lifting the inline E₇/E₆ blocks
into files — is **done**). Regenerate from scratch, in order (shared env
`source ~/pygame_env/bin/activate`):
1. `python e8_mandala.py` → `/tmp/e8_{R,P}.npy`.
2. `python e7_build.py` → reads e8 roots (or rebuilds), minuscule descent +
   E₇ Coxeter plane → `/tmp/e7_{W,P,r0}.npy`; asserts rings 2 + 3×18.
3. `python e6_build.py` → reads `/tmp/e7_W.npy`, descent + E₆ Coxeter plane →
   `/tmp/e6_{W,P}.npy`; asserts rings 3 + 2×12.
Shared machinery (roots, simple-root finder, Coxeter plane) is in `ladder.py`.
Verified 2026-11-25: the rebuilt weight-sets equal the originals exactly, and
the drawers reproduce the right edge counts (E₇ Gosset 756, E₆ Schläfli 135).
Descent vectors: **r₀ = e₁+e₂**, **w₀ = e₁+e₃**.

## Caveats (read before trusting)

- Drawings use **direct root/weight computation + Coxeter projection, NOT coset
  enumeration** of the Weyl groups. The G/H pairs are the Hecke-algebra
  *interpretation* of the same objects.
- The E₆ roots must be found as "E₈ roots orthogonal to **both** r₀ and w₀"
  (gives exactly 72). An earlier ad-hoc "reflections preserving the 27-set"
  test overcounted to 74 / 7 simple — do not use it.
- These Hecke algebras are **commutative (Gelfand pairs)** — Bose–Mesner
  algebras. That is *different* from Robin's main `H(PSL(2,q), C₃)` work, which
  is deliberately **non-Gelfand**. The connection is thematic (orbital Hecke
  algebras from groups), not that these are the same objects.
- Recurring gotcha: `E_8` in a LaTeX node needs math mode (`$E_8$`) or the
  compile dies with "Missing $ inserted".

## ✅ Resolved 2026-11-25 (thread #1 — the fold question)

**Done, and the conjecture was corrected.** The ladder IS a uniform tower of
antipodal double covers, but the quotients are only *sometimes* named SRGs:
- E₈ 4₂₁ graph (240v, diam 3) folds to **SRG(120,56,28,24)** — the O⁺(8,2):2
  root graph. ✓
- E₇ Gosset graph (56v, diam 3) folds only to **K₂₈** (complete) — **NOT**
  SRG(28,12,6,4). The conjectured bitangents-SRG is not there.
- E₆ Schläfli (27v, odd order) has no antipodal pairs — **terminal**. ✓

Spectral punchline (answers "anything interesting about the spectra?"): the
quotient's spectrum is exactly the σ-even half of the cover's; the **σ-odd part
is the reflection representation**, one eigenvalue (28, 9, −5) with multiplicity
**8 → 7 → 6** = the rank descent. All committed in `fold_analysis.py`; written
up in `README.md` (thread #5 — the SRG identification — also done there).

## Open threads / next steps (ranked)

1. **Full multiplication tables.** The complete rank-5 E₈ table (all
   intersection numbers for {I,A,B,C,D}) is now written out in
   `spectral-fold.tex` §1 and verified by valency row-sums. Note D acts as the
   involution I↔D, A↔C, B↦B. STILL TODO: the E₇ (rank-4) and E₆ (rank-3) tables.
2. **The quantum-isomorphism angle** — the E₈ graph is quantum- but not
   classically isomorphic to a Brouwer–Ivanov–Klin graph. Worth understanding
   as a connection to Robin's interests.
3. **A matched-scale triptych** (E₆, E₇, E₈ mandalas side by side) for a print,
   and the remaining sibling faces if wanted.
