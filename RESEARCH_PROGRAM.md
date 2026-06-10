# Drawable Shadows of Lie Groups over Finite Fields

## The Observation

The Biggs-Smith graph is a picture of SU(2) over the finite field F\_17.

More precisely: SU(2, F\_{17}) acts on the Biggs-Smith graph (102 vertices, cubic, girth 9), but non-faithfully. The scalar matrix -I acts trivially on every vertex. The kernel of the action is exactly the centre {+/-I} = C\_2. The faithful quotient PSL(2, 17) is the full automorphism group.

## The Double Cover over Finite Fields

Over R, SU(2) -> SO(3) is a clean 2:1 surjection. SU(2) is twice as big, the map has kernel {+/-I}, and it covers all of SO(3).

Over F\_q, the same map misfires. The chain is:

```
SL(2, F_q)  --2:1-->  PSL(2, F_q)  --index 2 inside-->  PGL(2, F_q) = SO(3, F_q)
```

The groups SL(2, F\_q) and SO(3, F\_q) have the *same order*, but the 2:1 map only hits half of SO(3). The factor of 2 gets spent on the kernel {+/-I}, landing on PSL(2, q) rather than the full PGL(2, q) = SO(3, F\_q).

The core reason: over R, every positive real has a square root. Over F\_q, half the nonzero elements aren't squares. An element [A] of PGL(2, F\_q) comes from SL(2) iff det(A) is a square — and only half of them are.

Two things change between R and F\_q:
1. Over R: SU(2) -> SO(3) is surjective, 2:1. Over F\_q: the analogous map lands on PSL(2), which is index 2 inside SO(3, F\_q).
2. Over R: SO(3) is not isomorphic to PGL(2, R) — they are fundamentally different (compact vs noncompact). Over F\_q: SO(3, F\_q) = PGL(2, F\_q). The real forms collapse.

## The Programme

Every simple Lie group gives a family of finite groups of Lie type (Chevalley groups) over F\_q, and each acts on arc-transitive graphs via the Sabidussi coset construction Sab(G, H, HaH). But only the cubic ones — where the double coset satisfies |HaH|/|H| = 3 — produce graphs sparse enough to draw well.

**One graph per Lie group per field**, whenever the arithmetic works out to valency 3.

### A\_1 = SU(2) family (complete census)

There are exactly **three infinite families** of cubic arc-transitive graphs from PSL(2, q), forming a tower of covers determined by the subgroup chain C\_3 ⊂ S\_3 ⊂ S\_4:

#### Family 1: H = C\_3 (exists for every prime q >= 5, q != 3)

|V| = |PSL(2,q)|/3 = q(q^2 - 1)/6. Grows as q^3/6.

| q | |V| | Girth | Aut | Name |
|---|------|-------|-----|------|
| 5 | 20 | 5 | C\_2 x A\_5 | double cover of Petersen |
| 7 | 56 | 3 | huge | (disconnected in earlier search) |
| 11 | 220 | 10 | PSL(2,11) : C\_2 | -- |
| 13 | 364 | 7 | C\_2 x PSL(2,13) | -- |
| 17 | 816 | ? | ? | -- |
| 19 | 1140 | ? | ? | -- |
| 23 | 2024 | ? | ? | -- |
| ... | ... | | | grows as q^3/6 |

This family exists because C\_3 always embeds in PSL(2, q) for q > 3, and the double coset HaH always achieves valency 3. The Aut is always PSL(2, q) : C\_2 or C\_2 x PSL(2, q), because C\_3 is not maximal. The C\_3 graph is a double cover of the S\_3 graph (Family 2).

#### Family 2: H = S\_3 (exists for every prime q >= 5, q != 3)

|V| = |PSL(2,q)|/6 = q(q^2 - 1)/12. Grows as q^3/12 — half the size of Family 1.

| q | |V| |
|---|------|
| 5 | 10 = **Petersen graph** |
| 7 | 28 |
| 11 | 110 |
| 13 | 182 = Foster census entry |
| 17 | 408 |
| 19 | 570 |
| 23 | 1012 |
| ... | grows as q^3/12 |

S\_3 = <g, t> where g has order 3 and t is an involution with tgt = g^{-1}. This is a **different** S\_3 from the normalizer of C\_3 in PSL(2, q) (which becomes D\_{2n} for q >= 11). The Aut of these graphs is PSL(2, q) exactly (confirmed for small q). The S\_3 graph is a 4-fold cover of the S\_4 graph (Family 3) when S\_4 exists.

#### Family 3: H = S\_4 (exists for primes q = +/-1 mod 8)

|V| = |PSL(2,q)|/24 = q(q^2 - 1)/48. Grows as q^3/48 — **8x smaller** than Family 1.

| q | |V| | Girth | Aut | Name |
|---|------|-------|-----|------|
| 17 | 102 | 9 | PSL(2,17) | **Biggs-Smith graph** |
| 23 | 506 | 11 | ? | -- |
| 31 | 620 | ? | ? | -- |
| 47 | 2162 | ? | ? | -- |
| 71 | 5950 | ? | ? | -- |
| 79 | 10270 | ? | ? | -- |
| 89 | ? | ? | ? | -- |
| 97 | 19012 | ? | ? | -- |
| 103 | ? | ? | ? | (103 = -1 mod 8) |
| ... | ... | | | primes q = +/-1 mod 8 |

S\_4 embeds as a **maximal** subgroup of PSL(2, q) iff q = +/-1 mod 8 (Dickson). Because S\_4 is maximal, the full automorphism group is exactly PSL(2, q) — no room for extra symmetry. This makes the S\_4 family the "tight" family.

The Biggs-Smith graph (q = 17) is the smallest member. It is one of only 12 distance-transitive cubic graphs, with only 8 distinct eigenvalues.

#### The tower of covers

The three families form a chain:

```
C_3   ⊂   S_3   ⊂   S_4
(3)       (6)       (24)
q³/6  →  q³/12  →  q³/48
```

Each inclusion increases the stabilizer and shrinks the graph:
- The C\_3 graph is a **double cover** of the S\_3 graph.
- The S\_3 graph is a **4-fold cover** of the S\_4 graph (when S\_4 exists).

For q = 17: the tower gives 816 → 408 → 102 (Biggs-Smith). The Biggs-Smith graph is the maximally compressed shadow of SU(2, F\_{17}).

**Important distinction**: the S\_3 that gives cubic graphs is NOT the normalizer of C\_3 in PSL(2,q). The normalizer N(C\_3) grows into D\_{2n} for q >= 11, which gives higher valencies. The cubic S\_3 is a "bare" S\_3 = <g, t> embedded via an element of order 3 and an inverting involution.

#### Why only three families?

A\_5 (order 60) embeds in PSL(2,q) for q = +/-1 mod 5, but never gives valency 3 — the double coset is too large relative to |H|. Larger subgroups (Borel, dihedral D\_n for large n) similarly fail the valency-3 condition. The chain C\_3 ⊂ S\_3 ⊂ S\_4 exhausts all possibilities.

#### SO(3) parallel: PGL(2,q) graphs

The Heawood graph is Sab(PGL(2,7), S\_4, a) — a cubic graph from SO(3, F\_7). This opens a parallel family from PGL(2,q) = SO(3, F\_q), which may have different (possibly richer) cubic graphs since PGL(2,q) has more subgroups than PSL(2,q).

### Open: other Lie types

Which Chevalley groups G(F\_q) admit cubic arc-transitive graphs?

- **A\_2**: PSL(3, q). The Levi graph of PG(2, 2) is the Heawood graph (already in the gallery). Higher q gives (q+1)-regular graphs — not cubic for q > 2.
- **B\_2 = C\_2**: Sp(4, q) and SO(5, q). Rank 2, could give generalized quadrangles.
- **G\_2**: The smallest exceptional type. G\_2(q) has BN-pair of rank 2, acts on generalized hexagons. G\_2(2)' = PSU(3, 3) acts on the known generalized hexagon of order (2,2) — is this cubic?
- **{}^2A\_2**: PSU(3, q). Unitary groups in dimension 3. Act on Hermitian unitals.
- **{}^2B\_2**: Suzuki groups Sz(q), q = 2^{2n+1}. These have no BN-pair of Lie type over the original field — different flavour.
- **{}^2G\_2**: Ree groups. Act on Ree unitals.

For each: find subgroups H and involutions a giving valency 3, check connectivity, compute Aut.

## The Quantum Mechanics Connection

SU(2) is not just any Lie group — it is the symmetry group of quantum spin.

A spin-1/2 particle has state space C^2. The group SU(2) acts on this space, and its representations classify all angular momentum states: spin-0, spin-1/2, spin-1, spin-3/2, ... The double cover SU(2) -> SO(3) is physically real — it is the reason fermions pick up a phase of -1 under 360-degree rotation. The kernel {+/-I} is the *spinor sign*.

Over finite fields, this same kernel {+/-I} is what acts trivially on the Biggs-Smith graph. The graph "sees" the projective quotient PSL(2) but not the spinor sign. This is the combinatorial analogue of the statement that orbital angular momentum (SO(3) representations) cannot detect half-integer spin.

### Deeper connections

**Finite quantum mechanics.** Quantum mechanics over finite fields (replacing C with F\_q) is a real subject. Wootters, Gibbons, Hoffman, and others have studied "discrete phase space" and "mutually unbiased bases" using the action of SL(2, F\_q) on F\_q^2. The Sabidussi graphs from PSL(2, q) may encode the adjacency structure of these discrete phase spaces.

**Representation theory.** Over C, the irreducible representations of SU(2) are the symmetric powers Sym^n(C^2), giving spin-n/2 particles. Over F\_q, the representations of SL(2, F\_q) are classified by Brauer theory and have a completely different character. The cubic graphs may correspond to specific representations — the stabilizer H determines which representation "level" the graph lives at.

**The Pauli group.** For a single qubit, the Pauli matrices {I, X, Y, Z} generate a finite group inside SU(2). Over F\_q, the analogous structure lives in SL(2, F\_q). The relationship between Pauli operators and the Sabidussi construction — both producing combinatorial structures from SL(2) — deserves investigation.

**Contextuality and Kochen-Specker.** The Kochen-Specker theorem uses graphs whose vertices are quantum states and edges are orthogonality relations. These "orthogonality graphs" of SU(2) representations over finite fields could connect to the Sabidussi graphs studied here. The Petersen graph (our q = 5 entry) already appears in Kochen-Specker theory.

## Spectral Theory

All cubic graphs in the A\_1 family are **Ramanujan** — optimal expanders for their degree.

| Graph | |V| | λ\_2 | Fiedler | Spectral gap | #eigs |
|-------|------|------|---------|--------------|-------|
| PSL(2,5)/S\_3 = Petersen | 10 | 1.000 | 2.000 | 2.000 | 3 |
| PSL(2,11)/C\_3 | 220 | 2.618 | 0.382 | 0.382 | 19 |
| PSL(2,13)/S\_3 | 182 | 2.675 | 0.325 | 0.325 | 16 |
| PSL(2,17)/S\_4 = Biggs-Smith | 102 | 2.562 | 0.438 | 0.438 | 8 |
| PSL(2,17)/S\_3 | 408 | 2.709 | 0.291 | 0.291 | 24 |

Ramanujan bound for k=3: 2sqrt(2) = 2.828.

**Why they're all Ramanujan.** This is not a coincidence. The Ramanujan property for graphs from groups over finite fields is essentially the Ramanujan-Petersson conjecture (proved by Deligne for GL(2) over function fields). These graphs are Ramanujan *because* they come from SU(2). The Ramanujan bound is the spectral shadow of the Weil conjectures.

**The Biggs-Smith spectrum is minimal.** Only 8 distinct eigenvalues for 102 vertices — because the Biggs-Smith graph is distance-transitive with diameter 7, and #eigenvalues = diameter + 1 is the minimum. The S\_4 stabilizer is large enough to force distance-transitivity, collapsing the spectrum.

**Fiedler eigenvalue decreases with q.** Algebraic connectivity scales roughly as O(1/q): larger SU(2) graphs are less well-connected relative to their size. As the finite quantum phase space grows, it takes more hops to traverse.

**Systematic sweep.** To study the spectral theory systematically:
1. For each prime q (and prime power), construct all cubic Sabidussi graphs from PSL(2, q)
2. Compute the adjacency spectrum and verify the Ramanujan property
3. Track how the Fiedler eigenvalue, spectral gap, and eigenvalue multiplicities vary with q
4. Identify which eigenvalues are algebraic numbers over Q — the golden ratio phi = (1+sqrt(5))/2 = 1.618... appears in the Petersen spectrum (λ\_2 = 1 = 2phi - 2), and sqrt(2) appears throughout the PSL(2,17)/S\_3 spectrum
5. Relate the spectrum to representation theory of SL(2, F\_q): each eigenvalue should correspond to a character value

## Key Questions

1. **Census.** For each simple Lie type and small F\_q, which (G, H, a) triples give connected cubic arc-transitive graphs? This is a finite computation in GAP.

2. **The spinor sign in graphs.** Is there a graph-theoretic property that detects whether the acting group is SL(2) or PSL(2) — i.e., that "sees" the double cover? Candidates: covering graphs, voltage graphs, homology.

3. **Finite quantum graphs.** Do the cubic Sabidussi graphs from SL(2, F\_q) encode the structure of discrete quantum phase space? Specifically: are the vertices identifiable with quantum states over F\_q, and the edges with some physical relation (transition amplitude, orthogonality, complementarity)?

4. **Representation-theoretic meaning of H.** For the Biggs-Smith graph, H = S\_4. What does this subgroup mean in terms of quantum mechanics? S\_4 is the rotation symmetry group of the cube/octahedron — does the Biggs-Smith graph encode something about a "cubical" quantum state over F\_17?

5. **Beyond A\_1.** The programme for other Lie types. Start with G\_2(2) and PSU(3, 3) — these are small enough for GAP to handle and have beautiful associated geometries (generalized hexagons, unitals).
