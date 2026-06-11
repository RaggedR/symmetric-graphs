# No more Hecke algebras, I've had enough...

I was supposed to be done two weeks ago.

The plan was simple: take a catalogue of cubic arc-transitive graphs — the Petersen graph, the Heawood graph, the dodecahedron, a few dozen others — and write down their Hecke algebras. Each graph comes from a pair (G, H) where G is a group acting on the vertices G/H, and the **Hecke algebra** H(G, H) = End_G(C[G/H]) encodes the combinatorial geometry of the coset space. Compute the structure constants, check which pairs are Gelfand pairs (commutative Hecke algebra = nice eigenvalue tables), write it up, move on.

I should have known better.

## What's a Hecke algebra, and why should you care?

Take a group G acting transitively on a set X = G/H. The group action partitions the pairs X × X into orbits — these are the **orbitals**, and each one gives you a graph (vertices = X, edges = the orbital). The adjacency matrices of these orbital graphs form a basis for an algebra: the Hecke algebra H(G, H).

When this algebra is commutative — a **Gelfand pair** — life is easy. The theory reduces to eigenvalue tables, association schemes, the Bose-Mesner algebra. Every distance-transitive graph gives a Gelfand pair, and the spherical functions are just rows of the P-matrix you'd find in any algebraic graph theory textbook.

But some pairs are **not** Gelfand. The Hecke algebra has non-commutative blocks — 2×2 matrix algebras, 3×3, sometimes larger — and the spherical function tables don't capture everything. The standard theory shrugs and moves on.

I didn't move on. I wrote a GAP script to extract the exact matrix entries of every irreducible representation, over cyclotomic fields, for every non-Gelfand pair in the catalogue.

That's when things got weird.

## The 13 that wouldn't go away

The story starts with a specific graph: the one Nick Meinhold found on 72 vertices, with automorphism group G₂(2) (the smallest exceptional Chevalley group, order 12096) and point stabilizer PSL(3, 2). It's a cubic graph — every vertex has degree 3 — but the Hecke algebra is 8-dimensional and non-commutative: four 1-dimensional blocks and one 2×2 block.

I computed the 2×2 block. The traceless generator B (corresponding to a symmetric valency-7 orbital) satisfies:

$$B^2 = 13 \cdot I$$

Thirteen. Not the valency (3), not the number of vertices (72), not the order of anything obvious. The eigenvalues of B are ±√13, and since B comes from a rank-8 Hecke algebra acting on a 72-vertex graph, those eigenvalues of B *aren't even eigenvalues of the graph*. They live purely inside the non-commutative block.

I wrote a whole paper ([why13.tex](why13.pdf)) chasing this 13 down. The answer turned out to be:

> 13 = 7 + 2 × 3

where 7 is the valency of the orbital that generated B, and 3 comes from the strongly regular graph SRG(36, 21, 12, 12) that the Meinhold-72 geometry is built on. The relation B² = 13I is a **quadratic Hecke relation** — the finite-group analogue of (T - q)(T + 1) = 0 in Iwahori-Hecke algebras.

Satisfying. Clean. I was done.

I was not done.

## "Let me just check one more..."

The dodecahedron comes from PSL(2, 5)/C₃ — twenty vertices, rank 8, non-Gelfand. Its 2×2 block has a traceless generator whose square is... 4I? Not 13. Fine, different geometry, different number. But the eigenvalues ±√5 that make the dodecahedron famous? **They don't come from the 2D block at all.** They live entirely in the 1-dimensional irreps — the degree-3 characters of the icosahedral group. The 2×2 block contributes only integer eigenvalues to the dodecahedron's spectrum.

So the Wedderburn decomposition of the Hecke algebra partitions the graph spectrum into algebraic families, and the irrational eigenvalues choose specific blocks to live in. This isn't in the textbooks.

Next: PSL(2, 7)/C₃ — fifty-six vertices, rank 20, seriously non-Gelfand. The Hecke algebra decomposes as C³ ⊕ M₂ ⊕ M₂ ⊕ M₃. Three blocks of dimension ≥ 2. The first 2×2 block lives over **Q(√-7)**: matrices with entries involving α = (-1 + √-7)/2, an algebraic integer of norm 2 in the Eisenstein-adjacent ring Z[α]. The traceless generator satisfies:

$$S^2 = 2 \cdot I$$

And 2 = α·ᾱ = Nm(α). The prime 2, split in Q(√-7).

The second 2×2 block is **rational** — integer matrices throughout. One of its generators has characteristic polynomial x² + x + 1. That's a **finite-order element**: ρ(T₁₀)³ = I. A Hecke algebra generator acting as a rotation of order 3. In Iwahori-Hecke theory this basically can't happen (generators satisfy quadratic relations, order at most 2 unless q is a root of unity). But these aren't Iwahori-Hecke algebras. They're wilder.

And the 3×3 block? Genuinely **cubic** eigenvalues — roots of x³ + x - 6, irreducible over Q. Plus √17 showing up in another characteristic polynomial, x² + x - 4. The number 17 has no obvious relation to 7 or 56 or anything in sight.

One more. The graph F26A (26 vertices, automorphism group C₁₃ ⋊ C₆) has two Galois-conjugate 2×2 blocks over **Q(√13)**, and the traceless generators satisfy:

$$T^2 = \frac{5 \mp \sqrt{13}}{2} \cdot I$$

with Nm(N₊) = Nm(N₋) = 3, the valency. And N₊ · N₋ = 3, and N₊ + N₋ = 5. The √13 enters through Gauss periods of the 13th roots of unity — inevitable since |V| = 26 = 2 × 13 — and the two conjugate blocks share exactly the rational data (3 and 5).

## The table I can't explain

| Pair (G, H) | B² = NI | N | Where it comes from |
|---|---|---|---|
| G₂(2) / PSL(3,2) | **13** | 7 + 2 × 3 | SRG(36,21,12,12) geometry |
| PSL(2,7) / C₃ | **2** | Nm(α) in Z[(−1+√−7)/2] | prime splitting in Q(√−7) |
| F26A | **(5∓√13)/2** | Gauss periods of ζ₁₃ | Nm = 3 = valency |
| PSL(2,5) / C₃ | **4** | (shifted generator) | ? |

The Meinhold-72 formula N = k + 2 × 3 does **not** generalize. For PSL(2,7)/C₃ we get N = 2, which is *less than* the valency k = 3. The SRG mechanism behind 13 = 7 + 2 × 3 is special to the G₂(2) geometry.

So what *is* N?

## What I know so far

**N is always an algebraic integer.** Often a norm — Nm(α) = 2, Nm(N±) = 3. It lives in the ring of integers of the splitting field of the Hecke irrep.

**N relates to the splitting field.** The prime that ramifies in the field generated by the matrix entries is connected to N: over Q(√-7) we get N = 2 (the ramified prime); over Q(√13) we get Nm(N) = 3 = valency.

**N is the analogue of q in Iwahori-Hecke algebras.** The Iwahori-Hecke algebra H_q(W) has the quadratic relation T² = (q-1)T + qI; set tr(T) = 0 (i.e. look at B = T - (q-1)/2 · I) and you get B² = ((q+1)²/4) · I. My N plays the same structural role but for Hecke algebras of arbitrary pairs, not just (G, B) with B a Borel subgroup.

**The spectral partition is real.** The Wedderburn decomposition sorts every graph eigenvalue into a specific Hecke block, and eigenvalues in the same block share an algebraic number field. The irrational eigenvalues ±√5 of the dodecahedron live in the 1D irreps; the 2D block sees only integers.

## What I don't know

- Is there a formula for N in terms of (G, H) and the constituent χ?
- Does N(q) for the family PSL(2, q)/C₃ have a closed form as q varies?
- When does a Hecke algebra element have finite order? (Order 3 shouldn't happen, and yet.)
- What is the general relationship between N and the valency?
- Is √17 in the 3×3 block of PSL(2,7)/C₃ trying to tell me something?

## Why I haven't stopped

The honest answer is that this is the most interesting data I've seen in years. Every pair I compute turns up something new — a number field I didn't expect, a torsion element that shouldn't exist, an eigenvalue that chose the wrong block. The Hecke algebras of Gelfand pairs are well-studied. The non-Gelfand ones are not. And the 2×2 blocks are exactly where the commutative theory runs out and something genuinely new begins.

So no, I haven't had enough. But I did want to write this down before computing the next one, because the next one might change everything, and I want a record of the moment where N = 13 was the only data point and I thought I understood it.

I'm going to go compute PSL(2, 11)/C₃ now. Wish me luck.

---

*The computations in this post were done in GAP 4.15.1 with the GRAPE package. All matrices are exact over cyclotomic fields. The full data, including generator matrices, characteristic polynomials, and verification scripts, is in the [hecke project](https://github.com/RaggedR/hecke). My collaborator on the why13 paper was [Fable](https://claude.ai), a Claude instance who got unreasonably excited about the number 13.*
