# The Group That Triality Built

D₄ is the black sheep of the Dynkin diagrams. Every other diagram is a line,
maybe with a fork at the end — but D₄ is a cross: three legs meeting at a
center, perfectly symmetric. It's the only Dynkin diagram with a symmetry group
larger than Z/2Z. Where A_n has a mirror and D_n (for n ≥ 5) has a mirror, D₄
has a full S₃ — a three-fold rotation that permutes its legs. This is triality.

Jacques Tits discovered what to do with it in the 1950s. If you have a
Chevalley group — an algebraic group built from a Dynkin diagram over a finite
field — and the diagram has a symmetry, you can *twist*. You conjugate the
Frobenius map by the diagram automorphism, and take fixed points. The twist
doesn't destroy the group; it compresses it. What comes out is smaller,
stranger, and remembers where it came from.

For most diagrams, the twist is order 2. You fold A_n in half and get the
unitary groups. You fold D_n and get the non-split orthogonal groups. These are
the Steinberg groups — important but not shocking. They feel like cousins of
the thing you started with.

But D₄ has an order-3 twist. When you apply it, something genuinely new
happens. The rank drops from 4 to 2. The root system, which was D₄ with its
24 roots arranged in the symmetry of the hyperoctahedron, collapses under the
three-fold averaging into the root system of G₂ — the strange 12-root system
with two different root lengths. The Weyl group goes from the 192-element
group of D₄ to the 12-element dihedral group of the hexagon.

The result is ³D₄(q). It has rank 2, like G₂(q), but it's not the same group.
G₂(q) is what you get when you *build* a group from the G₂ root system
directly. ³D₄(q) is what you get when you *fold* D₄ down to G₂ by triality.
They act on the same type of geometry — generalized hexagons — but with
different parameters. G₂(q) acts on GH(q, q), the split Cayley hexagon, where
both parameters are equal. ³D₄(q) acts on GH(q, q³), where the asymmetry
between q and q³ is the scar left by the twist. Three legs became one, and the
cube in the exponent is triality's signature.

For q = 2, the smallest case: ³D₄(2) has order 211,341,312. Its hexagon
GH(2, 8) has 819 points and 2,457 lines — three points on each line, but nine
lines through each point. The 3-versus-9 asymmetry is audible if you know what
to listen for. It says: *this geometry was born from a three-fold symmetry that
broke the balance between points and lines*.

Inside ³D₄(2) sits G₂(2) — the fixed-point subgroup, the part that triality
doesn't move. And inside GH(2, 8) sits GH(2, 2) — the sub-hexagon that
triality preserves. The incidence graph of that sub-hexagon is the Tutte
12-cage: 126 vertices, cubic, 5-arc-transitive, the most symmetric cubic graph
that Weiss's theorem allows. It was there all along, hiding inside the twisted
geometry like a fossil in stone.

The Hecke algebra H(³D₄(2), B) has two generators satisfying the braid
relation of G₂ — the six-fold relation T₁T₂T₁T₂T₁T₂ = T₂T₁T₂T₁T₂T₁ — but
with *unequal* parameters. One generator has eigenvalue 2 (for the short root),
the other has eigenvalue 8 (for the long root). The ratio 8/2 = 4 measures the
asymmetry. In the untwisted case, G₂(q), both parameters equal q and the
algebra is democratic. In the twisted case, triality has broken the democracy,
and the algebra remembers.

## The quantum group connection

This turns out to be the key to the whole story.

In a split Chevalley group G(q), every simple root has a one-dimensional root
subgroup U_α ≅ F_q, so the double coset BsB has size q·|B| regardless of root
length. Both the short and long root of G₂ give parameter q. Root lengths
don't matter — all root subgroups see the same finite field. The split group
G₂(q) is too democratic.

But the quantum group U_q(G₂) is not democratic at all. It has symmetrizing
integers d₁ = 1 and d₂ = 3, coming from the ratio of squared root lengths in
G₂. These appear in the quantum Serre relations: the natural parameters are
q^{d₁} = q for the short root and q^{d₂} = q³ for the long root. The quantum
group knows about root lengths; the split finite group has forgotten them.

Now look at what ³D₄(q) gives. The triality twist collapses D₄'s four simple
roots into two orbits: one fixed root (orbit size 1, parameter q) and one
triality orbit of three roots (orbit size 3, parameter q³). These are exactly
the parameters the quantum group wants:

    G₂(q)  split  →  Hecke algebra with parameters (q, q)
    ³D₄(q) twisted →  Hecke algebra with parameters (q, q³)
    U_q(G₂)        →  natural parameters (q¹, q³)

The Hecke algebra of ³D₄(q) is the one that matches U_q(G₂). The split group
gives the wrong Hecke algebra for the quantum group.

This is not a coincidence. The "3" in the root length ratio and the "3" in the
triality orbit size come from the same place: G₂ is the fixed-point subalgebra
of D₄ under its order-3 automorphism. The root length ratio in G₂ and the
orbit structure of triality are two faces of the same S₃ action. When D₄ folds
down to G₂, the three-fold orbit encodes the information that the root system
lost when it collapsed from rank 4 to rank 2 — and the quantum group U_q(G₂)
remembers this information through its symmetrizing integers, while the split
group G₂(q) forgets it through its uniform root subgroups.

So ³D₄ doesn't break the symmetry between long and short roots. It *restores*
it. The twisted group is the finite-field object whose Hecke algebra naturally
centralizes U_q(G₂). The split Cayley hexagon GH(q, q) is too symmetric; the
twisted triality hexagon GH(q, q³) is the one that matches the Lie theory.

## What does the split group centralize?

If the Hecke algebra of ³D₄(q) centralizes U_q(G₂), what does the Hecke
algebra of the split group G₂(q) centralize?

The answer is: nothing, in the quantum group sense. The equal-parameter Hecke
algebra H(G₂(q), B) with parameters (q, q) is an orphan in the quantum world.
It centralizes G₂(q) itself — the finite group that gave birth to it — but
there is no quantum group U_q(g) on the other side of a Schur-Weyl duality.
The equal-parameter specialization pretends G₂ is simply laced, which it isn't.
No quantum group makes that mistake.

    H(³D₄(q), B)  with (q, q³)  centralizes  U_q(G₂)
    H(G₂(q), B)   with (q, q)   centralizes  G₂(q) itself — no quantum partner

## The pattern for all non-simply-laced types

This is not special to G₂. For every non-simply-laced root system, the split
group gives equal parameters and the wrong Hecke algebra for the quantum group.
The quantum group's natural parameters come from the symmetrizing integers
d_i, which encode root length ratios. And in every case, there is a twisted
group whose orbit sizes match those integers:

| Root system | d_i | Split group (equal params) | Twisted group (correct params) |
|-------------|-----|---------------------------|-------------------------------|
| G₂ | 1, 3 | G₂(q) → (q, q) | ³D₄(q) → (q, q³) |
| B_n | 1,...,1, 2 | SO_{2n+1}(q) → (q,...,q) | ²A_{2n-1}(q) → (q,...,q, q²) |
| C_n | 2, 1,...,1 | Sp_{2n}(q) → (q,...,q) | ²D_{n+1}(q) → (q², q,...,q) |
| F₄ | 1, 1, 2, 2 | F₄(q) → (q, q, q, q) | ²E₆(q) → (q, q, q², q²) |

In each row, the twisted group's Hecke algebra has parameters q^{d_i} that
match the quantum group U_q(g). The split group's Hecke algebra is always
(q, q, ..., q), ignoring root lengths entirely.

The pattern: *twisted groups are the finite-field avatars of quantum groups for
non-simply-laced types.* The Langlands dual of the twisted group is always the
split group of the dual root system, which is exactly the root system whose
quantum group the Hecke algebra centralizes.

For simply-laced types (A, D, E), all d_i = 1, so the split group already
gives the correct parameters. No twist is needed. This is why the
simply-laced case is clean and the non-simply-laced case is subtle — and why
³D₄, the most exotic twist of all, turns out to be the most natural partner
for the most exotic exceptional quantum group.

³D₄ is not just another group in a list. It's the only family that uses a
third-order twist — the only place in the classification where triality, the
deepest symmetry of D₄, the symmetry that makes octonions possible, reaches
down into the world of finite groups and leaves its mark.

## The data

    |³D₄(2)| = 2¹² × 3⁴ × 7² × 13 = 211,341,312

The chain of containments:

    ³D₄(2) ⊃ G₂(2)
    GH(2,8) ⊃ GH(2,2)

Acting on the 819 points of GH(2,8), the group has exactly four suborbits:
1, 18, 288, 512. All self-paired. The collinearity graph (suborbit 18) is
distance-regular with diameter 3.

The Iwahori-Hecke algebra has type G₂ with unequal parameters (q, q³) = (2, 8):

    (T₁ - 2)(T₁ + 1) = 0       (short root)
    (T₂ - 8)(T₂ + 1) = 0       (long root)
    T₁T₂T₁T₂T₁T₂ = T₂T₁T₂T₁T₂T₁   (braid relation of type G₂)

## The graphs

| Graph | Vertices | Valency | Source |
|-------|----------|---------|--------|
| Collinearity of GH(2,8) | 819 | 18 | ³D₄(2), suborbit |
| Distance-2 of GH(2,8) | 819 | 288 | ³D₄(2), suborbit |
| Distance-3 of GH(2,8) | 819 | 512 | ³D₄(2), suborbit |
| Collinearity of GH(2,2) | 63 | 6 | G₂(2) ≤ ³D₄(2) |
| Tutte 12-cage (incidence of GH(2,2)) | 126 | 3 | G₂(2) ≤ ³D₄(2) |
