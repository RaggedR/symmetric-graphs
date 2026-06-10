# Tessellation Graphs of Hurwitz and Non-Hurwitz Surfaces

Preliminary results and conjectures from investigating the relationship
between Riemann surface symmetry and tessellation graph symmetry.

## Setup

A hyperbolic triangle group Delta(p,q,r) with 1/p + 1/q + 1/r < 1 tiles the
Poincare disk. A torsion-free normal subgroup Gamma gives a Riemann surface
X = Disk/Gamma with automorphism group G = Delta/Gamma.

For the (2,q,r) family, the surface carries dual tessellations {q,r} and {r,q}
whose 1-skeleta are two dual graphs:

- The **r-regular graph** on |G|/r vertices (stabilizer C_r)
- The **q-regular graph** on |G|/q vertices (stabilizer C_q)
- Both share |G|/2 edges (stabilizer C_2)

These are Sabidussi coset graphs G/C_r and G/C_q.

## Surfaces Studied

### Klein quartic (genus 3, Hurwitz)

- Triangle group: Delta(2,3,7)
- Automorphism group: PSL(2,7), order 168
- {7,3} tessellation: 24 heptagons, 56 vertices, 84 edges
- {3,7} tessellation: 56 triangles, 24 vertices, 84 edges
- Cubic graph: F056A (Klein graph), 56 vertices, girth 7
- 7-regular graph: 24 vertices
- Both are simple graphs, both are Sabidussi graphs of PSL(2,7)
- |Aut(F056A)| = 336 = 2 x 168 = full surface symmetry (holo + anti-holo)
- H^1 = V_3 + V_3' (two conjugate 3-dim irreps), character (6,-2,0,2,-1,-1)
- CM type: {zeta, zeta^2, zeta^4} = quadratic residues mod 7

### Macbeath surface (genus 7, Hurwitz)

- Triangle group: Delta(2,3,7)
- Automorphism group: PSL(2,8), order 504
- {7,3} tessellation: 72 heptagons, 168 vertices, 252 edges
- {3,7} tessellation: 168 triangles, 72 vertices, 252 edges
- Cubic graph (Macbeath graph): 168 vertices, girth 9
- 7-regular graph: 72 vertices, girth 3
- Both are simple graphs, both are Sabidussi graphs of PSL(2,8)
- |Aut(Macbeath graph)| = 1008 = 2 x 504 = C_2 x PSL(2,8)
- Constructed by `macbeath_graph.g`

### Bolza surface (genus 2, non-Hurwitz)

- Triangle group: Delta(2,3,8)
- Automorphism group: GL(2,3), order 48
- {8,3} tessellation: 6 octagons, 16 vertices, 24 edges
- {3,8} tessellation: 16 triangles, 6 vertices, 24 edges
- Cubic graph: Mobius-Kantor (F016), 16 vertices, girth 6
- 8-regular graph: **cannot exist as simple graph** (degree 8 > n-1 = 5)
- The {3,8} 1-skeleton is a **multigraph** on K_{3,3} with edge multiplicity 4
- GL(2,3) acts **non-faithfully** on the 6 vertices (kernel = center {+/-I} = C_2)
- Effective action: PGL(2,3) = S_4 (order 24) on 6 vertices = edges of tetrahedron
- |Aut(K_{3,3})| = 72, but surface provides only 24 effective symmetries
- S_4 does NOT act edge-transitively on K_{3,3} (|S_4|=24 does not divide 9 edges)
- Constructed by `bolza_graphs.g`

## Conjecture: Hurwitz Surfaces Fully Account for Their Graph Symmetries

**Statement.** If X is a Hurwitz surface with tessellation graphs Gamma_q and
Gamma_r from the {q,r} and {r,q} tessellations, then

    |Aut(Gamma)| = |Aut_full(X)|

where Aut_full includes both holomorphic and anti-holomorphic automorphisms.
That is, every graph automorphism is realized by a surface isometry.

For non-Hurwitz surfaces, |Aut(Gamma)| >= |Aut_full(X)| with potentially
strict inequality: the graph may have **hidden symmetries** not realizable
as surface isometries.

**Evidence.**

| Surface        | Hurwitz? | |Aut(graph)| | |Aut(surface, full)| | Match? |
|----------------|----------|-------------|---------------------|--------|
| Klein quartic  | yes      | 336         | 336                 | yes    |
| Macbeath       | yes      | 1008        | 1008                | yes    |
| Bolza          | no       | 72          | 48 (effective)      | no     |

**Status.** Three examples tested. Need more Hurwitz surfaces (PSL(2,13) at
genus 14, etc.) and more non-Hurwitz surfaces (from Delta(2,3,8), Delta(2,4,5),
etc.) to build confidence.

## Taxonomy of Failures for Non-Hurwitz Surfaces

The Bolza surface exhibits four distinct failures:

1. **Multigraph degeneracy.** The dual graph has degree > n-1, forcing multiple
   edges between the same vertex pair. Condition: |G|/r < r, i.e., |G| < r^2.

2. **Non-faithful action.** The group center is contained in the vertex
   stabilizer, so the action on vertices has nontrivial kernel.

3. **Edge non-transitivity.** The effective symmetry group is too small to act
   transitively on all edges. |G_eff| does not divide |E|.

4. **Hidden symmetries.** |Aut(graph)| > |Aut(surface)|. The graph has
   symmetries that cannot be realized as surface isometries.

**Open question.** Do all four failures always occur together? Or can a
non-Hurwitz surface fail on some but not others? Larger non-Hurwitz surfaces
(where |G| > r^2) might avoid failure 1 while still exhibiting failures 2-4.

## Foster Census Connections

Three Foster census graphs arise as tessellation 1-skeleta:

| Foster  | Name            | Vertices | Surface        | Tessellation |
|---------|-----------------|----------|----------------|--------------|
| F006    | K_{3,3}         | 6        | Bolza (dual)   | {3,8}        |
| F016    | Mobius-Kantor   | 16       | Bolza          | {8,3}        |
| F056A   | Klein graph     | 56       | Klein quartic  | {7,3}        |

The Macbeath graph (168 vertices, cubic, girth 9) is not in the Foster census
(which covers n <= 30 in this collection) but was constructed in `macbeath_graph.g`.

## Code

- `klein_quartic_cohomology.py` -- Full cohomology pipeline for Klein quartic
- `verify_cohomology.g` -- GAP verification via character table
- `triangle_group.py` -- Explicit SL(2,R) generators for Delta(2,3,7)
- `macbeath_graph.g` -- Macbeath graph and dual from PSL(2,8)
- `bolza_graphs.g` -- All Sabidussi graphs of GL(2,3)
