# The Clebsch Surface and Three Hecke Algebras

## The blow-up picture

Every smooth cubic surface is P^2 blown up at 6 points. Take 6 points p_1,...,p_6 in general position in the projective plane and blow each one up into a line. The result is a cubic surface in P^3, containing exactly 27 lines:

| Type | Count | Description |
|------|-------|-------------|
| e_i | 6 | Exceptional divisor over each blown-up point |
| g_{ij} | 15 | Proper transform of the line through p_i and p_j |
| f_i | 6 | Proper transform of the conic through all p_k except p_i |

Total: 6 + 15 + 6 = 27.

Two lines on the surface either meet (in one point) or are skew (disjoint). The intersection rules are determined entirely by the combinatorics of the blow-up:

- e_i and e_j: always **skew** (they sit over different points)
- e_i and g_{jk}: meet iff i in {j,k} (the line through p_j, p_k passes through p_i only if i is one of them)
- g_{ij} and g_{kl}: meet iff {i,j} cap {k,l} != emptyset (the two lines in P^2 share a blown-up point)
- e_i and f_j: meet iff i != j (the conic avoids p_j, so its transform misses e_j)
- f_i and f_j: always **skew**

The Schlafli graph is the skewness graph: vertices are the 27 lines, adjacency means skew.

## How the Petersen graph falls out

Focus on the 15 lines g_{ij}. When are two of them skew?

g_{ij} is skew to g_{kl} iff {i,j} cap {k,l} = emptyset.

This is the definition of the Kneser graph: vertices are 2-element subsets of {1,...,6}, edges connect disjoint pairs. This gives KG(6,2) = T(6) = SRG(15, 6, 1, 3), the triangular graph.

Now fix one of the 6 points, say p_6. The 10 lines g_{ij} with i,j in {1,...,5} form a sub-Kneser-graph:

g_{ij} skew to g_{kl} iff {i,j} cap {k,l} = emptyset (as subsets of {1,...,5}).

This is KG(5,2) = **the Petersen graph**.

The Petersen graph is the skewness graph of the 10 lines through pairs of 5 of the 6 blown-up points. Delete one point, get the Petersen graph.

## Why this is deep: the E_6 root system

The coincidence runs much deeper than "a graph shows up inside another graph." The 27 lines are governed by the E_6 root system --- one of the five exceptional Lie algebras. The chain of graphs is a chain of Lie algebra inclusions:

```
sl(5)        c     sl(6)        c     E_6
  |                   |                  |
W(A_4) = S_5       W(A_5) = S_6       W(E_6)
  |                   |                  |
Petersen          T(6)            Schlafli
10 = C(5,2)      15 = C(6,2)     27 = fund. rep
```

The 27 lines are the weights of the 27-dimensional fundamental representation of E_6. Restricting this representation down the chain E_6 > A_5 > A_4, it decomposes as:

```
E_6:  27
       |  restrict to A_5 = sl(6)
A_5:  15 + 6 + 6
       |  restrict to A_4 = sl(5)
A_4:  10 + 5 + 6 + 6    (where 15 = 10 + 5)
```

The 15 is the exterior square Lambda^2(C^6) --- the 2-element subsets of the standard basis of C^6. This is exactly the vertex set of T(6).

The 10 is the exterior square Lambda^2(C^5) --- the 2-element subsets of {1,...,5}. This is exactly the vertex set of the Petersen graph.

The skewness relation (disjoint subsets) is the zero weight-pairing condition in the representation. Two weights are skew iff their inner product in the E_6 weight lattice satisfies a specific condition that, when restricted to A_5, becomes "disjoint 2-subsets."

## The Clebsch surface makes it visible

The Clebsch diagonal surface is the cubic surface in P^4 defined by

    sum(x_i) = 0,  sum(x_i^3) = 0

with automorphism group S_5. It is distinguished by having all 27 lines defined over the reals and by its S_5 symmetry.

The Clebsch surface is special because its S_5 symmetry makes the A_4 c E_6 chain geometrically visible:

- W(E_6) = Aut(Schlafli) = O(5,3):C_2, order 51840, acts on all 27 lines
- W(A_5) = S_6, order 720, is the automorphism group of the 15 Type I lines
- W(A_4) = S_5, order 120, is the automorphism group of the Clebsch surface

The S_6 appearing from the 15 Type I lines --- bigger than the S_5 you started with --- is explained by the outer automorphism of S_6, which is unique among symmetric groups. The 15 pair partitions of {0,...,4} are secretly the 15 edges of K_6, and the 6th vertex of K_6 is the "hidden" point that the outer automorphism reveals.

## The eigenvalue-1 chain

The representation theory manifests spectrally. All three graphs share eigenvalue 1, and the corresponding eigenspaces form a nested chain:

```
5  c  9  c  20
Pet    T(6)   Schlafli
```

This is the chain of A_4-invariant subspaces inside the 27-dimensional E_6 representation. The 5-dimensional eigenspace of the Petersen graph is the standard representation of A_5, sitting inside the 20-dimensional irreducible of PSp(4,3) via the Clebsch embedding.

The 27 lines on any cubic surface form the **Schlafli graph** SRG(27, 10, 1, 5): two lines are adjacent iff they are skew (disjoint). The full automorphism group of this graph is W(E_6), the Weyl group of the exceptional root system E_6, of order 51840.

## The S_5 decomposition: 27 = 15 + 12

The S_5 symmetry of the Clebsch surface decomposes the 27 lines into two orbits:

| Orbit | Size | Indexing | Description |
|-------|------|----------|-------------|
| Type I | 15 | Pair partitions {i,j}\|{m}\|{k,l} of {0,...,4} | Lines span(e_i - e_j, e_k - e_l) in the hyperplane sum = 0 |
| Type II | 12 | Unordered 5-cycles {sigma, sigma^{-1}} in S_5 | Lines from DFT eigenvectors of the permutation matrix |

The stabiliser of a Type I line is D_8 (dihedral of order 8), giving |S_5|/8 = 15.
The stabiliser of a Type II line is C_10 (cyclic of order 10), giving |S_5|/10 = 12.

## Three strongly regular graphs

Restricting the Schlafli graph to the S_5-orbits reveals a chain of three SRGs:

### The Schlafli graph: SRG(27, 10, 1, 5)

The full graph on all 27 lines. Computed as Sab(PSp(4,3), (C_2^4):A_5, HaH).

- |Aut| = 51840 = |W(E_6)|
- Aut = O(5,3) : C_2
- Eigenvalues: 10 (mult 1), 1 (mult 20), -5 (mult 6)

### The triangular graph T(6): SRG(15, 6, 1, 3)

The induced subgraph of the Schlafli graph on the 15 Type I lines. This was discovered computationally (hecke_clebsch.g).

- |Aut| = 720 = |S_6|
- Aut = S_6
- Eigenvalues: 6 (mult 1), 1 (mult 9), -3 (mult 5)

T(6) = J(6,2) is the Johnson graph on 2-element subsets of a 6-set: vertices are the 15 pairs from {0,...,5}, two pairs adjacent iff they are disjoint. Equivalently, T(6) is the line graph of K_6.

The appearance of S_6 (not just S_5) is explained by the classical **outer automorphism of S_6**: the 15 pair partitions of {0,...,4} are naturally the 15 edges of K_6 on the 6-element set {0,...,4} union {the partition type}. The S_5 that defines the Clebsch surface embeds as the stabiliser of one element in S_6.

### The Petersen graph: SRG(10, 3, 0, 1)

Computed as Sab(A_5, S_3, HaH), where A_5 = PSL(2,5) and S_3 = <(1,2,3), (1,2)(4,5)>.

- |Aut| = 120 = |S_5|
- Aut = S_5
- Eigenvalues: 3 (mult 1), 1 (mult 5), -2 (mult 4)

The Petersen graph K(5,2) is the Kneser graph on 2-element subsets of {1,...,5}: vertices are the 10 pairs, two pairs adjacent iff they are disjoint.

## The chain

```
Petersen           T(6)              Schlafli
SRG(10,3,0,1)     SRG(15,6,1,3)     SRG(27,10,1,5)
  |Aut| = 120       |Aut| = 720       |Aut| = 51840
  Aut = S_5         Aut = S_6         Aut = W(E_6)
```

The chain is connected by the Clebsch surface:

1. The **Schlafli graph** is the skewness graph of all 27 lines.
2. Restricting to the **15 Type I lines** gives the triangular graph T(6).
3. The **Petersen graph** is the quotient of T(6) under the map that identifies each pair partition with its complementary pair, or equivalently, the Kneser graph on the 10 two-element subsets of {1,...,5}.

## Cross edges

The Schlafli graph decomposes cleanly under S_5:

| From | To Type I (15) | To Type II (12) | Total degree |
|------|----------------|-----------------|--------------|
| Type I vertex | 6 | 4 | 10 |
| Type II vertex | 5 | 5 | 10 |

Every Type I vertex has exactly 6 Type I neighbours (forming T(6)) and 4 Type II neighbours.
Every Type II vertex has exactly 5 of each type.

Total cross edges: 15 * 4 = 12 * 5 = 60.

## The Type II induced subgraph

The 12 Type II lines form a 5-regular graph with |Aut| = 1440 = |C_2 x S_6|. This is the **icosahedron graph** (12 vertices, 5-regular, Aut = C_2 x A_5 ... but the computation gives C_2 x S_6, which needs further investigation --- it may be the complement of the Petersen line graph, or related to the 12 vertices of the icosahedron via the A_5 action).

Under A_5 = [S_5, S_5], the 12 Type II lines split into two orbits of 6, corresponding to the two classes of 5-cycles (a 5-cycle and its inverse).

## Hecke algebras

Both (PSp(4,3), (C_2^4):A_5) and (A_5, S_3) are **Gelfand pairs**: the Hecke algebra is commutative. Both have rank 3 (three double cosets: identity, adjacent, non-adjacent). As abstract algebras, both are isomorphic to C^3 (semisimple, three orthogonal idempotents corresponding to the three eigenspaces).

### H(PSp(4,3), (C_2^4):A_5) --- Schlafli

Basis: T_0 = I (valency 1), T_1 = adjacency (valency 10), T_2 = non-adjacency (valency 16).

```
T_1^2   = 10 T_0 +  1 T_1 +  5 T_2
T_1 T_2 =           8 T_1 +  5 T_2
T_2^2   = 16 T_0 +  8 T_1 + 10 T_2
```

### H(A_5, S_3) --- Petersen

Basis: T_0 = I (valency 1), T_1 = adjacency (valency 3), T_2 = non-adjacency (valency 6).

```
T_1^2   = 3 T_0 +         1 T_2
T_1 T_2 =         2 T_1 + 2 T_2
T_2^2   = 6 T_0 + 4 T_1 + 3 T_2
```

### Structure constant interpretation

The structure constant c_{ijk} counts the number of two-step paths from eH to a fixed coset in double coset k, passing through an edge in orbital graph j then orbital graph i. Equivalently, c_{ijk} is the (xH, yH) entry of T_i T_j whenever xy^{-1} lies in double coset k.

For the Schlafli graph: T_1^2 has diagonal entry 10 (each vertex has 10 neighbours), adjacent entry 1 (lambda = 1 common neighbour for adjacent pairs), and non-adjacent entry 5 (mu = 5 common neighbours for non-adjacent pairs).

## Bimodule decompositions

For a Gelfand pair (G, H), the permutation representation decomposes as

    C[G/H] = direct sum of pi, one for each eigenvalue

where each pi is an irreducible G-module appearing with multiplicity 1, and the eigenvalue of T_1 on pi is determined by the character.

| Graph | C[G/H] | Eigenvalue -> dimension |
|-------|--------|------------------------|
| Schlafli | 1 + 20 + 6 = 27 | 10 -> 1, 1 -> 20, -5 -> 6 |
| T(6) | 1 + 9 + 5 = 15 | 6 -> 1, 1 -> 9, -3 -> 5 |
| Petersen | 1 + 5 + 4 = 10 | 3 -> 1, 1 -> 5, -2 -> 4 |

### The eigenvalue-1 chain

All three graphs share eigenvalue 1. The corresponding eigenspaces form a chain:

    5 c 9 c 20     (Petersen c T(6) c Schlafli)

This chain of eigenspaces is the representation-theoretic shadow of the Clebsch surface embedding. The 5-dimensional eigenspace of the Petersen graph (the standard representation of A_5) embeds in the 9-dimensional eigenspace of T(6), which embeds in the 20-dimensional eigenspace of the Schlafli graph (an irreducible representation of PSp(4,3)).

## Connection to E_6

The Schlafli graph encodes the incidence geometry of the 27 lines on a cubic surface, which is governed by the E_6 root system:

- |Aut(Schlafli)| = 51840 = |W(E_6)|
- The 27 lines correspond to the 27 weights of the fundamental representation of E_6
- The Hecke algebra H(PSp(4,3), (C_2^4):A_5) is the Bose-Mesner algebra of the E_6 association scheme restricted to the 27-weight orbit

The chain Petersen -> T(6) -> Schlafli corresponds to the restriction from E_6 to A_4 (= sl(5)) symmetry via the Clebsch surface. The S_5 = W(A_4) acting on the 27 lines is the Weyl group of the A_4 subalgebra of E_6.

## Computation

All results computed in GAP 4 with the GRAPE package. Script: `hecke_clebsch.g` in this directory.

## Open questions

1. **The 12-vertex graph.** Is the induced subgraph on the 12 Type II lines a known graph? Its automorphism group C_2 x S_6 (order 1440) is unexpectedly large. Identify it.

2. **The Hecke algebra of T(6).** Compute H(S_6, S_2 x S_4) or whichever pair gives T(6) as a Sabidussi graph. What is the rank? T(6) has 4 distinct eigenvalues (6, 1, 1, -3 ... actually 3: 6, 1, -3), so it should also be rank 3.

3. **The restriction map.** Make the embedding of Hecke algebras explicit: the map H(Schlafli) -> H(T(6)) -> H(Petersen) induced by restricting the 27 vertices to 15 to 10. Is this an algebra homomorphism?

4. **The Iwahori-Hecke connection.** Compute H(GL(2, F_q), B) for small q. The quadratic relation T_s^2 = (q-1)T_s + qT_0 should appear. Compare with the SRG structure constants.

5. **Non-Gelfand pairs.** Find a pair (G, H) from the graphviz census where the permutation representation has an irrep of multiplicity > 1 (i.e., the Hecke algebra is non-commutative). This would be the first example where the centralizer algebra has genuinely new structure.
