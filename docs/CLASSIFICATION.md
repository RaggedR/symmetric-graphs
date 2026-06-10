# The Groups: A Classification Guide

## Finite Field Conjugation

Over **R**, complex conjugation z -> z̄ is the unique nontrivial automorphism of C fixing R. It's the generator of Gal(C/R) = Z/2. The unitary group SU(2) is defined using this automorphism: M†M = I, where † = conjugate transpose.

Over **F\_q**, the same structure exists. The degree-2 extension F\_{q²}/F\_q has Galois group Z/2 with generator the **Frobenius automorphism**:

```
σ: x ↦ x^q
```

This is "finite field conjugation." It satisfies all the same properties:

| Property | C/R | F_{q²}/F_q |
|----------|-----|------------|
| Generator of Gal | z ↦ z̄ | x ↦ x^q |
| Fixes base field | z̄ = z iff z ∈ R | x^q = x iff x ∈ F_q |
| Order 2 | z̄̄ = z | x^{q²} = x |
| Field automorphism | ✓ | ✓ |

So one defines:

```
SU(2, q) = { M ∈ SL(2, F_{q²}) : M^{(q)T} M = I }
```

where M^{(q)} means applying x ↦ x^q to every entry.

**But it doesn't give anything new.** Both |SL(2,q)| and |SU(2,q)| equal q(q²-1). At the projective level, PSU(2,q) ≅ PSL(2,q) — they are the same simple group. The split form (SL) and non-split form (SU) of type A\_1 collapse over finite fields.

Over R, the distinction is profound: SU(2) is compact (the 3-sphere S³), SL(2,R) is non-compact. This drives all of physics — compact = bound states, non-compact = scattering. Over F\_q, compactness is meaningless (everything is finite), and the distinction evaporates.


## The GL(2,q) Family

Four groups, related by two independent operations: **restrict to det = 1** and **quotient by scalar matrices**.

```
         GL(2,q)
        /       \
     det=1     mod F_q* · I
      /             \
   SL(2,q)       PGL(2,q)
      \             /
    mod {±I}   restrict
        \       /
        PSL(2,q)
```

### Orders

| Group | Definition | Order (odd q) |
|-------|-----------|---------------|
| GL(2,q) | All invertible 2×2 over F_q | q(q-1)²(q+1) |
| SL(2,q) | ker(det: GL → F_q\*) | q(q²-1) |
| PGL(2,q) | GL(2,q) / {aI : a ∈ F_q\*} | q(q²-1) |
| PSL(2,q) | SL(2,q) / {±I} | q(q²-1)/2 |

SL and PGL have the **same order** but are different groups. SL has a nontrivial centre ({±I}); PGL has trivial centre. They are both "halfway" between GL and PSL, reached by different routes.

### Exact sequences

**Determinant sequence:**
```
1 → SL(2,q) → GL(2,q) --det--> F_q* → 1
```

**Spinor sequence:**
```
1 → {±I} → SL(2,q) → PSL(2,q) → 1
```

For odd q, this is a Z/2 quotient: |PSL| = |SL|/2. For q = 2^n (characteristic 2), -I = I, so PSL = SL.

**Square class sequence:**
```
1 → PSL(2,q) → PGL(2,q) → F_q*/(F_q*)² → 1
```

For odd q, exactly half the nonzero field elements are squares, so F\_q\*/(F\_q\*)² ≅ Z/2 and [PGL : PSL] = 2. An element of PGL(2,q) lies in PSL(2,q) iff its determinant is a square in F\_q\*.

This is the finite field analogue of orientation: "even" vs "odd" transformations are detected by quadratic residuosity, not by sign.

### For prime q (the relevant case)

When q = p is prime, a simplification: Gal(F\_q/F\_p) is trivial, so

```
Aut(PSL(2,p)) = PGL(2,p)
```

No field automorphisms. PGL is the **full** automorphism group of the simple group PSL. For prime powers q = p^f with f > 1, there is an additional factor:

```
Aut(PSL(2,q)) = PΓL(2,q) = PGL(2,q) ⋊ Gal(F_q/F_p) = PGL(2,q) ⋊ Z/f
```

All primes in our census (q = 5, 7, 13, 17) have f = 1, so Aut = PGL.


## The Adjoint Representation: PSL(2,q) ≅ SO(3,q)

### The Lie algebra sl(2)

sl(2, F\_q) is the 3-dimensional vector space of traceless 2×2 matrices over F\_q. Standard basis:

```
e = (0 1)    f = (0 0)    h = (1  0)
    (0 0)        (1 0)        (0 -1)
```

with bracket [X, Y] = XY - YX satisfying:

```
[h, e] = 2e      [h, f] = -2f      [e, f] = h
```

### The trace form

The inner product on sl(2) is:

```
⟨X, Y⟩ = tr(XY)
```

Gram matrix in the basis {e, f, h}:

```
       e  f  h
  e  ( 0  1  0 )
  f  ( 1  0  0 )
  h  ( 0  0  2 )
```

with det = -2 ≠ 0 (non-degenerate for any odd q).

The Killing form B(X,Y) = tr(ad(X) ad(Y)) equals exactly 4 · tr(XY) for sl(2). Both define the same orthogonal group.

**Structure of the form.** The basis {e, f, h} is not orthogonal. The elements e and f are **null vectors** (⟨e,e⟩ = ⟨f,f⟩ = 0) paired with each other (⟨e,f⟩ = 1), while h is orthogonal to both with ⟨h,h⟩ = 2. This is a **hyperbolic plane plus a line** — the standard non-degenerate form in 3 dimensions. Over F\_q in odd dimension, there is only one such form up to similarity, so this is *the* SO(3,q).

The null vectors e, f are the upper and lower triangular nilpotent elements — the raising and lowering operators from quantum mechanics. Their nullity reflects the root space decomposition sl(2) = n⁺ ⊕ t ⊕ n⁻.

### The adjoint action

SL(2,q) acts on sl(2) by conjugation:

```
Ad(g)(X) = gXg⁻¹
```

This preserves the trace form:

```
⟨gXg⁻¹, gYg⁻¹⟩ = tr(gXg⁻¹ · gYg⁻¹) = tr(gXYg⁻¹) = tr(XY) = ⟨X, Y⟩
```

So Ad(g) ∈ O(3, q) for every g ∈ SL(2,q). Since Ad is a morphism of algebraic groups and Ad(I) = I has determinant 1, in fact Ad(g) ∈ SO(3, q).

### The explicit 3×3 matrices

For g = (a b; c d) ∈ SL(2,q) with ad - bc = 1, the matrix of Ad(g) in the basis {e, f, h} is:

```
             ( a²    -b²    2ab  )
Ad(g)   =   ( -c²    d²   -2cd  )
             ( -ac    bd   ad+bc )
```

One can verify Ad(g)^T G Ad(g) = G where G is the Gram matrix above.

### The isomorphism

The kernel of Ad: SL(2,q) → SO(3,q) is exactly the centre {±I} (the elements that commute with all of sl(2)). This gives:

```
PSL(2,q) = SL(2,q) / {±I}  ≅  SO(3,q)
```

For odd q in odd dimension, Ω(3,q) = SO(3,q) (the commutator subgroup equals the full special orthogonal group), so this is a clean isomorphism with no subtlety about Ω vs SO.


## SO(2,q): Finite Field Circles

Over R, SO(2) is the circle group — rotations of the plane. Over F\_q, there is no circle, but there are **two** finite analogues, because there are two non-equivalent non-degenerate symmetric bilinear forms on F\_q² (unlike odd dimension, where there is only one).

### Split form (hyperbolic)

Form matrix: B = (0 1; 1 0), i.e. B(x,y) = x₁y₂ + x₂y₁.

```
SO⁺(2,q) = { (a   0 ) : a ∈ F_q* } ≅ F_q*       (cyclic, order q-1)
             (0  a⁻¹)
```

Diagonal matrices scaling the two null directions inversely. The form has two isotropic lines, and SO⁺ just permutes them — these are not genuine "rotations."

### Non-split form (anisotropic)

Form matrix: B = (1 0; 0 -ε) where ε is a fixed non-square in F\_q\*.

```
SO⁻(2,q) ≅ ker(N: F_{q²}* → F_q*)              (cyclic, order q+1)
```

where N is the field norm x ↦ x^{q+1}. No null vectors exist — every nonzero vector has nonzero length. These are genuine "rotations."

Both are **cyclic** (abelian). Over R, SO(2) ≅ S¹ is also abelian but continuous; over F\_q the circle collapses to a finite cyclic group.

### SO(2,q) inside SO(3,q)

Over R, SO(2) embeds in SO(3) as rotations around an axis: pick a unit vector v, and rotate in the plane v⊥. The same works over F\_q.

In SO(3,q) acting on (sl(2), tr(XY)), pick a non-isotropic vector and look at its stabilizer. The two natural choices give the two types of torus:

**The h-axis.** The vector h has ⟨h,h⟩ = 2 (a nonzero square for q ≡ 1 mod 4, non-square for q ≡ 3 mod 4 — but always nonzero). Its orthogonal complement h⊥ = span{e, f} carries the restricted form:

```
       e  f
  e  ( 0  1 )
  f  ( 1  0 )
```

This is the split (hyperbolic) form. The stabilizer of the h-axis is SO⁺(2,q) ≅ F\_q\*, order q−1. In PSL(2,q), this is the image of the diagonal subgroup — matrices (a 0; 0 a⁻¹) modulo ±I, giving a cyclic group of order (q−1)/2.

**An anisotropic axis.** Pick w ∈ sl(2) with ⟨w,w⟩ a non-square. The complement w⊥ carries the non-split form, and the stabilizer is SO⁻(2,q), order q+1. In PSL(2,q), this is a cyclic subgroup of order (q+1)/2, arising from the multiplicative group of F\_{q²} via the embedding F\_{q²}\* → GL(2,q).

These are the **two conjugacy classes of maximal tori** in PSL(2,q) ≅ SO(3,q):

| Torus | Form on axis⊥ | Order in SO(3,q) | In PSL(2,q) | Source |
|-------|---------------|-------------------|-------------|--------|
| Split | hyperbolic | q−1 | (q−1)/2 | diagonal matrices mod ±I |
| Non-split | anisotropic | q+1 | (q+1)/2 | F\_{q²}\* mod ±I |

The factor of 2 between SO(2) orders and PSL torus orders comes from the isomorphism PSL ≅ SO(3) passing through the adjoint representation, which sends ±g ∈ SL(2) to the same element of SO(3).

**A new phenomenon.** Over R, there is one SO(2) up to conjugacy — all axes are equivalent. Over F\_q, there are two, because the orthogonal complement can carry either the split or non-split 2D form. There are two "shapes of circle" over a finite field.

**Connection to the graph drawings.** The orbits of these tori on the vertex set of a Sabidussi graph give decompositions into concentric rings. The split torus gives (q−1)/2 orbits, the non-split gives (q+1)/2. This connects to the orbit quotient structures that the drawing algorithm selects.


## Which Groups Act on the Graphs

For the three families of cubic arc-transitive graphs from PSL(2,q):

| Family | H | Tutte level | Construction group | Full Aut(Γ) |
|--------|---|-------------|-------------------|-------------|
| 1 | C₃ | 1-arc-regular | PSL(2,q) | PSL(2,q) × C₂ or similar |
| 2 | S₃ | 2-arc-transitive | PSL(2,q) | **PGL(2,q)** |
| 3 | S₄ | 4-arc-transitive | PSL(2,q) | **PGL(2,q)** |

The Sabidussi construction uses PSL(2,q) as the group G, but the **full automorphism group** of the resulting graph can be larger. For the S₃ and S₄ families, the extra Z/2 from PGL/PSL ≅ F\_q\*/(F\_q\*)² gives additional symmetry — it's this extra symmetry that promotes 1-arc-regularity to 2-arc-transitivity (S₃ family) or to 4-arc-transitivity (S₄ family).

This is why the Foster census records:
- Petersen: Aut = PGL(2,5) = S₅ (S₃ family, q=5)
- Heawood: Aut = PGL(2,7) (S₄ family, q=7)
- Desargues: Aut = PGL(2,5) × C₂ (S₃ family, q=5, bipartite doubling)
- Biggs-Smith: Aut = PGL(2,17) (S₄ family, q=17)
- F182A: Aut = PGL(2,13) (S₃ family, q=13)

In each case, PSL constructs the graph, but PGL is its full symmetry group.


## Summary: the diamond for prime q

```
GL(2,q)                          — all invertible 2×2 matrices
   |  \
   |   SL(2,q)                   — det = 1, centre = {±I}
   |      |
   |      | mod {±I}             — spinor quotient
   |      |
PGL(2,q) = Aut(PSL(2,q))        — full graph automorphism group
   |
   | index 2                     — square classes F_q*/(F_q*)²
   |
PSL(2,q) ≅ SO(3,q)              — simple group, construction group
```

For prime q, this is the complete picture. No field automorphisms, no forms to distinguish, no complications. PGL(2,q) is as big as it gets.
