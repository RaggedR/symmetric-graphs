# Mathematics of the Bolza Surface Uniformization

## The Problem

We have TWO descriptions of the same object:

1. **Hyperbolic picture**: The Bolza surface is H/Gamma, where H is the upper half-plane and Gamma is a torsion-free Fuchsian group. We have the Poincare disk tiling (48 triangles forming a regular hyperbolic octagon) and the group GL(2,3) acting as symmetries.

2. **Algebraic picture**: The Bolza surface is the curve y^2 = x^5 - x = x(x-1)(x+1)(x-i)(x+i) in C^2. This has 6 branch points (the roots plus infinity) and genus 2 by Riemann-Hurwitz.

We want an EXPLICIT MAP between them: given a point z in the Poincare disk (inside the fundamental octagon), compute the corresponding point (x(z), y(z)) on the algebraic curve. This map is the **uniformization** — one of the central objects in Riemann surface theory.

## The Groups

### The (2,3,8) triangle group Delta

Delta = <a, b | a^2 = b^3 = (ab)^8 = 1> acts on H via SL(2,R):

    a_real = [[0, 1], [-1, 0]]       trace = 0,   order 2 in PSL(2,R)
    b_real = [[0.5, -f*sin(pi/3)],    trace = 1,   order 3 in PSL(2,R)
              [sin(pi/3)/f, 0.5]]

where f = (rhs + sqrt(rhs^2 - 4))/2, rhs = 2*cos(pi/8)/sin(pi/3).

c_real = a_real @ b_real has trace 2*cos(pi/8), order 8 in PSL(2,R).

The three fixed points z_a, z_b, z_c form the vertices of a (pi/2, pi/3, pi/8) hyperbolic triangle that tiles H.

### GL(2,3) — the automorphism group

GL(2, F_3) has order 48. It consists of all invertible 2x2 matrices over the field with 3 elements. There is a surjection Delta -> GL(2,3) sending a -> a_gl, b -> b_gl where a_gl has order 2 and b_gl has order 3 in GL(2,3).

### Gamma — the surface group

Gamma = ker(Delta -> GL(2,3)) is a torsion-free normal subgroup of Delta with index 48. The quotient H/Gamma is the Bolza surface. The quotient GL(2,3) = Delta/Gamma acts as the automorphism group of the surface.

**Key property**: Gamma is isomorphic to the fundamental group of a genus-2 surface:
    
    Gamma = pi_1(Bolza) = <a1, b1, a2, b2 | [a1,b1][a2,b2] = 1>

The generators a1, b1, a2, b2 are the **side-pairing transformations** of the octagon fundamental domain. They are the primary objects we need to compute.

### Relationship between the groups

    1 -> Gamma -> Delta -> GL(2,3) -> 1     (exact sequence)

Every element of Delta has:
- An SL(2,R) matrix (the Mobius transformation on H)
- A GL(2,3) element (its image under the surjection)

Two Delta elements with the SAME GL(2,3) image differ by a Gamma element:
    
    delta_1 * delta_2^{-1} in Gamma  iff  they map to the same GL(2,3) element

This is the master principle: **two triangles with the same GL(2,3) label are Gamma-translates of each other.**

## The Fundamental Domain

The 48 (2,3,8) triangles closest to the center (in hyperbolic distance) form a Dirichlet fundamental domain for Gamma. This domain is (approximately) a regular hyperbolic octagon.

Each triangle carries a GL(2,3) label (one of 48 elements) and an SL(2,R) matrix. The domain has exactly one triangle per label.

### Important subtlety: SL(2,R) vs PSL(2,R)

In SL(2,R), the matrix -I (negative identity) acts trivially on H: it gives the same Mobius transformation as I. So M and -M represent the same hyperbolic isometry but different SL(2,R) elements. In GL(2,3), -I maps to [[2,0],[0,2]] mod 3, which is the central involution (the hyperelliptic involution).

The BFS through the triangle group generates both M and -M as separate elements. They produce the same triangle geometry. When computing side pairings, the ratio of two matrices that are -I translates gives gamma = -I (trace = -2), which is TRIVIAL in PSL(2,R). These must be filtered out.

## Side-Pairing Transformations

### The algebraic method (Attempt 10 breakthrough)

For a fund domain triangle with GL(2,3) label g and SL(2,R) matrix M, its neighbor via generator 'a' in the tiling has:
- GL(2,3) label: g' = g * a_gl
- Expected SL(2,R) matrix (if in the same Gamma copy): M_expected = M @ a_real
- Actual fund domain matrix (the representative for label g'): M_actual = M_{fund(g')}

The side-pairing element is:

    gamma = M_actual @ M_expected^{-1} = M_{fund(g')} @ (M @ a_real)^{-1}

If gamma = +/- I, the neighbor is in the same Gamma copy (interior edge).
If gamma has |tr| > 2, it's a HYPERBOLIC element (genuine side-pairing generator).

### Results

14 distinct Gamma elements found, ALL with |trace| = 4.8284 = 2 + 2*sqrt(2).

This trace corresponds to the systole of the Bolza surface: the shortest closed geodesic has length 2*arccosh(1 + sqrt(2)) = 2*arccosh(2.414) = 3.057.

The 14 elements include the 4 generators and their inverses, plus some conjugates. (The Bolza surface has 12 systole geodesics due to its high symmetry; each appears as a side-pairing element when crossed by the octagon boundary.)

## The Hyperelliptic Involution and Weierstrass Points

### The involution

The hyperelliptic involution iota is the unique central involution of GL(2,3):

    iota = -I = [[2,0],[0,2]] mod 3

In the (2,3,8) presentation: iota = c^4 = (ab)^4.

In SL(2,R): c_real^4 is a rotation by pi about z_c (the order-8 fixed point). In the recentered Poincare disk, this acts as z -> -z.

On the algebraic curve y^2 = x^5 - x, the involution is (x, y) -> (x, -y) (sheet swap).

### Weierstrass points

The 6 Weierstrass points are the fixed points of iota on the surface. They are the branch points of the hyperelliptic projection x: Bolza -> CP^1.

On the algebraic curve, they are at (x, 0) where x^5 - x = 0, i.e., x in {0, 1, -1, i, -i}, plus the point at infinity.

On the universal cover H, a point z is a Weierstrass point iff:

    exists gamma in Gamma such that gamma(iota_tilde(z)) = z

where iota_tilde = c_real^4 (the SL(2,R) lift of iota, acting as z -> -z on the recentered disk).

Equivalently: z is a fixed point of the composition gamma * c_real^4 for some gamma in Gamma.

For gamma = I: fixed point is z = 0 (the center). This gives ONE Weierstrass point.

For other gamma: the fixed point of gamma * c_real^4 is an order-2 elliptic fixed point in H. Each such point is a Weierstrass point. The 6 Weierstrass points correspond to 6 Gamma-orbits.

### UNSOLVED: correct deduplication

The computation produces 17 candidate points (for 14 Gamma elements plus identity). These must be reduced to 6 by identifying Gamma-equivalent points. Two candidate points z1, z2 are the same Weierstrass point iff there exists gamma in Gamma with gamma(z1) = z2.

In the fundamental domain:
- z = 0 (center): 1 Weierstrass point
- z = octagon vertices (all at |w| = 0.6436): these are ALL identified to a single point P on the surface (all 8 vertices merge to one point). So {V0, ..., V7} is the 2nd Weierstrass point.
- Points at |w| = 0.91 (outside the octagon): these are Gamma-translates of interior points. Need to map them back into the fundamental domain to identify which surface point they represent.

**Question**: are the 8 octagon vertices really one Weierstrass point, or do they split into multiple surface points? All 8 are identified by the side pairings, so they are ONE point P. And c^4 fixes P (since c^4 permutes the vertices cyclically by 4 positions, and they're all identified). So yes, P is a Weierstrass point.

That gives 2 Weierstrass points so far (center and P). The remaining 4 must come from the |w| = 0.91 candidates, properly deduplicated. With 8 candidates at |w| = 0.91, and the octahedral symmetry, these should form 1 orbit of size 8 under GL(2,3), giving 8/2 = 4 Weierstrass points (each with a Gamma-translate pair)... but this needs verification.

**Conjecture**: The 6 Weierstrass points in the fundamental domain are:
1. The center (z = 0)
2. The octagon vertex identification point P
3-6. Four points at |w| = 0.91, identified in pairs under the octahedral symmetry

## The Holomorphic Differentials

### Theory

On a genus-2 Riemann surface, the space of holomorphic 1-forms (abelian differentials of the first kind) has dimension g = 2. On the algebraic curve y^2 = x^5 - x, these are:

    omega_1 = dx/(2y)
    omega_2 = x * dx/(2y)

Their ratio is: omega_2 / omega_1 = x (the hyperelliptic projection).

### As automorphic forms

Pulling back to H via the uniformization, omega_k = phi_k(z) dz where phi_k is a weight-2 automorphic form for Gamma:

    phi_k(gamma(z)) = (c*z + d)^2 * phi_k(z)   for gamma = [[a,b],[c,d]] in Gamma

The FUNCTION x(z) = phi_2(z) / phi_1(z) is weight-0 (Gamma-invariant):

    x(gamma(z)) = x(z)   for all gamma in Gamma

This is the key: x is a FUNCTION on the surface, with no automorphic twist.

### Computing x: three approaches attempted

**Approach 1 (Attempt 8): Poincare series**

    Theta_p(z) = sum_{gamma in Gamma} 1/((a_g - p*c_g)*z + (b_g - p*d_g))^2

Two series with different poles p1, p2 give weight-2 forms whose ratio is x (up to Mobius).

Result: 5% automorphy error with 1101 Gamma elements. The error is position-dependent: at gamma_0(z), the series needs Gamma elements near gamma_0^{-1}(pole neighborhood), which aren't in the truncation.

**Approach 2 (Attempt 9): Discrete Hodge theory (twisted BCs)**

Solve for phi as a weight-2 form: Delta(phi) = 0 with phi(gamma(z)) = j(gamma,z)^2 * phi(z) on the boundary. The kernel of the twisted Laplacian gives the holomorphic differentials.

Result: Failed due to incorrect boundary identification (geometric vertex matching doesn't work; the octagon boundary doesn't align with triangle edges).

**Approach 3 (Attempt 10): Discrete Laplacian with algebraic boundary identification**

Since x is Gamma-INVARIANT (not twisted), solve Delta(x) = 0 with x(z) = x(z') for identified boundary vertices. The identification is computed ALGEBRAICALLY from GL(2,3) labels (the Attempt 10 breakthrough).

Result: 14 side-pairing generators found, 80 boundary pairs matched, Laplacian solve succeeds. Eigenvalue spectrum shows lambda_0 = 0 (constant kernel), then pairs at ~1e-3.

### The eigenvalue spectrum and what it means

The reduced Laplacian (after merging identified boundary vertices) has spectrum:

    lambda_0 = 0            (the constant function)
    lambda_1,2 = ~1e-3      (conjugate pair)
    lambda_3,4 = ~1e-3      (conjugate pair)
    lambda_5   = ~1e-3

On a properly closed genus-2 surface, the discrete Laplacian should have:
- A 1-dimensional kernel (constants)
- First nonzero eigenvalue proportional to 1/area (the spectral gap)
- For the Bolza surface (maximal symmetry), the first nonzero eigenvalue has multiplicity 5 (the 5-dimensional representation of GL(2,3))

The fact that we see lambda_1 through lambda_5 all at ~1e-3 is CONSISTENT with the Bolza surface spectrum. The nonzero gap means the surface IS "closed" (all boundary identifications are working).

However, the Laplacian eigenvalues are NOT the holomorphic differentials. The Laplacian kernel on a closed surface is only the constants. The holomorphic differentials live in the kernel of the d-bar operator, not the Laplacian.

### CRITICAL INSIGHT: Laplacian kernel != holomorphic differentials

On a compact Riemann surface:
- The Laplacian Delta = d*d (Hodge Laplacian on functions) has kernel = {constants}
- The Hodge Laplacian on 1-FORMS has kernel = harmonic 1-forms, dimension 2g = 4
- The holomorphic differentials are the (1,0)-part of the harmonic forms, dimension g = 2

So the right approach is NOT to solve Delta(f) = 0 on functions. Instead:
- Either solve for the Hodge Laplacian on 1-forms (requires a 1-form discretization)
- Or solve for x directly as a MEROMORPHIC function (has a pole, so Delta(x) != 0 everywhere)
- Or use the Poincare series (Attempt 8) with better convergence

### The meromorphic function approach

The function x is harmonic everywhere EXCEPT at its pole (one Weierstrass point, x = infinity). On a genus-2 surface, a meromorphic function with a double pole at P is determined (up to scaling and additive constant) by:
- Delta(x) = 0 away from P
- x ~ 1/(z - z_P)^2 near P (in local coordinates)
- x(z) = x(z') for identified boundary vertices

This is a Laplacian with a POINT SOURCE at P (a Green's function problem):

    Delta(x) = c * delta_P   (Dirac delta at the pole)

where c is determined by the normalization. This is solvable as a sparse linear system with one modified equation (at the pole vertex).

Alternatively: remove the pole vertex from the mesh, solve Delta(x) = 0 on the remaining mesh with:
- Fixed values at the 5 finite Weierstrass points: x(W_k) = known values
- Free boundary at the removed vertex (Neumann condition)
- Identified boundary pairs: x(z) = x(z')

But this solves for the REAL PART of x only. The imaginary part requires the Cauchy-Riemann equations, which couple the two components.

### Resolution: solve for x as a complex harmonic function

On the Poincare disk (conformally equivalent to C), x is holomorphic = Re(x) and Im(x) are both harmonic. But they're COUPLED by the Cauchy-Riemann equations. On a simply connected domain (the octagon before identification), knowing Re(x) at a boundary determines Im(x) up to a constant (the harmonic conjugate). But on the SURFACE (genus 2), the harmonic conjugate has PERIODS around the handles.

The correct formulation: solve for x in C as a complex-valued function satisfying:
1. The DISCRETE HOLOMORPHICITY condition (not just harmonicity) at interior vertices
2. x(z) = x(z') at identified boundary pairs (no twist — x is a function)
3. Fixed values at 5 Weierstrass points
4. One vertex removed or unconstrained (the pole)

The discrete holomorphicity condition on a triangulated surface requires the d-bar operator, not the Laplacian. For a triangulation of a domain in C, the discrete d-bar condition at a vertex v is:

    sum_{triangles T around v} (area(T))^{-1} * sum_{edges e of T} x(head(e)) * (tail(e)_perp - head(e)_perp) = 0

This is more involved than the cotangent Laplacian and depends on the specific triangulation.

## What the Next Claude Needs to Do

### Option A: Fix the Poincare series (Attempt 8 improvement)

The series works but converges slowly. Two fixes:
1. **Adaptive summation**: For each evaluation point z, include Gamma elements whose translates bring z near the pole. This requires a spatial index on Gamma-translates.
2. **Use more elements**: Increase BFS from 50K to 500K Delta elements (~10K Gamma elements). The 5% error should drop to ~0.5%.

The Mobius calibration also needs improvement. Try all 48 assignments (2 * 4!) more carefully, using the automorphy check as the quality criterion (not just the Weierstrass point fit).

### Option B: Discrete d-bar operator (Attempt 10 improvement)

Replace the scalar Laplacian with the discrete d-bar operator. This directly encodes holomorphicity, not just harmonicity. The kernel of the discrete d-bar (with boundary identification) should give the holomorphic differentials.

References: Mercat (2001) "Discrete Riemann surfaces", Bobenko-Suris (2008) "Discrete differential geometry".

### Option C: Period matrix route

1. Compute the period matrix Omega of y^2 = x^5 - x numerically (integrate dx/(2y) and x*dx/(2y) along homology cycles — standard numerical quadrature with Gauss-Jacobi for the sqrt singularities at branch points)
2. The Riemann theta function Theta(z; Omega) is computable (fast convergence for genus 2)
3. The inverse Abel-Jacobi map (Jacobian -> curve) is given by derivatives of theta
4. Matching the Fuchsian group to the Jacobian: the Abel-Jacobi image of the surface group generators determines the lattice in C^2

This is the most well-established numerical approach (Deconinck-van Hoeij algorithm, implemented in Maple/SageMath).

### Option D: Schottky uniformization

Represent Gamma as a Schottky group (free group on 2 generators acting on CP^1 by Mobius transformations pairing circles). The Schottky-Klein prime function converges EXPONENTIALLY fast. The automorphic functions are products over the Schottky group, not sums over a Fuchsian group.

Converting between Fuchsian and Schottky representations requires the retrosection problem (finding the Schottky generators from the period matrix).

## Available Code and Data

All in `/Users/robin/git/graphviz/bolza-surface/`:

| File | What it provides |
|------|------------------|
| `bolza_group.py` | GL(2,3) class, TriangleGroup class, BolzaTiling class |
| `bolza_r4.py` | Attempt 8: Poincare series, 1101 Gamma elements, Mobius calibration |
| `bolza_r4_v2.py` | Attempt 9: Discrete Hodge (failed on boundary identification) |
| `bolza_r4_v3.py` | Attempt 10: Algebraic side pairings (14 generators found), Laplacian solve |
| `build_viz.py` | Original tiling visualization (Poincare disk + 3D torus) |
| `bolza_cohomology.py` | Cell complex, boundary maps, H^1 decomposition |
| `octagon-folding-attempts.md` | Detailed log of all 10 attempts |

### Key computed data (reusable)

- **SL(2,R) generators**: a = [[0,1],[-1,0]], b = [[0.5,-1.246],[0.602,0.5]], c = a@b
- **Fixed points**: z_a = i, z_b = 1.438i, z_c = 0.644 + 0.765i
- **Disk recentering**: w_c_raw = to_disk(z_c), recenter(w) = (w - w_c_raw)/(1 - conj(w_c_raw)*w)
- **14 Gamma generators**: all with |tr| = 4.828, det = 1 (in bolza_r4_v3.py output)
- **Systole**: length 2*arccosh(1+sqrt(2)) = 3.057, trace 2+2*sqrt(2) = 4.828
- **48 fund domain triangles** with GL(2,3) labels and SL(2,R) matrices

### Key mathematical facts

- The Bolza surface is the UNIQUE genus-2 surface with 48 automorphisms
- Its automorphism group GL(2,3) has the (2,3,8) presentation
- The curve equation y^2 = x^5 - x has branch points at {0, 1, -1, i, -i, infinity}
- The hyperelliptic involution is c^4 (rotation by pi in the disk, sheet swap on the curve)
- The space of holomorphic differentials is 2-dimensional (genus = 2)
- x = omega_2/omega_1 = ratio of holomorphic differentials = the hyperelliptic projection
- x has degree 2, with 6 branch points (the Weierstrass points)
