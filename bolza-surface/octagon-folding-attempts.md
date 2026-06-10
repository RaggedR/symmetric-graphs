# Folding the Octagon into a Double Torus

## The problem

We have 48 triangles in the Poincaré disk forming the Bolza surface's Dirichlet fundamental domain (a regular hyperbolic octagon centered at the order-8 vertex). We want to map each triangle's disk coordinates to 3D positions on a double torus.

## The octagon identification

Sides labeled counterclockwise: a₁ b₁ a₁⁻¹ b₁⁻¹ a₂ b₂ a₂⁻¹ b₂⁻¹

Identification pairs (with reversed orientation):
- a₁ (side 1, V₀→V₁) ↔ a₁⁻¹ (side 3, V₃→V₂)
- b₁ (side 2, V₁→V₂) ↔ b₁⁻¹ (side 4, V₄→V₃)
- a₂ (side 5, V₄→V₅) ↔ a₂⁻¹ (side 7, V₇→V₆)
- b₂ (side 6, V₅→V₆) ↔ b₂⁻¹ (side 8, V₀→V₇)  ← note: wraps back to V₀

All 8 vertices merge to ONE point P on the surface.
χ = 1 vertex − 4 edges + 1 face = −2, genus 2. ✓

## Approach: Cut into two tori

**Cut:** Geodesic from V₀ to V₄ through the octagon center. Since V₀ = V₄ = P on the surface, this is a **closed curve** — it's a simple closed loop separating the two handles.

This gives two quadrilaterals:
- **Quad 1** (sides a₁, b₁, a₁⁻¹, b₁⁻¹, cut): after identification → torus 1 with hole
- **Quad 2** (sides a₂, b₂, a₂⁻¹, b₂⁻¹, cut): after identification → torus 2 with hole

The cut edges become the **neck** connecting the two tori.

## Coordinate system

### Within each quadrilateral: (s, t)

Each quad has 4 identified sides forming a torus, plus the cut side. Ignoring the cut for a moment, the 4 sides define a torus with coordinates:
- **s ∈ [0, 1]**: the "a direction" — going from side aₖ to side aₖ⁻¹ (around the big circle)
- **t ∈ [0, 1]**: the "b direction" — going from side bₖ to side bₖ⁻¹ (around the tube)

### 3D embedding: (s, t, handle) → (x, y, z)

Each handle is a standard torus:

```
Torus 1: center at (−d, 0, 0)
  x = −d + (R + r·cos(2πt)) · cos(2πs)
  y = (R + r·cos(2πt)) · sin(2πs)
  z = r · sin(2πt)

Torus 2: center at (+d, 0, 0), same formula with +d
```

The neck connects them where the cut was: at a specific (s₀, t₀) position on each torus, a hole is cut and a smooth tube bridges them.

## How to compute (s, t) from disk coordinates

This is the hard part. Each point in the octagon has Poincaré disk coordinates (x, y). I need to map to (s, t) within its quadrilateral.

### Attempt 1: Bilinear interpolation (Euclidean)

Find the 4 corners of the quadrilateral in disk coords. Use bilinear interpolation to compute (s, t). Simple but ignores hyperbolic geometry — will distort triangles near the boundary.

### Attempt 2: Hyperbolic barycentric coordinates

Use the hyperbolic metric to compute "how far along" each side the point is. More correct but more complex.

### Attempt 3: Use the (ρ, φ) polar coordinates directly

Each point has (ρ, φ) centered at the order-8 vertex.
- φ determines **which sector** (and hence handle + s position)
- ρ determines **t position** (tube angle)

The 8 sectors (one per vertex of the octagon) map as:
- Sectors 0,1,2,3 → handle 1, s = φ/(π)
- Sectors 4,5,6,7 → handle 2, s = (φ−π)/π

And ρ → t: center (ρ=0) = inner tube (t=0), boundary (ρ=max) = outer tube (t=0.5)?

**Problem with this:** ρ=0 is a single point (the center), but t=0 should be a circle (the tube inner edge). So ρ can't map linearly to t — we need something like t = f(ρ) where f expands the center point into a circle.

Actually: at the center of the octagon, we're at the **neck** between the two tori. The neck is a circle, not a point. So the mapping must "blow up" the center into a circle. This is topologically a surgery — replacing a point with a circle.

In practice: map ρ to v (tube angle) with v ranging from -π to π as ρ goes from the cut side to the opposite boundary. Points near the cut (ρ near 0 from one direction, or equivalently ρ near 0 from the other direction across the cut) map to the neck circle.

## What could go wrong

1. **Triangles straddling the cut**: some of the 48 triangles may have vertices on both sides of the V₀-V₄ geodesic. These need to be split or assigned to one side.

2. **The center singularity**: the (ρ, φ) → (u, v) map has a singularity at ρ=0 (one point maps to a circle). Need smooth regularization.

3. **Side identification mismatch**: when side aₖ is glued to aₖ⁻¹ with reversed orientation, the s coordinate must reverse. If the bilinear interpolation doesn't respect this, the resulting mesh will have twisted triangles.

4. **The Dirichlet domain might not be exactly the standard octagon**: the 48-triangle Dirichlet domain from the distance-based selection might differ from the canonical fundamental domain with side-identification a₁b₁a₁⁻¹b₁⁻¹a₂b₂a₂⁻¹b₂⁻¹. The identification pattern depends on WHICH generators we choose, not just the geometry.

5. **Non-uniform triangle sizes in 3D**: hyperbolic triangles are all the same hyperbolic size, but in 3D they'll have varying Euclidean size (large on the outer torus, small at the neck). This is unavoidable.

6. **Which vertex is V₀?**: the octagon has 8 vertices, and I need to identify which one is V₀ in my disk tiling to determine the cut direction. The 8 vertices are all at the same distance from center (they're vertices of the regular octagon). I need to pick one and compute the cut geodesic from it to the opposite vertex V₄.

## Plan of attack

1. **Identify the octagon vertices** in the disk tiling (the 8 order-8 points closest to the boundary of the fundamental domain)
2. **Choose V₀ and V₄** (opposite vertices, defining the cut)
3. **Classify each of the 48 triangles** as handle 1 or handle 2 (which side of the cut?)
4. **Compute (s, t) for each triangle vertex** using polar coordinates (φ → s, ρ → t)
5. **Map to 3D torus coordinates** with neck blending
6. **Render**

## First attempt strategy

Start with the simplest version:
- Use angle φ to assign handle and compute s
- Use distance ρ to compute t (with some nonlinear mapping to avoid the center singularity)
- Don't worry about the neck — just put the two tori side by side and see if the triangles tile them correctly
- Fix the neck in a second pass

The first attempt will probably produce a mesh with gaps or overlaps at the cut. That's fine — it tells us where the mapping needs adjustment.

## Attempt 1 result: FAILED

**What happened:** Two crumpled fans, one per handle. No toroidal shape.

**Why:** The (ρ, φ) → (u, v) mapping is fundamentally wrong:
- ρ goes 0 → max (one-way), but v needs to be periodic (0 → 2π → 0). The periodicity comes from the b ↔ b⁻¹ side identification, which polar coordinates ignore entirely.
- φ creates a radial fan, not a tube wrap.

**Lesson:** Polar coordinates don't encode the torus topology. The torus has TWO independent periodicities (u and v), both created by side identifications. Polar coordinates have one periodicity (φ) and one non-periodic direction (ρ).

## Attempt 2 plan: Side-aligned (s, t) coordinates

Use coordinates aligned with the quadrilateral sides:
- **s**: interpolation from side aₖ to side aₖ⁻¹ → maps to u (big circle)  
- **t**: interpolation from side bₖ to side bₖ⁻¹ → maps to v (tube)

For each point P in a quadrilateral with corners A, B, C, D:
1. Project P onto the a-direction (A→B) to get s
2. Project P onto the b-direction (A→D) to get t
3. Map (s, t) → (2πs, 2πt) → torus

The corners of each quadrilateral are:
- Quad 1: V₀, V₁, V₂, V₃ (with V₀=V₄ cut edge) → but these are octagon vertices
  Wait — the quadrilateral has 5 sides (4 octagon sides + cut), not 4.
  
Hmm, the quadrilateral from cutting the octagon has sides: a₁, b₁, a₁⁻¹, b₁⁻¹, and the cut edge. That's 5 sides, not 4. It's a pentagon, not a quadrilateral.

Actually no: the cut goes from V₀ to V₄, and the octagon sides from V₀ to V₄ (going clockwise through V₁, V₂, V₃) are: a₁ (V₀→V₁), b₁ (V₁→V₂), a₁⁻¹ (V₂→V₃), b₁⁻¹ (V₃→V₄). So the quad has sides: cut (V₄→V₀), a₁ (V₀→V₁), b₁ (V₁→V₂), a₁⁻¹ (V₂→V₃), b₁⁻¹ (V₃→V₄). That's 5 sides.

Wait, but in the identification: a₁ ↔ a₁⁻¹ and b₁ ↔ b₁⁻¹. After identification of just these two pairs, the result is a torus-with-hole where the cut edge is the hole boundary.

So the relevant coordinates are:
- **s**: "distance" from the a₁ side to the a₁⁻¹ side (perpendicular to both)
- **t**: "distance" from the b₁ side to the b₁⁻¹ side (perpendicular to both)

In a regular octagon, these directions are well-defined. The a₁ side is at angle (let's say) 0°, and the a₁⁻¹ side is at 90° (two sides away). So s increases as you go from 0° toward 90°.

For a regular octagon with vertices at 45° intervals (-180°, -135°, ..., 135°):
- V₀ at -180° (= 180°)
- V₁ at -135°
- V₂ at -90°
- V₃ at -45°
- V₄ at 0°
- V₅ at 45°
- V₆ at 90°
- V₇ at 135°

Side a₁: V₀→V₁ (from -180° to -135°), midpoint at -157.5°
Side b₁: V₁→V₂ (from -135° to -90°), midpoint at -112.5°
Side a₁⁻¹: V₂→V₃ (from -90° to -45°), midpoint at -67.5°
Side b₁⁻¹: V₃→V₄ (from -45° to 0°), midpoint at -22.5°

The a₁↔a₁⁻¹ identification pairs angles -157.5° with -67.5° (90° apart).
The b₁↔b₁⁻¹ identification pairs angles -112.5° with -22.5° (also 90° apart).

So s measures "how far from a₁ toward a₁⁻¹" and t measures "how far from b₁ toward b₁⁻¹". Both are roughly angular, separated by 90°.

**Key idea**: in the hyperbolic octagon, s and t are like rotated/skewed coordinates. For a first approximation, I can define:
- s = projection onto the direction from midpoint(a₁) to midpoint(a₁⁻¹)
- t = projection onto the direction from midpoint(b₁) to midpoint(b₁⁻¹)

These directions are perpendicular (90° apart in the octagon). So (s,t) is just a rotation of (x,y) by some angle!

## Attempt 2 result: FAILED (same as attempt 1)

**What I see:** Two crumpled fans of triangles, one left and one right. No toroidal shape. Looks identical to attempt 1.

**Why it failed:** I projected disk coordinates onto two axis directions that are 45° apart (not 90°). Worse, I projected INDEPENDENTLY onto each axis and then mapped to torus angles. This is wrong because:

1. The s-axis (from midpoint(a₁) to midpoint(a₁⁻¹)) and t-axis (from midpoint(b₁) to midpoint(b₁⁻¹)) are at -22.5° and +22.5° — only 45° apart, not perpendicular.
2. Projecting onto non-perpendicular axes independently gives correlated coordinates — points cluster along a diagonal, not filling the (u,v) rectangle.
3. The linear projection doesn't account for the HYPERBOLIC geometry — the octagon's interior is hyperbolic, not Euclidean.

**Key insight:** The translation vectors τ_a = midpoint(a₁⁻¹) - midpoint(a₁) and τ_b = midpoint(b₁⁻¹) - midpoint(b₁) form a LATTICE. Points in the quadrilateral should be expressed in this lattice basis: p = p₀ + s·τ_a + t·τ_b, then solve for (s,t).

## Attempt 3 plan: Lattice-based (s,t) via matrix solve

For each handle's quadrilateral:
1. Compute τ_a and τ_b as complex vectors in the disk
2. Choose reference point p₀ = midpoint of side a₁
3. For each point w: solve w - p₀ = s·τ_a + t·τ_b (2×2 linear system)
4. Map (s,t) → (2πs, 2πt) → torus (u,v) → 3D

This properly decomposes the point into the two independent lattice directions.

## Attempt 3 result: FAILED (identical to attempt 2 — same crumpled fans)

The lattice solve gives s ∈ [-1, 0.414], t ∈ [0.29, 1.71]. Even after normalization, the shape won't be toroidal because the approach is fundamentally wrong: it treats the octagon as a parallelogram, ignoring the identification topology.

## Robin's insight: Map the 8 loops first

"A 2-torus has eight special loops which each correspond to a pair of edges. Try mapping them first then determining a relative coordinate system."

The 8 octagon sides become 4 loops on the double torus (each traversed twice, once in each direction):
- a₁, a₁⁻¹ → loop around big circle of torus 1
- b₁, b₁⁻¹ → loop around tube of torus 1
- a₂, a₂⁻¹ → loop around big circle of torus 2
- b₂, b₂⁻¹ → loop around tube of torus 2

All 8 octagon vertices merge to P (the neck between the two tori).

## Attempt 4 plan: Loop anchoring + Laplacian interpolation

**Step 1: Identify boundary vs interior vertices.**
- Boundary: vertices lying on the 8 geodesic sides of the octagon
- Interior: all other vertices

**Step 2: Map boundary vertices to 3D loops.**
Each octagon side becomes a loop from P back to P on the double torus:
- Side a₁ (V₀→V₁): parameterize position t ∈ [0,1] along the side
  → 3D curve: (-d + (R+r)cos(2πt), (R+r)sin(2πt), 0) (big circle of torus 1)
  At t=0 and t=1: position = P = (0,0,0) ✓
- Side b₁ (V₁→V₂): 
  → 3D curve: (-d + (R+r·cos(2πt)), 0, r·sin(2πt)) (tube of torus 1)
- Similarly for a₂, b₂ on torus 2.

Identified pairs (a₁ ↔ a₁⁻¹) share the same loop → merge vertices.

**Step 3: Solve Laplacian for interior vertices.**
Fix boundary vertices at their 3D loop positions.
For each interior vertex: position = average of neighbors' positions.
This is a sparse linear system — solve once.

**Why this should work:** Tutte's theorem guarantees that fixing the boundary of a planar graph to a convex polygon and solving the Laplacian gives a crossing-free embedding. Here the "boundary" is the octagon (mapped to loops), and the "interior" is filled by harmonic interpolation. The result should be a smooth, continuous surface.

## Attempt 4 result: FAILED (collapsed to origin)

**What happened:** All 3D positions = (0,0,0). The Laplacian converged in 1 iteration.

**Why:** Zero intermediate vertices on octagon sides! The octagon boundary has ONLY the 8 corner vertices (all identified to P = origin). With all boundary conditions = (0,0,0), the Laplacian gives all interior vertices = (0,0,0).

The (2,3,8) triangles don't subdivide the octagon sides — each side is a single edge from V_k to V_{k+1} with no intermediate points.

## Attempt 5 plan: Lattice coordinates with modular wrapping

Go back to lattice (s,t) coordinates but fix the wrapping.

**Key insight:** The octagon and the parallelogram are DIFFERENT TILES of the SAME lattice. Points in the octagon that have (s,t) outside [0,1]² are "borrowed" from adjacent parallelogram tiles. After taking s % 1 and t % 1, they wrap to the correct torus position.

**Critical fix for genus-2:** The octagon identification has 4 generators (a₁,b₁,a₂,b₂), creating a genus-2 surface — NOT a torus. So I CAN'T use a single (s,t) for the whole octagon. Each HANDLE has its own lattice:
- Handle 1: lattice from τ_{a₁}, τ_{b₁}
- Handle 2: lattice from τ_{a₂}, τ_{b₂}

Split the octagon by angle (upper vs lower semicircle), compute (s,t) per handle, wrap, and map to the corresponding torus.

## Attempts 5-7: Bilinear map + wounded tori

**What works:**
- Bilinear map from V₀V₁V₂V₃ → (0,0)(1,0)(1,1)(0,1) is exact at corners ✓
- Handle classification (24+24) is correct ✓  
- Gash curve computed correctly: V₀→V₄ maps to (0,0)→(-1.414, 0.773) ✓
- Wounded tori visualization with gash looks great ✓

**What still breaks:**
- Vertices near V₄ (the cut vertex) get (u,v) like (-1.4, 0.77). After mod 1: (0.586, 0.773).
- A triangle with one vertex at (0.0, 0.0) and another at (0.586, 0.773) has vertices on OPPOSITE sides of the torus — the triangle stretches across the whole surface as a giant shard.
- The mod 1 wrapping creates discontinuities: vertices that are close in the octagon end up far apart on the torus because they're on different "sheets" of the wrapping.

**Root cause:** The pentagon is NOT a parallelogram. It extends BEYOND the unit square [0,1]² in the bilinear parameterization. Vertices in the "overhang" (near V₄, with u < 0 or u > 1) get wrapped to the wrong side of the torus.

**The fix needed:** When a triangle has vertices that span a wrapping discontinuity (e.g., one vertex at u=0.05 and another at u=0.85 after mod 1), the triangle should be drawn with the CLOSER copy — i.e., one vertex should use u=0.05 and the other u=0.85-1=-0.15, so they're close together. This is the same problem as drawing on a periodic domain (like a world map near the date line).

## Attempt 8: R⁴ algebraic curve embedding via Poincaré series

**New approach:** Abandon R³ entirely. The Bolza surface is the algebraic curve y² = x⁵ − x, which lives naturally in C² = R⁴ with no singularities. Map each point z in the octagon to (x(z), y(z)) ∈ C², giving R⁴ coordinates (Re x, Im x, Re y, Im y). Visualize by projecting R⁴ → R³.

**Key mathematical idea:** For the hyperelliptic curve y² = x⁵ − x, the two holomorphic differentials are ω₁ = dx/(2y) and ω₂ = x·dx/(2y). Their ratio is x. So x(z) = Θ₂(z)/Θ₁(z) (up to Möbius), where Θ₁, Θ₂ are weight-2 Poincaré series over the Fuchsian group Γ = ker(Δ(2,3,8) → GL(2,3)).

**The Poincaré series:** For a pole point p ∈ H (outside the fundamental domain):

    Θ_p(z) = Σ_{γ ∈ Γ} 1 / ((a_γ − p·c_γ)z + (b_γ − p·d_γ))²

where γ = [[a,b],[c,d]] ∈ SL(2,R). This converges absolutely because the 1/(z−p)² seed provides extra decay beyond the weight-2 automorphy factor. Two series with different poles p₁, p₂ give linearly independent weight-2 forms; their ratio is a degree-2 meromorphic function on the surface (the hyperelliptic projection x, up to Möbius).

**Calibration:** The 6 Weierstrass points (= 6 order-8 vertices of the {8,3} cell complex = cosets of ⟨c⟩ in GL(2,3)) should map to x ∈ {0, 1, −1, i, −i, ∞}. Evaluate the ratio Θ₂/Θ₁ at these 6 points, then fit a Möbius transformation to match the known x-values. Try all 2 × 4! = 48 possible assignments, pick the one with lowest error.

**Implementation:** Files `bolza_group.py` (shared GL(2,3) + SL(2,R) infrastructure, extracted from the 4 duplicate copies) and `bolza_r4.py` (main pipeline).

**Results:**
- 1101 Γ-elements enumerated from 50,000 Δ-elements via BFS
- Systole trace: 4.8284 = 2 + 2√2, exact match ✓
- 6 Weierstrass points correctly identified ✓
- Automorphy check: max |x(γz) − x(z)| = 7.5 × 10⁻², mean 4.7 × 10⁻² (∼5% error)
- Degree check: preimages come in pairs (distances 0.0177, 0.0177), confirming degree 2 ✓
- Möbius calibration total error: 1.73 (summed over 5 finite Weierstrass points)
- R⁴ coordinate ranges: Re(x) ∈ [−11, 7], Re(y) ∈ [0, 707] (huge spike near x = ∞)

**What works:**
- The Poincaré series machinery is correct — automorphy holds to ∼5%
- The degree-2 structure of x is confirmed (preimage pairs)
- y² = x⁵ − x holds to machine precision (but this is tautological: y = √(x⁵−x))
- No NaN/Inf vertices; all 2160 mesh points are finite
- Subdivided mesh: 48 × 64 = 3072 sub-triangles, 2160 vertices

**What doesn't work:**
- Möbius calibration error of 1.73 is too high — the assignment of x-values to Weierstrass points may be wrong, or the series needs more terms
- R⁴ coordinates are extremely non-uniform: y ∼ |x|^{5/2}, so near x = ∞ the coordinates blow up to ∼700. arcsinh normalization + centering + rescaling helps but produces a distorted shape
- Visualization (Three.js) not yet verified — headless Playwright can't render WebGL
- The 5% automorphy error suggests 1101 Γ-elements may not be enough. Convergence is slow for weight-2 series (marginal exponent)

**Open questions:**
1. Is 5% automorphy error acceptable, or does the series need fundamentally more terms?
2. The assignment W2 → x = ∞ puts the pole at the center of the octagon — every triangle has large y-values near the center. Would a different normalization of the curve (e.g., degree-6 model y² = f₆(x) with no point at infinity) avoid the coordinate blowup?
3. Should we use weight-4 Poincaré series (which converge faster) instead of weight-2, then extract x from the Sym²(S₂) decomposition?
4. Would the Schottky uniformization approach (products over a free group on 2 generators) converge faster than the Poincaré series approach?

**Root cause analysis:** The 5% error is NOT because we need more terms. With 1101 Γ-elements, the tail of the series contributes < 0.1% (terms decay as |γ|⁻⁴, counting grows as |γ|², so convergence goes as 1/N). The actual problem: the Poincaré series converges to different accuracy at different points because evaluating at γ₀(z) requires terms near γ₀⁻¹·(neighborhood of pole), which may not be in our 1101-element truncation. The error is in the AUTOMORPHY CHECK, not the series itself. An "adaptive summation" (including the most relevant Γ-elements for each evaluation point) would fix this, but is complex to implement.

## Attempt 9 plan: Discrete Hodge theory (finite element, no infinite sums)

**Core idea:** Instead of reconstructing the holomorphic differentials from an infinite Poincaré series, solve for them directly as a finite element problem on the octagon mesh. A holomorphic differential φ(z) satisfies:
- Δφ = 0 in the octagon interior (holomorphic ⟹ harmonic in the Euclidean metric, since the Poincaré disk is conformally flat)
- φ(T_k z) = (c_k z + d_k)² · φ(z) on identified boundary pairs

where T_k are the 4 side-pairing transformations (generators of Γ acting as SL(2,R) matrices).

This is a LINEAR eigenvalue-0 problem for the Laplacian with twisted boundary conditions. The solution space has dimension g = 2 (the two holomorphic differentials). No infinite sums, no group enumeration, no convergence issues. Accuracy depends only on mesh resolution.

**Why this works:**
- Holomorphicity = harmonicity on a conformally flat domain (the Poincaré disk). The Cauchy-Riemann equations become Δu = 0, Δv = 0 where φ = u + iv.
- The octagon boundary identifications with their automorphy twists (c_k z + d_k)² turn this into a standard FEM problem with periodic boundary conditions and a complex phase twist.
- The kernel of the resulting sparse linear system is exactly the 2D space of holomorphic differentials.
- Their ratio φ₂/φ₁ gives x (up to Möbius), same as attempt 8.

**What's needed that we don't have:**
1. **Side-pairing transformations** T_{a₁}, T_{b₁}, T_{a₂}, T_{b₂} as SL(2,R) matrices. These map each octagon side to its identified partner. They are generators of Γ = π₁(Bolza). Can find them by:
   - Searching the 1101 Γ-elements from attempt 8 for the ones that map V₀ → V₃ and V₁ → V₂ (for the a₁ pairing), etc.
   - Or computing them directly from the octagon vertex geometry.
2. **Fine mesh with identified boundary vertices.** Subdivide the 48 triangles, then MERGE identified boundary vertex pairs (rather than keeping them separate as in attempt 8). The twist factors relate φ-values at merged pairs.
3. **Cotangent-weight Laplacian** on the subdivided mesh. Standard FEM construction.

**Advantages over Poincaré series (attempt 8):**
- No infinite sum — finite linear algebra (sparse eigenvalue solve)
- No Fuchsian group enumeration beyond the 4 generators
- Error controlled by mesh refinement, which is straightforward
- Naturally handles the interior of the octagon uniformly (no position-dependent convergence)

**Risks:**
- Finding the side-pairing transformations correctly (orientation matters)
- The twist factors (c_k z + d_k)² couple Re(φ) and Im(φ), making the linear system complex
- The pole of x (at one Weierstrass point) creates large values on the mesh — same visualization issue as attempt 8
- Need to verify that the discrete Laplacian's kernel is actually 2-dimensional (mesh quality matters)

**The coordinate blowup problem (shared with attempt 8):**
The curve y² = x⁵ − x has a point at infinity where |x| → ∞, |y| → ∞. Any R⁴ embedding using (Re x, Im x, Re y, Im y) must deal with this. Options:
- Stereographic projection on x: map x ∈ C∪{∞} to S² via (2Re x/(1+|x|²), 2Im x/(1+|x|²), (|x|²−1)/(|x|²+1)), then embed using (S² coords, bounded y). This gives a BOUNDED embedding.
- Two-chart approach: use (x, y) where |x| < C and (1/x, y/x³) where |x| ≥ C, with smooth interpolation.
- Accept the blowup and use arcsinh normalization (attempt 8's approach — ugly but functional).

## Attempt 9 result: FAILED

Ran but hit three interconnected issues:
1. Side-pairing matrices for a₁ and b₁ had det = 3 (should be 1) — the assumed octagon vertex labeling was wrong
2. Only 16/2145 vertices classified as boundary (used radius threshold, but octagon sides curve inward)
3. Two different Weierstrass points (W0, W5) mapped to the same mesh vertex (W5 was outside the octagon at |w| = 0.897)

**Root cause:** All three issues stem from incorrect assumptions about the fundamental domain geometry, specifically about which octagon vertices correspond to which side-pairing identifications.

## Attempt 10: GL(2,3) labels as the single source of truth

**Diagnosis of why attempts 8–9 failed:**

The fundamental domain, side pairings, and Weierstrass points were all derived from geometric heuristics (closest 48 triangles, vertex angle ordering, coset orbits). These heuristics produced inconsistent data:
- 9 distinct order-8 vertex positions in the disk (should map to 6 on the surface)
- 16 triangles meeting at the center (should be 8 for an order-8 vertex)
- Weierstrass points at |w| > R_oct (outside the octagon entirely)

**The fix:** Use the GL(2,3) labels — already computed for each triangle — as the SOLE source of structural data. No geometric guessing.

### Key mathematical insights

**Side pairings from GL(2,3):**
A boundary edge separates a fund domain triangle (label g₁, matrix M₁) from an exterior triangle (label g₂, matrix M₂). The side-pairing γ maps the exterior triangle to its fund domain equivalent:

    γ = M_fund(g₂) · M_ext(g₂)⁻¹

Since both matrices map to the same GL(2,3) element g₂, their ratio maps to the identity — confirming γ ∈ Γ. This requires keeping the SL(2,R) matrices of boundary-adjacent exterior triangles, which the BFS already generates but currently discards.

**Weierstrass points are NOT order-8 vertices:**
The 6 Weierstrass points are fixed points of the hyperelliptic involution ι = c⁴ ∈ GL(2,3). On the disk, c⁴ acts as z → −z. The center (z = 0) is the only fixed point of c⁴ itself. The other 5 are at positions where some γ ∈ Γ satisfies γ(−z) = z. Each such γ gives a Weierstrass point at:

    z_W = fixed_point(M_γ · c_real⁴)

This is an order-2 isometry with a unique fixed point in H.

**The 16-triangle anomaly at the center:**
The "closest 48" heuristic grabs triangles from both the identity Gamma-coset AND the c⁴-coset at the center (since c⁴ fixes the center, both cosets have triangles at distance 0). This corrupts the fundamental domain. Fix: use GL(2,3) labels to select exactly one triangle per label, breaking ties consistently.

### Implementation plan

1. **BFS with full data retention:** Generate ~4000 triangles, keeping SL(2,R) matrices AND GL(2,3) labels for ALL of them (not just the 48 closest)

2. **Proper fundamental domain:** For each GL(2,3) label, select the triangle whose ORDER-2 VERTEX (wa) is closest to the origin. This gives exactly 48 triangles, one per label, with no duplication. Verify: exactly 8 triangles share each order-8 vertex position.

3. **Side pairings from boundary adjacency:** For each edge on the octagon boundary, find the adjacent exterior triangle, look up its GL(2,3) label, find the corresponding fund domain triangle, and compute γ = M_fund · M_ext⁻¹. This gives 4 generators of Γ (one per identified side pair).

4. **Weierstrass points:** Compute M_γ · c_real⁴ for the 4 generators (and their products), find fixed points, check which are inside the fundamental domain. Should give 6 distinct points.

5. **Subdivide + FEM:** With correct identifications, build the discrete Laplacian and solve for holomorphic differentials (attempt 9's approach, now with correct boundary data).

6. **Stereographic embedding:** Map x → S² (bounded), avoiding the coordinate blowup.

### Attempt 10 result: PARTIAL SUCCESS

The algebraic side-pairing approach works:
- **14 Gamma generators found** (all with |tr| = 4.8284 = systole trace, det = 1) ✓
- **80 boundary pairs matched** via the Gamma elements ✓
- **Laplacian eigenvalue solve** ran successfully; spectrum shows lambda_0 ≈ 0 (constant kernel) ✓

Issues remaining:
- 17 "Weierstrass points" found (should be 6) — the Gamma-orbit deduplication doesn't account for identified octagon vertices and points outside the fundamental domain
- Some fund domain triangles still have wc outside the octagon (8/48 at |w| = 0.84)
- The eigenvalue gap (1e-3) between the constant and the next eigenvector suggests the boundary identification may be incomplete
- The Laplacian kernel should be 2D on a properly identified genus-2 surface (one constant + one from the topology); getting only 1D kernel means the surface isn't fully "closed"

**Key insight that broke the logjam:** Computing neighbors ALGEBRAICALLY from the GL(2,3) multiplication table, not geometrically from vertex proximity. For fund triangle (label g, matrix M), the neighbor via generator 'a' has label g·a and EXPECTED matrix M·a_real. The ACTUAL fund domain triangle with that label has matrix M'. The ratio gamma = M' · (M·a_real)⁻¹ is the side-pairing element. This is exact — no tolerance matching needed.

**Critical mathematical error discovered:** The scalar Laplacian on functions has kernel = {constants} on any compact surface. The holomorphic differentials live in the kernel of the d-bar operator (or the Hodge Laplacian on 1-forms), NOT the scalar Laplacian. The eigenvalue spectrum (lambda_1 through lambda_5 all at ~1e-3 with multiplicity 5) is actually the FIRST LAPLACIAN EIGENVALUE of the Bolza surface — confirming the boundary identification works, but solving the wrong equation.

Full mathematical analysis written to MATHEMATICS.md.

## Attempt 11 plan: Period matrix from the algebraic side

**Strategy shift:** Stop fighting the hyperbolic side. Start from the algebraic curve y² = x⁵ − x, where everything is explicit, and work TOWARD the Poincaré disk.

**Step 1: Compute the period matrix Omega.**
The holomorphic differentials omega_1 = dx/(2y) and omega_2 = x·dx/(2y) are explicit. Choose a canonical homology basis {A1, A2, B1, B2} by connecting branch points {0, 1, -1, i, -i, infinity} with branch cuts. Integrate omega_k along each cycle using Gauss-Jacobi quadrature (handles the sqrt singularities at branch points). The period matrix is Omega = A^{-1} B where A_{jk} = integral_{A_k} omega_j, B_{jk} = integral_{B_k} omega_j.

**Step 2: Compute the Riemann theta function.**
Theta(z; Omega) = sum_{n in Z^2} exp(pi*i * n^T Omega n + 2*pi*i * n^T z). For genus 2, this is a double sum that converges exponentially fast (each term decays as exp(-pi * |n|^2 * Im(Omega))).

**Step 3: Abel-Jacobi map on the algebraic curve.**
For each point P = (x, y) on the curve, compute u(P) = integral_{P0}^{P} (omega_1, omega_2) by numerical path integration along the curve. This gives a map from the curve to the Jacobian J = C^2 / (Z^2 + Omega Z^2).

**Step 4: Inverse Abel-Jacobi map via theta functions.**
Given a point u in the Jacobian, recover x(u) using the Thomae-type formula or by solving Theta(u - Abel(P); Omega) = 0 for P. This inverts Step 3.

**Step 5: Connect to the Poincaré disk.**
The Abel-Jacobi map u: curve -> Jacobian lifts to a map u_tilde: H -> C^2 on the universal cover. On H, u_tilde is holomorphic with u_tilde(gamma(z)) = u_tilde(z) + period(gamma) for gamma in Gamma. The Gamma generators (found in Attempt 10) map to LATTICE TRANSLATIONS in C^2.

To find u_tilde(z) for a disk point z:
- Evaluate u_tilde at the Weierstrass points (known in both pictures)
- Integrate the holomorphic differentials (as weight-2 forms) along paths in H from one Weierstrass point to z
- Or: use the Poincaré series from Attempt 8 to evaluate phi_1(z), phi_2(z), then u_tilde(z) = integral phi(z) dz

**Step 6: Compose the maps.**
For a disk point z: z -> u_tilde(z) -> (x, y) via inverse Abel-Jacobi. This gives the R^4 embedding.

**Why this should work:**
- Steps 1–4 are entirely on the algebraic side — no Fuchsian group needed, standard numerical methods
- Step 5 is the only place where the two pictures meet, and it can be done pointwise
- The period matrix for y² = x⁵ − x may have a known closed form (the Bolza surface has CM by Z[i])
- Genus 2 means 2×2 matrices and double sums — computationally tiny

**Risk:** Step 5 still requires evaluating holomorphic differentials on H (the same Poincaré series problem). But now we have the period matrix as a CROSS-CHECK: the series must reproduce the known periods. This constrains the error and gives a validation criterion that Attempt 8 lacked.

### Attempt 11 result: STUCK on canonical homology basis

Computed period integrals numerically (trapezoidal rule with complex branch tracking). Got finite A-periods and B-periods, but the resulting Omega is NOT symmetric (error 0.98) and Im(Omega) has a negative eigenvalue. Theta function explodes to 10^38.

**Root cause:** The cycles I chose don't form a canonical homology basis satisfying A_j · B_k = delta_{jk}. Getting this right requires:
1. Correct branch point ordering
2. Branch cuts that don't cross each other
3. A-cycles that encircle exactly one cut
4. B-cycles that connect the correct pairs of cuts WITHOUT crossing other cuts
5. Correct intersection numbers verified

This is a well-known difficulty in computational algebraic geometry. Standard algorithms (Deconinck-van Hoeij 2001, Tretkoff-Tretkoff) handle it automatically, but implementing them from scratch is non-trivial.

**Options:**
1. Use SageMath's `hyperelliptic_curves` module (has period matrix computation built in)
2. Implement Deconinck-van Hoeij's algorithm (proper cycle construction + Gauss quadrature)
3. Fix my cycles by computing intersection numbers and applying a symplectic correction
4. Abandon the period matrix route and go back to fixing Attempt 8's visualization
