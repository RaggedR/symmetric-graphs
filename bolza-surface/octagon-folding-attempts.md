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
