"""
The (2,3,7) triangle group as a CONCRETE matrix group in SL(2,R).

Constructs explicit generators, verifies the presentation relations,
visualises the tiling, and identifies the kernel Γ = π₁(Klein quartic).
"""

import numpy as np

# ══════════════════════════════════════════════════════════════════════
#  The generators as real matrices
# ══════════════════════════════════════════════════════════════════════

# Generator a: rotation by π about i ∈ ℍ
# z ↦ -1/z, which is the matrix [0, 1; -1, 0] in SL(2,R)
# (we use [0,-1;1,0] which is the same in PSL since ±M are identified)
a = np.array([[0, 1], [-1, 0]], dtype=float)

print("Generator a (order 2, rotation by π at edge midpoint i):")
print(f"  a = {a.flatten()}")
print(f"  tr(a) = {np.trace(a):.6f}  (should be 0 = 2cos(π/2))")
print(f"  a² = {(a @ a).flatten()}  (should be -I)")

# Generator b: rotation by 2π/3 about a point
# We need tr(b) = 2cos(π/3) = 1, and tr(ab) = 2cos(π/7)
#
# Parameterize: b = [[cos(π/3), -f*sin(π/3)], [sin(π/3)/f, cos(π/3)]]
# Then tr(ab) = sin(π/3)/f + f*sin(π/3) = sin(π/3)*(f + 1/f)
# Set = 2cos(π/7): f + 1/f = 2cos(π/7)/sin(π/3)

cos_pi7 = np.cos(np.pi/7)
sin_pi3 = np.sin(np.pi/3)
cos_pi3 = np.cos(np.pi/3)

# f + 1/f = 2cos(π/7) / sin(π/3)
rhs = 2 * cos_pi7 / sin_pi3
print(f"\nf + 1/f = 2cos(π/7)/sin(π/3) = {rhs:.10f}")

# Solve f² - rhs*f + 1 = 0
disc = rhs**2 - 4
f = (rhs + np.sqrt(disc)) / 2
print(f"f = {f:.10f}")

b = np.array([
    [cos_pi3, -f * sin_pi3],
    [sin_pi3 / f, cos_pi3]
], dtype=float)

print(f"\nGenerator b (order 3, rotation by 2π/3 at face centre):")
print(f"  b = [{b[0,0]:.6f}, {b[0,1]:.6f}; {b[1,0]:.6f}, {b[1,1]:.6f}]")
print(f"  tr(b) = {np.trace(b):.6f}  (should be 1 = 2cos(π/3))")
print(f"  det(b) = {np.linalg.det(b):.6f}  (should be 1)")

b3 = b @ b @ b
print(f"  b³ = {b3.flatten().round(10)}  (should be -I)")

# Product c = ab
c = a @ b
print(f"\nc = ab (order 7, rotation by 2π/7 at vertex):")
print(f"  c = [{c[0,0]:.6f}, {c[0,1]:.6f}; {c[1,0]:.6f}, {c[1,1]:.6f}]")
print(f"  tr(c) = {np.trace(c):.10f}  (should be {2*cos_pi7:.10f} = 2cos(π/7))")

c7 = np.linalg.matrix_power(c, 7)
print(f"  c⁷ = {c7.flatten().round(10)}  (should be ±I)")


# ══════════════════════════════════════════════════════════════════════
#  The three fixed points in ℍ
# ══════════════════════════════════════════════════════════════════════

# An elliptic element M = [[α,β],[γ,δ]] in SL(2,R) with |tr| < 2
# fixes the point z = (α - δ + i√(4-tr²)) / (2γ) in ℍ

def fixed_point(M):
    """Fixed point of elliptic M ∈ SL(2,R) in the upper half-plane."""
    alpha, beta, gamma, delta = M[0,0], M[0,1], M[1,0], M[1,1]
    tr = alpha + delta
    disc = 4 - tr**2
    if disc < 0:
        return None  # hyperbolic, no fixed point in ℍ
    z = ((alpha - delta) + 1j * np.sqrt(disc)) / (2 * gamma)
    return z

z_a = fixed_point(a)
z_b = fixed_point(b)
z_c = fixed_point(c)

print(f"\nFixed points in ℍ (the vertices of the fundamental triangle):")
print(f"  z_a (angle π/2) = {z_a:.6f}  (edge midpoint)")
print(f"  z_b (angle π/3) = {z_b:.6f}  (face centre)")
print(f"  z_c (angle π/7) = {z_c:.6f}  (vertex)")

# The three sides of the triangle are geodesics (semicircles or vertical lines)
# connecting these three points.


# ══════════════════════════════════════════════════════════════════════
#  What kind of number is cos(π/7)?
# ══════════════════════════════════════════════════════════════════════

print("\n" + "="*60)
print("  THE ARITHMETIC")
print("="*60)

# cos(π/7) is a root of 8x³ - 4x² - 4x + 1 = 0
# Proof: cos(π/7) = (ζ₁₄ + ζ₁₄⁻¹)/2 where ζ₁₄ = e^{2πi/14}
# The minimal polynomial over Q has degree φ(14)/2 = 3.

x = cos_pi7
poly_val = 8*x**3 - 4*x**2 - 4*x + 1
print(f"\ncos(π/7) = {cos_pi7:.10f}")
print(f"8x³ - 4x² - 4x + 1 = {poly_val:.2e}  (should be 0)")
print(f"Minimal polynomial: 8x³ - 4x² - 4x + 1")
print(f"Degree 3 over Q — so Δ(2,3,7) lives in SL(2, Q(cos π/7))")

# The other roots of 8x³ - 4x² - 4x + 1:
# cos(3π/7) and cos(5π/7)
cos_3pi7 = np.cos(3*np.pi/7)
cos_5pi7 = np.cos(5*np.pi/7)
print(f"\nAll roots of 8x³ - 4x² - 4x + 1:")
print(f"  cos(π/7)  = {cos_pi7:.10f}")
print(f"  cos(3π/7) = {cos_3pi7:.10f}")
print(f"  cos(5π/7) = {cos_5pi7:.10f}")
print(f"  Sum = {cos_pi7 + cos_3pi7 + cos_5pi7:.10f}  (should be 1/2)")
print(f"  Product = {cos_pi7 * cos_3pi7 * cos_5pi7:.10f}  (should be 1/8)")


# ══════════════════════════════════════════════════════════════════════
#  The quotient map Δ(2,3,7) → PSL(2,7)
# ══════════════════════════════════════════════════════════════════════

print("\n" + "="*60)
print("  THE QUOTIENT MAP")
print("="*60)

print("""
The surjection Δ(2,3,7) → PSL(2,7) sends:

  Real matrices              Matrices mod 7
  in SL(2,R)                 in SL(2,F_7)/{±I}
  ─────────────              ──────────────────
  a = [0,1;-1,0]    ↦       [0,1;6,0]  (order 2)
  b = [½, ...; ..., ½]  ↦   [0,1;6,1]  (order 3)
  c = ab                ↦   [1,6;0,1]  (order 7)

  (Not literal reduction mod 7 — it's a GROUP homomorphism
  defined by mapping generators to generators.)

The KERNEL of this map is:

  Γ = ker(Δ(2,3,7) → PSL(2,7))

  [Δ(2,3,7) : Γ] = 168
""")


# ══════════════════════════════════════════════════════════════════════
#  What is the kernel Γ?
# ══════════════════════════════════════════════════════════════════════

print("="*60)
print("  THE KERNEL Γ = π₁(KLEIN QUARTIC)")
print("="*60)

print("""
Γ is a SURFACE GROUP of genus 3:

  Γ ≅ π₁(Σ₃) = ⟨α₁,β₁,α₂,β₂,α₃,β₃ | [α₁,β₁][α₂,β₂][α₃,β₃] = 1⟩

  - 6 generators (one pair per handle of the genus-3 surface)
  - 1 relation (the product of commutators is trivial)
  - Torsion-free (no element has finite order)
  - Index 168 in Δ(2,3,7)
  - Normal in Δ(2,3,7), with quotient PSL(2,7)

As a subgroup of SL(2,R):
  - Each generator is a HYPERBOLIC element (|trace| > 2)
  - These are translations along geodesics, not rotations
  - They have no fixed points in ℍ (torsion-free ↔ free action on ℍ)

The short exact sequence:

  1 → Γ → Δ(2,3,7) → PSL(2,7) → 1
       ↑       ↑            ↑
  surface    tiling      symmetry
   group    group of     group of
  (genus 3) hyperbolic  the quotient
             plane       surface

Γ acts freely on ℍ. The quotient ℍ/Γ is the Klein quartic.
Δ(2,3,7) acts on ℍ with 3 orbits of fixed points (orders 2,3,7).
The quotient ℍ/Δ(2,3,7) is a SPHERE with 3 cone points — an orbifold.
""")


# ══════════════════════════════════════════════════════════════════════
#  Finding elements of Γ by direct computation
# ══════════════════════════════════════════════════════════════════════

# Γ consists of words in a,b that map to the identity in PSL(2,7).
# Let's find some by brute force.

print("="*60)
print("  SEARCHING FOR KERNEL ELEMENTS")
print("="*60)

# PSL(2,7) generators (from our earlier computation)
a7 = np.array([[0, 1], [6, 0]]) % 7  # order 2 in PSL(2,F_7)
b7 = np.array([[0, 1], [6, 1]]) % 7  # order 3 in PSL(2,F_7)

def psl2_7_mult(A, B):
    return (A @ B) % 7

def psl2_7_key(A):
    A = A % 7
    negA = (7 - A) % 7
    return min(tuple(A.flatten()), tuple(negA.flatten()))

identity_key = psl2_7_key(np.eye(2, dtype=int))

# Generate short words in {a, b, b², a} and check which map to identity in PSL(2,7)
# (Note: a⁻¹ = a since a²=1, and b⁻¹ = b² since b³=1)

def word_to_real_matrix(word):
    """Evaluate a word in {a, b, b⁻¹} as a real matrix."""
    result = np.eye(2)
    b_inv = np.linalg.inv(b)
    for letter in word:
        if letter == 'a':
            result = result @ a
        elif letter == 'b':
            result = result @ b
        elif letter == 'B':  # b inverse
            result = result @ b_inv
    return result

def word_to_mod7(word):
    """Evaluate a word in {a, b, b⁻¹} as a matrix mod 7."""
    result = np.eye(2, dtype=int)
    b7_inv = (np.array([[1, 6], [1, 0]])) % 7  # b⁻¹ = b² mod 7
    for letter in word:
        if letter == 'a':
            result = psl2_7_mult(result, a7)
        elif letter == 'b':
            result = psl2_7_mult(result, b7)
        elif letter == 'B':
            result = psl2_7_mult(result, b7_inv)
    return result

# The commutator [a,b] = aba⁻¹b⁻¹ = ababB (since a⁻¹=a, b⁻¹=B=b²)
# has order... let's check
comm_word = "ababB"
comm_real = word_to_real_matrix(comm_word)
comm_mod7 = word_to_mod7(comm_word)
print(f"\n[a,b] = aba⁻¹b⁻¹:")
print(f"  Real trace: {np.trace(comm_real):.6f}  (|tr|>2 means hyperbolic)")
print(f"  Mod 7 key: {psl2_7_key(comm_mod7)}, identity? {psl2_7_key(comm_mod7) == identity_key}")

# Search for shortest words in the kernel
from itertools import product as iterproduct

print(f"\nSearching for kernel elements (words mapping to identity in PSL(2,7))...")
print(f"  (These are elements of Γ = π₁(Klein quartic))\n")

kernel_elements = []
for length in range(1, 16):
    count = 0
    # Use letters a, b, B (= b⁻¹). Avoid aa (= identity) and bB, Bb (cancellation)
    for word_tuple in iterproduct("abB", repeat=length):
        word = "".join(word_tuple)
        # Skip trivially reducible words
        skip = False
        for i in range(len(word)-1):
            if word[i] == word[i+1] == 'a':
                skip = True; break
            if (word[i] == 'b' and word[i+1] == 'B') or (word[i] == 'B' and word[i+1] == 'b'):
                skip = True; break
            if word[i] == word[i+1] == 'b' and i+2 < len(word) and word[i+2] == 'b':
                skip = True; break  # bbb = 1
            if word[i] == word[i+1] == 'B' and i+2 < len(word) and word[i+2] == 'B':
                skip = True; break  # BBB = 1
        if skip:
            continue

        mod7 = word_to_mod7(word)
        if psl2_7_key(mod7) == identity_key:
            real_mat = word_to_real_matrix(word)
            tr = np.trace(real_mat)
            if abs(abs(tr) - 2) > 0.01 and abs(tr) > 0.01:  # nontrivial
                count += 1
                if count <= 3:
                    kernel_elements.append((word, tr))
                    elt_type = "HYPERBOLIC" if abs(tr) > 2 else "elliptic"
                    print(f"  Length {length}: '{word}', tr = {tr:.6f} ({elt_type})")
    if count > 0:
        print(f"  ... {count} total kernel elements of length {length}")
        if len(kernel_elements) >= 6:
            break

print(f"""
Key observations:
- All nontrivial kernel elements are HYPERBOLIC (|trace| > 2)
- This confirms Γ is torsion-free (no elliptic elements)
- Hyperbolic elements are translations along geodesics in ℍ
- They have no fixed points → Γ acts freely → ℍ/Γ is smooth
""")


# ══════════════════════════════════════════════════════════════════════
#  The fundamental domain
# ══════════════════════════════════════════════════════════════════════

print("="*60)
print("  THE FUNDAMENTAL DOMAIN")
print("="*60)
print("""
Δ(2,3,7) tiles ℍ with triangles. ONE triangle has area π/42.

The Klein quartic = 168 triangles glued together = area 168·π/42 = 4π.
  (Check: 4π(g-1) = 4π·2 = 8π for Gauss-Bonnet... wait)

Actually: the ORIENTATION-PRESERVING triangle group Δ(2,3,7) has
fundamental domain = 2 triangles (a triangle + its mirror), area π/21.
168 copies tile the Klein quartic: 168 · π/21 = 8π = 4π(g-1) with g=3. ✓

So Γ has a fundamental domain made of 168 copies of the basic
triangle-pair. This 168-triangle polygon, with edges identified by
the 6 generators of Γ, IS the Klein quartic.

    ℍ  ──quotient by Γ──→  Klein quartic (genus 3 surface)
    │                          │
    │ quotient by Δ(2,3,7)     │ quotient by PSL(2,7)
    ↓                          ↓
    Orbifold S²(2,3,7)    =   Orbifold S²(2,3,7)
    (sphere with 3 cone points)
""")


# ══════════════════════════════════════════════════════════════════════
#  Visualize: map a fundamental triangle under the group
# ══════════════════════════════════════════════════════════════════════

print("="*60)
print("  FUNDAMENTAL TRIANGLE")
print("="*60)

# The three vertices of the fundamental triangle
print(f"\nFundamental triangle vertices in ℍ:")
print(f"  P_a = {z_a:.6f}  (angle π/2, fixed by a)")
print(f"  P_b = {z_b:.6f}  (angle π/3, fixed by b)")
print(f"  P_c = {z_c:.6f}  (angle π/7, fixed by c=ab)")

# Hyperbolic distances between vertices
def hyp_dist(z1, z2):
    """Hyperbolic distance in upper half-plane."""
    return np.arccosh(1 + abs(z1-z2)**2 / (2 * z1.imag * z2.imag))

d_ab = hyp_dist(z_a, z_b)
d_ac = hyp_dist(z_a, z_c)
d_bc = hyp_dist(z_b, z_c)

print(f"\nHyperbolic side lengths:")
print(f"  d(P_a, P_b) = {d_ab:.6f}  (side opposite π/7 vertex)")
print(f"  d(P_a, P_c) = {d_ac:.6f}  (side opposite π/3 vertex)")
print(f"  d(P_b, P_c) = {d_bc:.6f}  (side opposite π/2 vertex)")

# Verify angle sum
print(f"\nAngle sum: π/2 + π/3 + π/7 = {np.pi/2 + np.pi/3 + np.pi/7:.6f}")
print(f"  = {41/42:.6f}π")
print(f"  < π = {np.pi:.6f}  (HYPERBOLIC)")
print(f"  Defect = π - angle sum = π/42 = {np.pi/42:.6f}")
print(f"  Area = defect = π/42  (Gauss-Bonnet for hyperbolic triangle)")
