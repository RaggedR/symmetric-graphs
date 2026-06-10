#!/usr/bin/env python3
"""Verify the algebraic characterization of the Langer graph (collinearity graph of GH(2,2)).

Points of Q(6,2) are nonzero (a,b) in F_2^3 x F_2^3, with x_3 = a·b.
Two points are collinear in GH(2,2) iff the Zorn vector matrix product M1*M2 = 0.

Zorn multiplication (Baez, "The Octonions"):
  [[α,a],[b,α]] · [[γ,c],[d,γ]] =
    [[ αγ + a·d,   γa + αc + b×d ],
     [ γb + αd + a×c,   αγ + b·c ]]

where α = a·b, γ = c·d (norm-zero / imaginary condition).

Conditions for product = 0:
  S1: αγ + a·d = 0           (scalar)
  V1: γa + αc + b×d = 0      (vector, F₂³)
  V2: γb + αd + a×c = 0      (vector, F₂³)
  S2: αγ + b·c = 0           (scalar; with S1, implies a·d = b·c = perpendicularity)
"""

from itertools import product as cart_product

F2 = [0, 1]

points = []
for bits in cart_product(F2, repeat=6):
    if all(b == 0 for b in bits):
        continue
    a = bits[:3]
    b = bits[3:]
    points.append((a, b))

assert len(points) == 63

def dot(a, b):
    return sum(a[i]*b[i] for i in range(3)) % 2

def cross(a, b):
    return (
        (a[1]*b[2] + a[2]*b[1]) % 2,
        (a[2]*b[0] + a[0]*b[2]) % 2,
        (a[0]*b[1] + a[1]*b[0]) % 2,
    )

def vec_add(a, b):
    return tuple((a[i] + b[i]) % 2 for i in range(3))

def scalar_mul(s, v):
    return tuple((s * v[i]) % 2 for i in range(3))

def is_zero_vec(v):
    return all(x == 0 for x in v)

def zorn_collinear(p1, p2):
    """Two points are collinear in GH(2,2) iff Zorn product M1*M2 = 0."""
    a, b = p1
    c, d = p2
    if p1 == p2:
        return False

    alpha = dot(a, b)  # = a·b
    gamma = dot(c, d)  # = c·d

    # S1: αγ + a·d = 0
    if (alpha * gamma + dot(a, d)) % 2 != 0:
        return False

    # S2: αγ + b·c = 0
    if (alpha * gamma + dot(b, c)) % 2 != 0:
        return False

    # V1: γa + αc + b×d = 0  (NOTE: coeff of a is γ, coeff of c is α)
    v1 = vec_add(vec_add(scalar_mul(gamma, a), scalar_mul(alpha, c)), cross(b, d))
    if not is_zero_vec(v1):
        return False

    # V2: γb + αd + a×c = 0  (NOTE: coeff of b is γ, coeff of d is α)
    v2 = vec_add(vec_add(scalar_mul(gamma, b), scalar_mul(alpha, d)), cross(a, c))
    if not is_zero_vec(v2):
        return False

    return True


# Build adjacency and check regularity
degrees = []
edges = 0
adj_matrix = [[False]*63 for _ in range(63)]
for i in range(63):
    deg = 0
    for j in range(63):
        if zorn_collinear(points[i], points[j]):
            adj_matrix[i][j] = True
            deg += 1
            if i < j:
                edges += 1
    degrees.append(deg)

print(f"Degree sequence: {sorted(set(degrees))}")
print(f"Min degree: {min(degrees)}, Max degree: {max(degrees)}")
print(f"Undirected edges: {edges}")
print(f"Expected: 6-regular, 189 edges")

if sorted(set(degrees)) == [6] and edges == 189:
    print("\n✓ CONFIRMED: Zorn product formula gives 6-regular, 189 edges")

    # Check distance-regularity: intersection array {6,4,4;1,1,3}
    # For adjacent vertices: should have exactly 1 common neighbor (λ=1)
    lambda_vals = set()
    for i in range(63):
        for j in range(i+1, 63):
            if adj_matrix[i][j]:
                common = sum(1 for k in range(63) if adj_matrix[i][k] and adj_matrix[j][k])
                lambda_vals.add(common)
    print(f"λ values (common neighbors of adjacent pairs): {lambda_vals}")

    # For non-adjacent vertices at distance 2: should have μ=1
    mu_vals = set()
    for i in range(63):
        for j in range(i+1, 63):
            if not adj_matrix[i][j]:
                common = sum(1 for k in range(63) if adj_matrix[i][k] and adj_matrix[j][k])
                if common > 0:  # distance 2
                    mu_vals.add(common)
    print(f"μ values (common neighbors of non-adj distance-2 pairs): {mu_vals}")

    if lambda_vals == {1} and mu_vals == {1, 3}:
        print("\n✓ CONFIRMED: Distance-regular with intersection array {6,4,4;1,1,3}")
        print("This is the collinearity graph of GH(2,2) = the Langer graph!")
    else:
        print(f"\n? Distance-regularity check: λ={lambda_vals}, μ={mu_vals}")
else:
    print("\n✗ FAILED: formula does not match expected Langer graph")
