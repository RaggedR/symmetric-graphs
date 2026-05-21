"""
Explore the CSS code structure of GH(2,2) via the parabolic quadric Q(6,2).

Points of Q(6,2): nonzero vectors x in F_2^7 with Q(x) = x0*x4 + x1*x5 + x2*x6 + x3^2 = 0
These are the 63 points of the split Cayley hexagon, with Aut = G_2(2).
"""
import numpy as np

def quadratic_form(x):
    """Q(x) = x0*x4 + x1*x5 + x2*x6 + x3^2 over F_2"""
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

# Enumerate all nonzero vectors in F_2^7
points = []
for i in range(1, 128):  # 1 to 127 (skip zero vector)
    x = [(i >> j) & 1 for j in range(7)]
    if quadratic_form(x) == 0:
        points.append(x)

print(f"Number of points on Q(6,2): {len(points)}")
assert len(points) == 63, f"Expected 63, got {len(points)}"

# Form the 7 x 63 parity-check matrix H (points as columns)
H = np.array(points, dtype=int).T  # 7 x 63
print(f"H shape: {H.shape}")
print(f"Rank of H over F_2: {np.linalg.matrix_rank(H % 2)}")

# Check H * H^T mod 2
HHT = (H @ H.T) % 2
print(f"\nH * H^T mod 2:")
print(HHT)
print(f"\nH * H^T = 0 mod 2? {np.all(HHT == 0)}")

# If not zero, try the adjacency matrix of the collinearity graph instead
# Two points are collinear in GH(2,2) iff they are on a common line of the hexagon
# In Q(6,2), two points x,y are collinear on the quadric iff B(x,y) = 0
# where B is the associated bilinear form

def bilinear_form(x, y):
    """B(x,y) = x0*y4 + x4*y0 + x1*y5 + x5*y1 + x2*y6 + x6*y2 over F_2
    (polarization of Q, noting x3^2 term vanishes in char 2)"""
    return (x[0]*y[4] + x[4]*y[0] + x[1]*y[5] + x[5]*y[1] + x[2]*y[6] + x[6]*y[2]) % 2

# Build adjacency on Q(6,2) via perpendicularity (B(x,y) = 0)
n = len(points)
perp = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if i != j and bilinear_form(points[i], points[j]) == 0:
            perp[i][j] = 1

print(f"\nPerp graph: {n} vertices")
degrees = perp.sum(axis=1)
print(f"Degree sequence: min={degrees.min()}, max={degrees.max()}")
if np.all(degrees == degrees[0]):
    print(f"Regular of degree {degrees[0]}")

# Check perp adjacency matrix squared mod 2
A2 = (perp @ perp.T) % 2
print(f"\nA * A^T mod 2 = 0? {np.all(A2 == 0)}")
if not np.all(A2 == 0):
    print(f"Nonzero entries: {np.count_nonzero(A2)}")
    print(f"Diagonal of A^2 mod 2: {np.diag(A2)}")

# Also check: lines of the hexagon
# In GH(2,2), lines are specific totally isotropic lines on Q(6,2)
# (not ALL t.i. lines, only a G_2(2)-orbit of them)
# For now, check what we get from the full perp relation
print(f"\n--- Summary ---")
print(f"Q(6,2) has {len(points)} points (expected 63)")
print(f"Perp graph is {degrees[0]}-regular")
print(f"H (7x63 coordinate matrix): H*H^T = 0? {np.all(HHT == 0)}")
print(f"A (63x63 perp adjacency):   A*A^T = 0? {np.all(A2 == 0)}")
