#!/usr/bin/env python3
"""
Compute the gash curve: where does the octagon cut V₀→V₄ land on each torus?

For each pentagon (half of the octagon):
  - The 4 octagon sides form a torus (via a↔a⁻¹, b↔b⁻¹ identification)
  - The cut becomes a curve on this torus — the "gash"
  - Sample points along the cut, compute their torus (u,v) coordinates
"""

import numpy as np

# Regular octagon centered at origin, vertices at 45° intervals
# V₀ at -180°, V₁ at -135°, ..., V₇ at 135°
R = 1.0  # vertex distance from center
verts = [R * np.exp(1j * np.radians(-180 + k*45)) for k in range(8)]

print("Octagon vertices:")
for k, v in enumerate(verts):
    print(f"  V{k} = ({v.real:.4f}, {v.imag:.4f}) at {np.degrees(np.angle(v)):.0f}°")

# Lower pentagon: V₀, V₁, V₂, V₃, V₄ (the a₁b₁ half)
# Quadrilateral for bilinear map: V₀, V₁, V₂, V₃
Q = [verts[0], verts[1], verts[2], verts[3]]  # corners → (0,0), (1,0), (1,1), (0,1)

def bilinear_inverse(w, Q):
    """Find (s,t) such that w = (1-s)(1-t)Q0 + s(1-t)Q1 + st*Q2 + (1-s)t*Q3."""
    q = [np.array([c.real, c.imag]) for c in Q]
    p = np.array([w.real, w.imag])
    s, t = 0.5, 0.5
    for _ in range(30):
        fwd = (1-s)*(1-t)*q[0] + s*(1-t)*q[1] + s*t*q[2] + (1-s)*t*q[3]
        F = fwd - p
        if np.linalg.norm(F) < 1e-14: break
        dFs = -(1-t)*q[0] + (1-t)*q[1] + t*q[2] - t*q[3]
        dFt = -(1-s)*q[0] - s*q[1] + s*q[2] + (1-s)*q[3]
        J = np.column_stack([dFs, dFt])
        try: delta = np.linalg.solve(J, -F)
        except: break
        s += delta[0]; t += delta[1]
    return s, t

# Sample the cut: V₀ → V₄, straight line
print("\nGash curve (V₀ → V₄) in torus coordinates:")
print(f"  {'t':>5}  {'disk_x':>8} {'disk_y':>8}  {'u':>8} {'v':>8}")

N = 20
gash_uv = []
for i in range(N+1):
    t = i / N
    w = (1-t) * verts[0] + t * verts[4]  # point along cut
    s_torus, t_torus = bilinear_inverse(w, Q)
    gash_uv.append((s_torus, t_torus))
    print(f"  {t:5.2f}  {w.real:8.4f} {w.imag:8.4f}  {s_torus:8.4f} {t_torus:8.4f}")

print(f"\nGash starts at ({gash_uv[0][0]:.3f}, {gash_uv[0][1]:.3f})")
print(f"Gash passes through ({gash_uv[N//2][0]:.3f}, {gash_uv[N//2][1]:.3f}) at midpoint")
print(f"Gash ends at ({gash_uv[-1][0]:.3f}, {gash_uv[-1][1]:.3f})")
print(f"\nOn the torus, start and end should both be at (0,0) = P")
print(f"The midpoint (center of octagon) maps to the middle of the gash")
