#!/usr/bin/env python3
"""
R^4 embedding of the Bolza surface via the algebraic curve y^2 = x^5 - x.

Computes the uniformization map from the Poincare disk to C^2 = R^4
using weight-2 Poincare series over the Fuchsian group Gamma.

The key insight: for a genus-2 hyperelliptic curve y^2 = f(x), the two
holomorphic differentials are dx/(2y) and x*dx/(2y). Their ratio is x.
So x(z) = Theta_2(z) / Theta_1(z) (up to Mobius transformation), where
Theta_1, Theta_2 are weight-2 Poincare series.
"""

import numpy as np
import json
import sys
from collections import deque
from bolza_group import GL23, TriangleGroup, BolzaTiling


# ══════════════════════════════════════════════════════════════════════
#  Step 1: Enumerate Gamma = ker(Delta -> GL(2,3))
# ══════════════════════════════════════════════════════════════════════

def enumerate_gamma(gl23, tg, max_elements=50000):
    """BFS through Delta(2,3,8), keeping elements that map to identity in GL(2,3).

    Returns list of SL(2,R) matrices for elements of Gamma (the surface group).
    """
    gen_sl = [('a', tg.a_real), ('b', tg.b_real), ('B', tg.b_inv_real)]
    gen_gl = {'a': gl23.a_gen, 'b': gl23.b_gen, 'B': gl23.b_inv_gen}

    gamma_mats = []
    visited = set()
    queue = deque([(np.eye(2), gl23.identity)])

    while queue and len(visited) < max_elements:
        M_sl, g_gl = queue.popleft()
        mk = tg.mat_key(M_sl, tol=1e-8)
        if mk in visited:
            continue
        visited.add(mk)

        # Is this element in Gamma?
        if g_gl == gl23.identity:
            gamma_mats.append(M_sl.copy())

        for gn, gmat in gen_sl:
            new_M = M_sl @ gmat
            new_g = gl23.mult(g_gl, gen_gl[gn])
            queue.append((new_M, new_g))

    return gamma_mats


# ══════════════════════════════════════════════════════════════════════
#  Step 2: Weight-2 Poincare series
# ══════════════════════════════════════════════════════════════════════

def poincare_series(gamma_mats, z_points, pole):
    """Compute weight-2 Poincare series with pole seed at `pole`.

    Theta_p(z) = sum_{gamma in Gamma} 1 / ((a - p*c)*z + (b - p*d))^2

    where gamma = [[a,b],[c,d]] in SL(2,R) and p = pole (in the upper half-plane).

    Args:
        gamma_mats: list of SL(2,R) 2x2 arrays
        z_points: array of complex points in the upper half-plane
        pole: complex number in H (the seed pole)

    Returns:
        array of complex values, shape (len(z_points),)
    """
    z = np.asarray(z_points, dtype=complex)

    # Extract (a, b, c, d) for all gamma elements
    N_gamma = len(gamma_mats)
    abcd = np.zeros((N_gamma, 4), dtype=complex)
    for i, M in enumerate(gamma_mats):
        abcd[i] = [M[0,0], M[0,1], M[1,0], M[1,1]]

    a, b, c, d = abcd[:, 0], abcd[:, 1], abcd[:, 2], abcd[:, 3]

    # Vectorized: denom[i, j] = (a[i] - pole*c[i]) * z[j] + (b[i] - pole*d[i])
    coeff_z = (a - pole * c)[:, None]   # (N_gamma, 1)
    coeff_1 = (b - pole * d)[:, None]   # (N_gamma, 1)
    denom = coeff_z * z[None, :] + coeff_1  # (N_gamma, N_z)

    # Sum 1/denom^2 over gamma
    result = np.sum(1.0 / (denom * denom), axis=0)
    return result


# ══════════════════════════════════════════════════════════════════════
#  Step 3: Calibrate x via Weierstrass points
# ══════════════════════════════════════════════════════════════════════

def find_mobius_transform(r_values, x_targets):
    """Find Mobius transformation x = (A*r + B) / (C*r + D) mapping r_values to x_targets.

    Uses 3 point-pairs. Returns (A, B, C, D) normalized with AD - BC = 1.

    For r -> infinity (Theta_1 = 0), pass r = np.inf and x = the target value.
    For x -> infinity (a Weierstrass point at infinity), we handle it separately.
    """
    # Filter out infinities — use 3 finite pairs
    finite_pairs = [(r, x) for r, x in zip(r_values, x_targets)
                    if np.isfinite(r) and np.isfinite(x)]

    if len(finite_pairs) < 3:
        raise ValueError(f"Need 3 finite pairs, got {len(finite_pairs)}")

    (r1, x1), (r2, x2), (r3, x3) = finite_pairs[:3]

    # x = (A*r + B) / (C*r + D) => x*(C*r + D) = A*r + B
    # => x*C*r + x*D - A*r - B = 0
    # => r*(x*C - A) + (x*D - B) = 0
    #
    # Three equations:
    # r1*(x1*C - A) + (x1*D - B) = 0
    # r2*(x2*C - A) + (x2*D - B) = 0
    # r3*(x3*C - A) + (x3*D - B) = 0
    #
    # Set D = 1 (free normalization), solve for A, B, C:
    # r1*x1*C - r1*A + x1 - B = 0
    # r2*x2*C - r2*A + x2 - B = 0
    # r3*x3*C - r3*A + x3 - B = 0
    #
    # Matrix form: M @ [A, B, C]^T = [x1, x2, x3]^T
    M = np.array([
        [-r1, -1, r1*x1],
        [-r2, -1, r2*x2],
        [-r3, -1, r3*x3],
    ], dtype=complex)
    rhs = np.array([-x1, -x2, -x3], dtype=complex)
    ABC = np.linalg.solve(M, rhs)
    A, B, C = ABC
    D = 1.0 + 0j

    # Normalize: AD - BC = 1
    det = A * D - B * C
    scale = np.sqrt(det)
    return A/scale, B/scale, C/scale, D/scale


def apply_mobius(r, A, B, C, D):
    """Apply Mobius transformation x = (A*r + B) / (C*r + D)."""
    return (A * r + B) / (C * r + D)


# ══════════════════════════════════════════════════════════════════════
#  Step 4: Subdivide triangles
# ══════════════════════════════════════════════════════════════════════

def subdivide_triangle(w_a, w_b, w_c, N=8):
    """Subdivide a triangle in the Poincare disk using Euclidean barycentrics.

    Returns (vertices, sub_triangles) where vertices[i] is complex
    and sub_triangles[j] is (i0, i1, i2).
    """
    verts = []
    idx_map = {}
    for i in range(N + 1):
        for j in range(N + 1 - i):
            k = N - i - j
            w = (i * w_a + j * w_b + k * w_c) / N
            idx_map[(i, j)] = len(verts)
            verts.append(w)

    tris = []
    for i in range(N):
        for j in range(N - i):
            # Upward triangle
            v0 = idx_map[(i, j)]
            v1 = idx_map[(i + 1, j)]
            v2 = idx_map[(i, j + 1)]
            tris.append((v0, v1, v2))
            # Downward triangle (if it exists)
            if i + j + 1 < N:
                v3 = idx_map[(i + 1, j + 1)]
                tris.append((v1, v3, v2))

    return np.array(verts, dtype=complex), np.array(tris, dtype=int)


# ══════════════════════════════════════════════════════════════════════
#  Main pipeline
# ══════════════════════════════════════════════════════════════════════

def main():
    print("="*60)
    print("  Bolza Surface: R^4 Embedding via y^2 = x^5 - x")
    print("="*60)

    # --- Build group infrastructure ---
    print("\n[1/7] Building GL(2,3) and (2,3,8) triangle group...")
    gl23 = GL23()
    tg = TriangleGroup()
    tiling = BolzaTiling(gl23, tg)
    print(f"  48 triangles, {len(tiling.vertices)} unique vertices")
    print(f"  z_a = {tg.z_a:.6f} (order-2, in H)")
    print(f"  z_b = {tg.z_b:.6f} (order-3, in H)")
    print(f"  z_c = {tg.z_c:.6f} (order-8, in H)")

    # --- Enumerate Gamma ---
    print("\n[2/7] Enumerating Gamma = ker(Delta -> GL(2,3))...")
    MAX_DELTA = 50000
    gamma_mats = enumerate_gamma(gl23, tg, max_elements=MAX_DELTA)
    print(f"  {len(gamma_mats)} Gamma elements (from {MAX_DELTA} Delta elements)")

    # Verify: identity and -I should be in Gamma
    traces = [abs(np.trace(M)) for M in gamma_mats]
    n_identity = sum(1 for t in traces if abs(t - 2) < 1e-6)
    n_neg_I = sum(1 for t in traces if abs(t + 2) < 1e-6)
    print(f"  Identity-like: {n_identity}, -I-like: {n_neg_I}")
    print(f"  Trace range: [{min(traces):.4f}, {max(traces):.4f}]")

    # Systole check: smallest trace > 2 should be 2 + 2*sqrt(2) = 4.828
    hyperbolic_traces = sorted([t for t in traces if t > 2.1])
    if hyperbolic_traces:
        print(f"  Systole trace: {hyperbolic_traces[0]:.4f} (expected {2 + 2*np.sqrt(2):.4f})")

    # --- Weierstrass points ---
    print("\n[3/7] Identifying Weierstrass points...")
    wp_disk = tiling.weierstrass_points_disk()
    print(f"  {len(wp_disk)} Weierstrass points in the Poincare disk:")
    for i, w in enumerate(wp_disk):
        z = tg.disk_to_halfplane(w)
        print(f"    W{i}: disk = {w:.6f}, H = {z:.6f}")

    # Convert to upper half-plane for Poincare series evaluation
    wp_H = np.array([tg.disk_to_halfplane(w) for w in wp_disk])

    # --- Compute Poincare series at Weierstrass points ---
    print("\n[4/7] Computing weight-2 Poincare series at Weierstrass points...")

    # Choose poles well outside the fundamental domain
    # Pole 1: high on the imaginary axis (corresponds to disk point near boundary)
    pole1 = 19j  # corresponds to disk point ~0.9
    # Pole 2: different direction
    pole2 = -3 + 5j  # asymmetric, well outside

    theta1_W = poincare_series(gamma_mats, wp_H, pole1)
    theta2_W = poincare_series(gamma_mats, wp_H, pole2)

    print("  Theta_1 at Weierstrass points:")
    for i, (t1, t2) in enumerate(zip(theta1_W, theta2_W)):
        ratio = t2/t1 if abs(t1) > 1e-15 else np.inf
        print(f"    W{i}: |Theta_1| = {abs(t1):.6e}, |Theta_2| = {abs(t2):.6e}, ratio = {ratio:.6f}")

    # --- Calibrate the Mobius transform ---
    print("\n[5/7] Calibrating x = Mobius(Theta_2/Theta_1)...")

    # The 6 Weierstrass points map to {0, 1, -1, i, -i, infinity} on y^2 = x^5 - x.
    # We need to determine the assignment.
    #
    # Strategy: the hyperelliptic involution is c^4 (rotation by pi in the disk).
    # It fixes x (since x descends to CP^1 = C / <iota>).
    # The generator c acts as x -> i*x (4th root of unity, since c has order 8
    # and x -> -x is c^4 = iota, so c acts as a primitive 4th root on x).
    #
    # The two c-fixed Weierstrass points have x in {0, infinity} (the fixed points
    # of x -> i*x). The remaining 4 form a single orbit {1, i, -1, -i} under c.
    #
    # Identify c-fixed points: they are at the center (w=0) and at the
    # octagon vertex identification point.

    # Find the two Weierstrass points closest to the origin (c-fixed points)
    wp_dists = [(i, abs(wp_disk[i])) for i in range(6)]
    wp_dists.sort(key=lambda x: x[1])

    # The center (w ~ 0) is one c-fixed point
    center_idx = wp_dists[0][0]
    print(f"  Center Weierstrass point: W{center_idx} at |w| = {wp_dists[0][1]:.6f}")

    # The farthest point is the other c-fixed point (octagon vertex, x = infinity or 0)
    far_idx = wp_dists[-1][0]
    print(f"  Farthest Weierstrass point: W{far_idx} at |w| = {wp_dists[-1][1]:.6f}")

    # Determine which is x=0 and which is x=infinity by checking Theta_1.
    # At x=infinity, the holomorphic differential dx/(2y) has a zero
    # (since dx ~ -2/t^3 dt and 2y ~ 2/t^5 * ... near the point at infinity,
    # the differential dx/(2y) vanishes). This means Theta_1 should be SMALLER
    # there. But actually, our Poincare series with arbitrary poles aren't the
    # actual holomorphic differentials — they're a linear combination.
    # So we can't directly identify x=infinity this way.
    #
    # Instead, use a brute-force approach: try assigning x-values and check
    # which assignment makes y^2 = x^5 - x hold at other points.

    # Compute ratios at all Weierstrass points
    ratios = np.array([theta2_W[i] / theta1_W[i] if abs(theta1_W[i]) > 1e-15
                       else np.inf + 0j for i in range(6)])

    print("\n  Ratios Theta_2/Theta_1 at Weierstrass points:")
    for i, r in enumerate(ratios):
        print(f"    W{i}: {r:.6f}")

    # The 6 target x-values for the Bolza curve y^2 = x^5 - x = x(x-1)(x+1)(x-i)(x+i):
    x_targets_all = [0, 1, -1, 1j, -1j, np.inf]

    # Try all permutations of the c-orbit {1, i, -1, -i} among the 4 non-fixed points.
    # The c-fixed points get {0, infinity} in some order. That's 2 * 4! = 48 possibilities,
    # but symmetry reduces this.
    #
    # Practical approach: use 3 finite Weierstrass points to fit the Mobius map,
    # then validate at the remaining points by checking consistency.

    best_error = np.inf
    best_mobius = None
    best_assignment = None

    from itertools import permutations

    # The two c-fixed points get {0, infinity} — try both assignments
    non_fixed = [i for i in range(6) if i != center_idx and i != far_idx]

    for inf_choice in [center_idx, far_idx]:
        zero_choice = far_idx if inf_choice == center_idx else center_idx

        # The 4 non-fixed points get some permutation of {1, -1, i, -i}
        finite_targets = [1, -1, 1j, -1j]
        for perm in permutations(finite_targets):
            # Build assignment
            x_assign = {}
            x_assign[zero_choice] = 0
            x_assign[inf_choice] = np.inf
            for k, idx in enumerate(non_fixed):
                x_assign[idx] = perm[k]

            # Collect finite (ratio, x) pairs for Mobius fit
            finite_pairs_r = []
            finite_pairs_x = []
            for idx in range(6):
                if np.isfinite(ratios[idx]) and np.isfinite(x_assign[idx]):
                    finite_pairs_r.append(ratios[idx])
                    finite_pairs_x.append(x_assign[idx])

            if len(finite_pairs_r) < 3:
                continue

            try:
                A, B, C, D = find_mobius_transform(
                    finite_pairs_r[:3], finite_pairs_x[:3])
            except (np.linalg.LinAlgError, ValueError):
                continue

            # Validate: compute x at ALL finite Weierstrass points
            error = 0
            for idx in range(6):
                if not np.isfinite(ratios[idx]) or not np.isfinite(x_assign[idx]):
                    continue
                x_computed = apply_mobius(ratios[idx], A, B, C, D)
                error += abs(x_computed - x_assign[idx])

            if error < best_error:
                best_error = error
                best_mobius = (A, B, C, D)
                best_assignment = dict(x_assign)

    if best_mobius is None:
        print("  ERROR: Could not calibrate Mobius transform!")
        sys.exit(1)

    A, B, C, D = best_mobius
    print(f"\n  Best Mobius calibration (total error = {best_error:.6e}):")
    print(f"    x = ({A:.6f}*r + {B:.6f}) / ({C:.6f}*r + {D:.6f})")
    print(f"  Assignment:")
    for idx in range(6):
        x_val = best_assignment[idx]
        print(f"    W{idx} -> x = {x_val}")

    # --- Subdivide mesh and compute R^4 coordinates ---
    print("\n[6/7] Subdividing mesh and computing R^4 embedding...")

    SUB_N = 8  # subdivision level per triangle
    all_disk_pts = []     # complex disk coordinates
    all_sub_tris = []     # triangle indices (global)
    all_parent = []       # which of the 48 original triangles
    vertex_offset = 0

    for ti, (wa, wb, wc) in enumerate(tiling.fund_triangles):
        sub_verts, sub_tris = subdivide_triangle(wa, wb, wc, N=SUB_N)
        all_disk_pts.extend(sub_verts)
        all_sub_tris.extend(sub_tris + vertex_offset)
        all_parent.extend([ti] * len(sub_tris))
        vertex_offset += len(sub_verts)

    all_disk_pts = np.array(all_disk_pts)
    all_sub_tris = np.array(all_sub_tris)
    print(f"  {len(all_disk_pts)} vertices, {len(all_sub_tris)} triangles")

    # Convert disk points to upper half-plane
    print("  Converting to upper half-plane...")
    all_H_pts = np.array([tg.disk_to_halfplane(w) for w in all_disk_pts])

    # Evaluate Poincare series
    print("  Evaluating Poincare series (this may take a moment)...")
    theta1 = poincare_series(gamma_mats, all_H_pts, pole1)
    theta2 = poincare_series(gamma_mats, all_H_pts, pole2)

    # Compute x(z)
    ratio = theta2 / theta1
    x_vals = apply_mobius(ratio, A, B, C, D)

    # Compute y(z) = sqrt(x^5 - x), choosing branch by continuity
    f_x = x_vals**5 - x_vals
    y_vals = np.sqrt(f_x)

    # Check: flip branch where needed for continuity
    # (within each original triangle, the branch should be consistent)
    # For now, use numpy's default branch and verify

    # R^4 coordinates
    r4_coords = np.column_stack([x_vals.real, x_vals.imag, y_vals.real, y_vals.imag])
    print(f"  R^4 coordinate ranges:")
    for i, label in enumerate(['Re(x)', 'Im(x)', 'Re(y)', 'Im(y)']):
        finite_mask = np.isfinite(r4_coords[:, i])
        if finite_mask.any():
            vals = r4_coords[finite_mask, i]
            print(f"    {label}: [{vals.min():.4f}, {vals.max():.4f}]")

    # --- Validation ---
    # 1. Automorphy check: x(gamma(z)) should equal x(z) for gamma in Gamma
    print("\n  Validation 1: Automorphy check (x(gamma*z) == x(z))...")
    n_auto_tests = min(20, len(gamma_mats) - 1)
    # Pick a test point well inside the fundamental domain
    z_test = tg.disk_to_halfplane(0.1 + 0.1j)  # a generic interior point
    x_test = apply_mobius(
        poincare_series(gamma_mats, [z_test], pole2)[0] /
        poincare_series(gamma_mats, [z_test], pole1)[0],
        A, B, C, D)
    auto_errors = []
    for gi in range(1, n_auto_tests + 1):
        M = gamma_mats[gi]
        gz = tg.mobius(M, z_test)
        x_gz = apply_mobius(
            poincare_series(gamma_mats, [gz], pole2)[0] /
            poincare_series(gamma_mats, [gz], pole1)[0],
            A, B, C, D)
        auto_errors.append(abs(x_gz - x_test))
    print(f"    Tested {n_auto_tests} Gamma elements")
    print(f"    Max |x(gamma*z) - x(z)|: {max(auto_errors):.6e}")
    print(f"    Mean: {np.mean(auto_errors):.6e}")

    # 2. Degree check: x should take each generic value exactly twice
    print("\n  Validation 2: Degree check (x should have degree 2)...")
    # Count approximate preimages of x = 0.5 among mesh vertices
    for target in [0.5, 0.3+0.4j]:
        dists = np.abs(x_vals - target)
        n_close = np.sum(dists < 0.3)
        closest = np.sort(dists)[:4]
        print(f"    Preimages of x={target}: {n_close} within 0.3 (closest: {closest[:3]})")

    # Check for NaN/Inf
    n_nan = np.sum(~np.isfinite(r4_coords).all(axis=1))
    print(f"\n  Non-finite vertices: {n_nan} / {len(r4_coords)}")

    # --- Generate visualization ---
    print("\n[7/7] Generating Three.js visualization...")

    # Normalize R^4 coordinates using arcsinh to compress large values
    # arcsinh(x) ~ x for small x, ~ sign(x)*ln(2|x|) for large x
    scale = 2.0  # controls compression: larger = less compression
    r4_norm = np.arcsinh(r4_coords / scale) * scale
    r4_norm[~np.isfinite(r4_norm)] = 0.0

    # Center at origin and rescale so max radius ~ 3
    r4_center = r4_norm.mean(axis=0)
    r4_norm -= r4_center
    r4_max = np.max(np.linalg.norm(r4_norm, axis=1))
    if r4_max > 0:
        r4_norm *= 3.0 / r4_max

    print(f"  Normalized R^4 ranges (arcsinh + centered, scale={scale}):")
    for i, label in enumerate(['Re(x)', 'Im(x)', 'Re(y)', 'Im(y)']):
        vals = r4_norm[:, i]
        print(f"    {label}: [{vals.min():.4f}, {vals.max():.4f}]")

    # Color by arg(x) — hue wheel
    x_arg = np.angle(x_vals) / (2 * np.pi) % 1.0
    x_arg[~np.isfinite(x_arg)] = 0.0

    # Project R^4 -> R^3: start with orthographic (drop Im(y))
    # But also support 4D rotation via slider
    html = generate_html(r4_norm, all_sub_tris, x_arg, all_parent, all_disk_pts)

    outpath = '/Users/robin/git/graphviz/bolza-surface/bolza_r4.html'
    with open(outpath, 'w') as f:
        f.write(html)
    print(f"  Written to {outpath}")
    print("\nDone!")


# ══════════════════════════════════════════════════════════════════════
#  HTML/Three.js visualization
# ══════════════════════════════════════════════════════════════════════

def generate_html(r4_coords, triangles, hue_per_vertex, parent_tri, disk_pts):
    """Generate self-contained HTML with Three.js visualization.

    Features:
    - 4D rotation slider (rotating in the (Re(y), Im(y)) plane)
    - Color by arg(x)
    - Wireframe overlay for original (2,3,8) triangle boundaries
    - Toggle between projection modes
    """
    # Prepare data as JSON
    positions = r4_coords.tolist()
    tris = triangles.tolist()
    hues = hue_per_vertex.tolist()
    parents = [int(p) for p in parent_tri]
    disk = [(float(z.real), float(z.imag)) for z in disk_pts]

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Bolza Surface in R^4 — y² = x⁵ - x</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box }}
body {{ background:#0a0a12; overflow:hidden; font-family:'Courier New',monospace }}
#info {{
  position:absolute; top:12px; left:16px; color:#8899aa;
  font-size:13px; line-height:1.6; pointer-events:none; max-width:420px;
  text-shadow: 0 0 8px rgba(0,0,0,0.8);
}}
#info h2 {{ color:#ccd0d8; font-size:16px; margin-bottom:4px; font-weight:normal }}
#info .n {{ color:#e8b84d }}
#controls {{
  position:absolute; bottom:20px; left:50%; transform:translateX(-50%);
  background:rgba(20,20,30,0.85); padding:12px 24px; border-radius:8px;
  display:flex; gap:24px; align-items:center; color:#99aabb; font-size:13px;
  font-family:'Courier New',monospace;
  border: 1px solid rgba(100,120,140,0.3);
}}
#controls label {{ display:flex; align-items:center; gap:8px; cursor:pointer }}
#controls input[type=range] {{ width:160px; accent-color:#e8b84d }}
#controls button {{
  background:rgba(100,120,140,0.3); border:1px solid rgba(100,120,140,0.4);
  color:#aabbcc; padding:4px 12px; border-radius:4px; cursor:pointer;
  font-family:inherit; font-size:12px;
}}
#controls button:hover {{ background:rgba(100,120,140,0.5) }}
#controls button.active {{ background:rgba(232,184,77,0.3); border-color:#e8b84d; color:#e8b84d }}
</style></head><body>
<div id="info">
<h2>Bolza Surface — R⁴ Embedding</h2>
y² = x⁵ − x embedded in C² = R⁴<br>
<span class="n">{len(tris)}</span> triangles, <span class="n">{len(positions)}</span> vertices<br>
Color = arg(x) on the hyperelliptic projection<br>
4D rotation mixes Re(y) and Im(y)
</div>
<div id="controls">
  <label>4D angle: <input type="range" id="angle4d" min="0" max="314" value="0">
    <span id="angleVal">0°</span></label>
  <button id="btnWire" class="active" onclick="toggleWire()">wireframe</button>
  <button id="btnProj" onclick="cycleProj()">proj: drop Im(y)</button>
</div>
<script type="importmap">
{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}}}
</script>
<script type="module">
import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';

// --- Data ---
const R4 = {json.dumps(positions)};
const T = {json.dumps(tris)};
const H = {json.dumps(hues)};
const PAR = {json.dumps(parents)};

// --- Scene setup ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0a0a12);
const cam = new THREE.PerspectiveCamera(45, innerWidth/innerHeight, 0.01, 200);
cam.position.set(0, 5, 12);
const ren = new THREE.WebGLRenderer({{ antialias: true }});
ren.setSize(innerWidth, innerHeight);
ren.setPixelRatio(Math.min(devicePixelRatio, 2));
ren.toneMapping = THREE.ACESFilmicToneMapping;
ren.toneMappingExposure = 1.2;
document.body.appendChild(ren.domElement);
const ctl = new OrbitControls(cam, ren.domElement);
ctl.enableDamping = true; ctl.dampingFactor = 0.05;
ctl.autoRotate = true; ctl.autoRotateSpeed = 0.3;

// Lights
scene.add(new THREE.AmbientLight(0x556677, 0.6));
let l;
l = new THREE.DirectionalLight(0xffeedd, 1.2); l.position.set(5,8,5); scene.add(l);
l = new THREE.DirectionalLight(0x99aacc, 0.5); l.position.set(-4,2,-3); scene.add(l);
l = new THREE.DirectionalLight(0x8899bb, 0.3); l.position.set(0,-4,5); scene.add(l);

// --- Projection ---
let projMode = 0; // 0 = drop Im(y), 1 = drop Re(y), 2 = drop Im(x), 3 = drop Re(x)
const projNames = ['drop Im(y)', 'drop Re(y)', 'drop Im(x)', 'drop Re(x)'];

function project(r4, angle4d) {{
  // First apply 4D rotation in the (Re(y), Im(y)) plane
  const c = Math.cos(angle4d), s = Math.sin(angle4d);
  const x0 = r4[0], x1 = r4[1];
  const x2 = c * r4[2] - s * r4[3];
  const x3 = s * r4[2] + c * r4[3];
  // Then project to 3D
  switch(projMode) {{
    case 0: return [x0, x1, x2];       // drop Im(y)
    case 1: return [x0, x1, x3];       // drop Re(y)
    case 2: return [x0, x2, x3];       // drop Im(x)
    case 3: return [x1, x2, x3];       // drop Re(x)
  }}
}}

// --- HSL to RGB ---
function hsl2rgb(h, s, l) {{
  if (s === 0) return [l, l, l];
  const q = l < 0.5 ? l*(1+s) : l+s-l*s, p = 2*l - q;
  function f(t) {{
    if (t<0) t+=1; if (t>1) t-=1;
    if (t < 1/6) return p + (q-p)*6*t;
    if (t < 1/2) return q;
    if (t < 2/3) return p + (q-p)*(2/3-t)*6;
    return p;
  }}
  return [f(h+1/3), f(h), f(h-1/3)];
}}

// --- Build geometry ---
let meshGeo, wireMesh, surfMesh;
let showWire = true;

function buildMesh(angle4d) {{
  const ps = [], ns = [], cs = [];
  for (let i = 0; i < T.length; i++) {{
    const [a, b, c] = T[i];
    const pa = project(R4[a], angle4d);
    const pb = project(R4[b], angle4d);
    const pc = project(R4[c], angle4d);

    // Normal
    const u = [pb[0]-pa[0], pb[1]-pa[1], pb[2]-pa[2]];
    const v = [pc[0]-pa[0], pc[1]-pa[1], pc[2]-pa[2]];
    const nx = u[1]*v[2]-u[2]*v[1], ny = u[2]*v[0]-u[0]*v[2], nz = u[0]*v[1]-u[1]*v[0];
    const nl = Math.sqrt(nx*nx+ny*ny+nz*nz) || 1;

    for (const idx of [a, b, c]) {{
      const p = project(R4[idx], angle4d);
      ps.push(p[0], p[1], p[2]);
      ns.push(nx/nl, ny/nl, nz/nl);
      const [r, g, bl] = hsl2rgb(H[idx], 0.7, 0.55);
      cs.push(r, g, bl);
    }}
  }}

  if (meshGeo) meshGeo.dispose();
  meshGeo = new THREE.BufferGeometry();
  meshGeo.setAttribute('position', new THREE.Float32BufferAttribute(ps, 3));
  meshGeo.setAttribute('normal', new THREE.Float32BufferAttribute(ns, 3));
  meshGeo.setAttribute('color', new THREE.Float32BufferAttribute(cs, 3));
  return meshGeo;
}}

function buildWireframe(angle4d) {{
  // Only show edges of the ORIGINAL 48 triangles (boundaries between parent triangles)
  const wp = [];
  // Group subdivision triangles by parent, find boundary edges
  // Simple approach: draw all subdivision wireframe but thin
  for (let i = 0; i < T.length; i++) {{
    const [a, b, c] = T[i];
    const pa = project(R4[a], angle4d);
    const pb = project(R4[b], angle4d);
    const pc = project(R4[c], angle4d);
    wp.push(pa[0],pa[1],pa[2], pb[0],pb[1],pb[2]);
    wp.push(pb[0],pb[1],pb[2], pc[0],pc[1],pc[2]);
    wp.push(pc[0],pc[1],pc[2], pa[0],pa[1],pa[2]);
  }}
  const wg = new THREE.BufferGeometry();
  wg.setAttribute('position', new THREE.Float32BufferAttribute(wp, 3));
  return wg;
}}

let angle4d = 0;

function rebuild() {{
  if (surfMesh) {{ scene.remove(surfMesh); surfMesh.geometry.dispose(); }}
  if (wireMesh) {{ scene.remove(wireMesh); wireMesh.geometry.dispose(); }}

  const geo = buildMesh(angle4d);
  surfMesh = new THREE.Mesh(geo, new THREE.MeshPhongMaterial({{
    vertexColors: true, side: THREE.DoubleSide, shininess: 30,
    specular: new THREE.Color(0x445566)
  }}));
  scene.add(surfMesh);

  if (showWire) {{
    const wg = buildWireframe(angle4d);
    wireMesh = new THREE.LineSegments(wg, new THREE.LineBasicMaterial({{
      color: 0x223344, transparent: true, opacity: 0.15
    }}));
    scene.add(wireMesh);
  }}
}}

rebuild();

// --- Controls ---
const slider = document.getElementById('angle4d');
const angleLabel = document.getElementById('angleVal');
slider.addEventListener('input', () => {{
  angle4d = parseFloat(slider.value) / 100;
  angleLabel.textContent = Math.round(angle4d * 180 / Math.PI) + '°';
  rebuild();
}});

window.toggleWire = function() {{
  showWire = !showWire;
  document.getElementById('btnWire').classList.toggle('active');
  rebuild();
}};

window.cycleProj = function() {{
  projMode = (projMode + 1) % 4;
  document.getElementById('btnProj').textContent = 'proj: ' + projNames[projMode];
  rebuild();
}};

window.onresize = () => {{
  cam.aspect = innerWidth/innerHeight;
  cam.updateProjectionMatrix();
  ren.setSize(innerWidth, innerHeight);
}};

(function animate() {{
  requestAnimationFrame(animate);
  ctl.update();
  ren.render(scene, cam);
}})();
</script></body></html>"""


if __name__ == '__main__':
    main()
