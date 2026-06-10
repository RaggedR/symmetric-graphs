#!/usr/bin/env python3
"""
Attempt 9: R^4 embedding via discrete Hodge theory.

Instead of Poincare series (slow convergence), solve for the holomorphic
differentials directly on the mesh using the Laplacian with twisted
boundary conditions from the octagon side-pairing transformations.
"""

import numpy as np
import json
import sys
from collections import deque, defaultdict
from scipy import sparse
from scipy.sparse.linalg import eigsh, eigs
from bolza_group import GL23, TriangleGroup, BolzaTiling


# ══════════════════════════════════════════════════════════════════════
#  Step 1: Find the 4 side-pairing transformations
# ══════════════════════════════════════════════════════════════════════

def find_side_pairings(tiling, tg):
    """Find SL(2,R) matrices that pair opposite sides of the octagon.

    The octagon sides (counterclockwise): a1 b1 a1^-1 b1^-1 a2 b2 a2^-1 b2^-1
    Side-pairings: a1 <-> a1^-1, b1 <-> b1^-1, a2 <-> a2^-1, b2 <-> b2^-1

    For side a1 (V0->V1) paired with a1^-1 (V2<-V3, reversed):
    The transformation T maps V0->V3 and V1->V2.

    Returns dict: {'a1': M, 'b1': M, 'a2': M, 'b2': M} as SL(2,R) matrices.
    """
    ov = tiling.octagon_vertices()
    V = [tiling.vertices[i] for i in ov]  # 8 vertices in disk coords

    # Convert to upper half-plane for SL(2,R) computation
    Z = [tg.disk_to_halfplane(v) for v in V]

    # Side pairings (source side -> target side with reversed orientation):
    # a1: V0->V1 pairs with a1^-1: V3->V2 (reversed), so T maps V0->V3, V1->V2
    # b1: V1->V2 pairs with b1^-1: V4->V3 (reversed), so T maps V1->V4, V2->V3
    # a2: V4->V5 pairs with a2^-1: V7->V6 (reversed), so T maps V4->V7, V5->V6
    # b2: V5->V6 pairs with b2^-1: V0->V7 (reversed), so T maps V5->V0, V6->V7
    pairings = {
        'a1': (0, 1, 3, 2),  # V0->V3, V1->V2
        'b1': (1, 2, 4, 3),  # V1->V4, V2->V3
        'a2': (4, 5, 7, 6),  # V4->V7, V5->V6
        'b2': (5, 6, 0, 7),  # V5->V0, V6->V7
    }

    result = {}
    for name, (s1, s2, t1, t2) in pairings.items():
        # Find SL(2,R) matrix mapping Z[s1]->Z[t1] and Z[s2]->Z[t2]
        M = find_sl2r_mapping(Z[s1], Z[s2], Z[t1], Z[t2])
        result[name] = M

        # Verify
        err1 = abs(tg.mobius(M, Z[s1]) - Z[t1])
        err2 = abs(tg.mobius(M, Z[s2]) - Z[t2])
        det = abs(np.linalg.det(M) - 1)
        print(f"  {name}: V{s1}->V{t1} err={err1:.2e}, V{s2}->V{t2} err={err2:.2e}, det-1={det:.2e}")

    return result


def find_sl2r_mapping(z1, z2, w1, w2):
    """Find the SL(2,R) matrix mapping z1->w1 and z2->w2.

    Mobius: w = (az+b)/(cz+d), ad-bc=1.
    Two point conditions give 2 complex = 4 real equations for 3 real parameters
    (SL(2,R) has dimension 3). Overdetermined but consistent for isometries.
    """
    # Use the cross-ratio: the unique Mobius transformation mapping
    # z1->w1, z2->w2, and preserving the real axis (SL(2,R) condition).
    #
    # Approach: compose two Mobius transforms.
    # M1 maps z1->0, z2->1 (known formula)
    # M2 maps 0->w1, 1->w2 (known formula)
    # T = M2 @ M1

    # M1: z -> (z-z1)/(z-z2) * (z2-z1_bar)/(z1-z1_bar) ... no, for SL(2,R)
    # we need to map via the real axis.

    # Direct approach: solve the linear system.
    # az1 + b = w1(cz1 + d)  =>  az1 + b - w1*cz1 - w1*d = 0
    # az2 + b = w2(cz2 + d)  =>  az2 + b - w2*cz2 - w2*d = 0
    # ad - bc = 1
    #
    # Set d = 1 (can normalize later), solve for a, b, c:
    # a*z1 - w1*c*z1 + b - w1 = 0  =>  a*z1 + b + c*(-w1*z1) = w1
    # a*z2 - w2*c*z2 + b - w2 = 0  =>  a*z2 + b + c*(-w2*z2) = w2
    #
    # Two complex equations, 3 complex unknowns. But SL(2,R) means a,b,c,d are REAL.
    # So: 2 complex = 4 real equations for 3 real unknowns (a, b, c with d=1).
    # Overdetermined — solve least-squares.

    # Stack real and imaginary parts
    # Re: a*Re(z1) + b - c*Re(w1*z1) = Re(w1)
    # Im: a*Im(z1)     - c*Im(w1*z1) = Im(w1)
    # (b has no imaginary part since b is real and contributes only to Re)

    # Actually wait — a, b, c, d are all REAL for SL(2,R).
    # So from az + b = w(cz + d) with z, w complex and a,b,c,d real:
    # Re part: a*Re(z) + b = Re(w)*c*Re(z) - Im(w)*c*Im(z) + Re(w)*d
    #        = Re(w*(cz+d)) = c*Re(wz) + d*Re(w)
    # This is getting messy. Let me just use the real/imag decomposition.

    # From w = (az+b)/(cz+d) where a,b,c,d real:
    # w*(cz+d) = az+b
    # Expanding with z = x+iy, w = u+iv:
    # (u+iv)(c(x+iy)+d) = a(x+iy)+b
    # (u+iv)(cx+d + icy) = ax+b + iay
    # u(cx+d) - vcy + i[v(cx+d) + ucy] = ax+b + iay
    #
    # Real: u(cx+d) - vcy = ax + b
    # Imag: v(cx+d) + ucy = ay

    # For point 1: z1=x1+iy1, w1=u1+iv1:
    # u1(cx1+d) - v1*c*y1 = a*x1 + b        ... (R1)
    # v1(cx1+d) + u1*c*y1 = a*y1             ... (I1)

    # For point 2: z2=x2+iy2, w2=u2+iv2:
    # u2(cx2+d) - v2*c*y2 = a*x2 + b        ... (R2)
    # v2(cx2+d) + u2*c*y2 = a*y2             ... (I2)

    # Set d=1. Unknowns: a, b, c.
    # R1: a*x1 + b - c*(u1*x1 - v1*y1) = u1
    # I1: a*y1     - c*(-v1*x1 - u1*y1) = v1  =>  a*y1 + c*(v1*x1 + u1*y1) = v1
    # R2: a*x2 + b - c*(u2*x2 - v2*y2) = u2
    # I2: a*y2 + c*(v2*x2 + u2*y2) = v2

    x1, y1 = z1.real, z1.imag
    x2, y2 = z2.real, z2.imag
    u1, v1 = w1.real, w1.imag
    u2, v2 = w2.real, w2.imag

    # 4 equations, 3 unknowns [a, b, c]:
    A_mat = np.array([
        [x1, 1, -(u1*x1 - v1*y1)],
        [y1, 0,  (v1*x1 + u1*y1)],
        [x2, 1, -(u2*x2 - v2*y2)],
        [y2, 0,  (v2*x2 + u2*y2)],
    ], dtype=float)
    rhs = np.array([u1, v1, u2, v2], dtype=float)

    abc, residual, _, _ = np.linalg.lstsq(A_mat, rhs, rcond=None)
    a, b, c = abc
    d = 1.0

    # Normalize: ad - bc = 1
    det = a * d - b * c
    if det < 0:
        a, b, c, d = -a, -b, -c, -d
        det = -det
    scale = np.sqrt(abs(det))
    M = np.array([[a/scale, b/scale], [c/scale, d/scale]])
    return M


# ══════════════════════════════════════════════════════════════════════
#  Step 2: Build fine mesh with identified boundaries
# ══════════════════════════════════════════════════════════════════════

def build_fine_mesh(tiling, sub_n=12):
    """Subdivide the 48 triangles and track boundary vertices.

    Returns:
        verts: list of complex (disk coordinates)
        tris: list of (i,j,k) triangles
        boundary_pairs: list of (i, j, side_name) for identified boundary vertices
        octagon_sides: dict mapping side_name -> list of vertex indices on that side
    """
    verts = []
    tris = []
    vert_map = {}  # (round(re), round(im)) -> index
    tol = 1e-9

    def add_vert(w):
        key = (round(w.real / tol) * tol, round(w.imag / tol) * tol)
        if key not in vert_map:
            vert_map[key] = len(verts)
            verts.append(w)
        return vert_map[key]

    for ti, (wa, wb, wc) in enumerate(tiling.fund_triangles):
        idx_map = {}
        for i in range(sub_n + 1):
            for j in range(sub_n + 1 - i):
                k = sub_n - i - j
                w = (i * wa + j * wb + k * wc) / sub_n
                idx_map[(i, j)] = add_vert(w)

        for i in range(sub_n):
            for j in range(sub_n - i):
                v0 = idx_map[(i, j)]
                v1 = idx_map[(i + 1, j)]
                v2 = idx_map[(i, j + 1)]
                tris.append((v0, v1, v2))
                if i + j + 1 < sub_n:
                    v3 = idx_map[(i + 1, j + 1)]
                    tris.append((v1, v3, v2))

    return verts, tris


def classify_boundary(verts, tiling):
    """Classify vertices as interior or on specific octagon sides."""
    ov = tiling.octagon_vertices()
    V = [tiling.vertices[i] for i in ov]  # 8 octagon vertices

    R_oct = max(abs(v) for v in V)
    boundary_threshold = R_oct * 0.98  # vertices near octagon boundary

    # Octagon sides: side k goes from V[k] to V[(k+1)%8]
    side_names = ['a1', 'b1', 'a1_inv', 'b1_inv', 'a2', 'b2', 'a2_inv', 'b2_inv']

    # For each vertex near the boundary, determine which side it's on
    # by finding the nearest octagon edge
    vertex_side = {}  # vertex_index -> (side_name, parameter t in [0,1])

    for vi, w in enumerate(verts):
        if abs(w) < boundary_threshold:
            continue  # interior vertex

        # Find closest octagon side
        best_dist = float('inf')
        best_side = None
        best_t = 0

        for k in range(8):
            p0 = V[k]
            p1 = V[(k + 1) % 8]
            # Project w onto line segment p0-p1
            d = p1 - p0
            if abs(d) < 1e-12:
                continue
            t = ((w - p0) * np.conj(d)).real / (abs(d) ** 2)
            t = max(0, min(1, t))
            proj = p0 + t * d
            dist = abs(w - proj)

            if dist < best_dist:
                best_dist = dist
                best_side = k
                best_t = t

        if best_dist < 0.02:  # close enough to a side
            vertex_side[vi] = (side_names[best_side], best_t, best_side)

    return vertex_side


def find_boundary_pairs(verts, vertex_side, side_pairings, tg):
    """Match boundary vertices on paired sides.

    For side a1 (V0->V1) paired with a1_inv (V2<-V3):
    A vertex at parameter t on a1 matches the vertex at parameter 1-t on a1_inv.
    """
    # Group vertices by side
    side_verts = defaultdict(list)
    for vi, (sname, t, sidx) in vertex_side.items():
        side_verts[sname].append((vi, t))

    # Side pairing map (with reversed orientation)
    pairing_map = {
        'a1': 'a1_inv', 'a1_inv': 'a1',
        'b1': 'b1_inv', 'b1_inv': 'b1',
        'a2': 'a2_inv', 'a2_inv': 'a2',
        'b2': 'b2_inv', 'b2_inv': 'b2',
    }

    # For each side, sort vertices by parameter
    pairs = []  # (vi, vj, side_name, T_matrix)

    for src_side in ['a1', 'b1', 'a2', 'b2']:
        dst_side = pairing_map[src_side]
        src_list = sorted(side_verts.get(src_side, []), key=lambda x: x[1])
        dst_list = sorted(side_verts.get(dst_side, []), key=lambda x: x[1])

        if not src_list or not dst_list:
            print(f"  WARNING: no vertices on side {src_side} or {dst_side}")
            continue

        T = side_pairings[src_side]

        # Match by reversed parameter: src at t matches dst at 1-t
        for vi, t_src in src_list:
            t_target = 1.0 - t_src
            # Find closest dst vertex to t_target
            best_vj = None
            best_dt = float('inf')
            for vj, t_dst in dst_list:
                dt = abs(t_dst - t_target)
                if dt < best_dt:
                    best_dt = dt
                    best_vj = vj

            if best_vj is not None and best_dt < 0.05:
                pairs.append((vi, best_vj, src_side, T))

    return pairs


# ══════════════════════════════════════════════════════════════════════
#  Step 3: Cotangent Laplacian with twisted boundary conditions
# ══════════════════════════════════════════════════════════════════════

def build_cotangent_laplacian(verts, tris):
    """Build the cotangent-weight Laplacian matrix.

    L[i,j] = (cot alpha_ij + cot beta_ij) / 2 for neighbors i,j
    L[i,i] = -sum_j L[i,j]

    Returns sparse matrix of shape (N, N).
    """
    N = len(verts)
    W = defaultdict(float)  # edge weights

    for i0, i1, i2 in tris:
        z0 = verts[i0]
        z1 = verts[i1]
        z2 = verts[i2]

        # Edges as complex numbers
        e01 = z1 - z0
        e02 = z2 - z0
        e12 = z2 - z1

        # Cotangent of angle at each vertex
        # At v0: angle between e01 and e02
        cross01_02 = (e01.real * e02.imag - e01.imag * e02.real)
        dot01_02 = (e01.real * e02.real + e01.imag * e02.imag)
        cot0 = dot01_02 / cross01_02 if abs(cross01_02) > 1e-15 else 0

        e10 = -e01
        e12_ = e12
        cross10_12 = (e10.real * e12_.imag - e10.imag * e12_.real)
        dot10_12 = (e10.real * e12_.real + e10.imag * e12_.imag)
        cot1 = dot10_12 / cross10_12 if abs(cross10_12) > 1e-15 else 0

        e20 = -e02
        e21 = -e12
        cross20_21 = (e20.real * e21.imag - e20.imag * e21.real)
        dot20_21 = (e20.real * e21.real + e20.imag * e21.imag)
        cot2 = dot20_21 / cross20_21 if abs(cross20_21) > 1e-15 else 0

        # Cotangent weight for edge (i,j) = cot of angle at opposite vertex / 2
        W[(i1, i2)] += cot0 / 2
        W[(i2, i1)] += cot0 / 2
        W[(i0, i2)] += cot1 / 2
        W[(i2, i0)] += cot1 / 2
        W[(i0, i1)] += cot2 / 2
        W[(i1, i0)] += cot2 / 2

    rows, cols, vals = [], [], []
    for (i, j), w in W.items():
        rows.append(i)
        cols.append(j)
        vals.append(w)

    L = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))
    # Diagonal: L[i,i] = -sum of off-diagonal entries in row i
    diag = -np.array(L.sum(axis=1)).flatten()
    L = L + sparse.diags(diag)
    return L


def solve_holomorphic_differentials(verts, tris, boundary_pairs, tg):
    """Solve for the 2D space of holomorphic differentials.

    Build the Laplacian with twisted boundary conditions and find its kernel.
    """
    N = len(verts)
    print(f"  Building cotangent Laplacian ({N} vertices, {len(tris)} triangles)...")
    L = build_cotangent_laplacian(verts, tris)

    # Identify which vertices are "dependent" (determined by a partner via twist)
    # and which are "independent"
    dependent = {}  # vj -> (vi, twist_factor)
    for vi, vj, side_name, T in boundary_pairs:
        # phi(vj) = (c*z_i + d)^2 * phi(vi)  where T = [[a,b],[c,d]]
        # z_i is the H-coordinate of vertex vi
        z_i = tg.disk_to_halfplane(verts[vi])
        c, d = T[1, 0], T[1, 1]
        twist = (c * z_i + d) ** 2
        if vj not in dependent:
            dependent[vj] = (vi, twist)

    independent = [i for i in range(N) if i not in dependent]
    N_ind = len(independent)
    ind_map = {v: k for k, v in enumerate(independent)}

    print(f"  {N_ind} independent vertices, {len(dependent)} dependent")

    # Build the reduced system: for each independent vertex i,
    # the Laplacian equation is:
    # sum_{j neighbor of i} w_ij * (phi_j - phi_i) = 0
    #
    # If j is independent: phi_j contributes directly
    # If j is dependent on k: phi_j = twist_jk * phi_k

    rows, cols, vals = [], [], []
    for i in independent:
        k_i = ind_map[i]
        # Row k_i: Laplacian at vertex i
        # Diagonal contribution
        diag_val = L[i, i]

        # Off-diagonal contributions
        row_start = L.indptr[i]
        row_end = L.indptr[i + 1]
        for ptr in range(row_start, row_end):
            j = L.indices[ptr]
            w_ij = L.data[ptr]
            if j == i:
                continue  # handle diagonal separately

            if j in ind_map:
                # j is independent
                k_j = ind_map[j]
                rows.append(k_i)
                cols.append(k_j)
                vals.append(complex(w_ij))
            elif j in dependent:
                # j is dependent on some independent vertex
                v_ind, twist = dependent[j]
                if v_ind in ind_map:
                    k_j = ind_map[v_ind]
                    rows.append(k_i)
                    cols.append(k_j)
                    vals.append(complex(w_ij) * twist)
                # else: dependent on another dependent — shouldn't happen

        # Diagonal
        rows.append(k_i)
        cols.append(k_i)
        vals.append(complex(diag_val))

    L_reduced = sparse.csr_matrix(
        (vals, (rows, cols)), shape=(N_ind, N_ind), dtype=complex)

    print(f"  Reduced system: {N_ind}x{N_ind} complex sparse matrix")
    print(f"  Non-zeros: {L_reduced.nnz}")

    # Find the kernel (smallest eigenvalues)
    print("  Finding smallest eigenvalues...")
    try:
        eigenvalues, eigenvectors = eigs(L_reduced, k=6, sigma=0, which='LM')
    except Exception as e:
        print(f"  eigs failed ({e}), trying eigsh on Hermitian part...")
        # Make it Hermitian: (L + L^H) / 2
        L_herm = (L_reduced + L_reduced.conj().T) / 2
        eigenvalues, eigenvectors = eigsh(L_herm, k=6, which='SM')

    # Sort by magnitude
    order = np.argsort(np.abs(eigenvalues))
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    print("  Smallest eigenvalues:")
    for i, ev in enumerate(eigenvalues):
        print(f"    lambda_{i} = {ev:.6e}")

    # The kernel should be ~2-dimensional (the holomorphic differentials)
    # Take the two eigenvectors with smallest eigenvalue
    phi1_reduced = eigenvectors[:, 0]
    phi2_reduced = eigenvectors[:, 1]

    # Expand back to full vertex set
    phi1_full = np.zeros(N, dtype=complex)
    phi2_full = np.zeros(N, dtype=complex)
    for i in independent:
        k = ind_map[i]
        phi1_full[i] = phi1_reduced[k]
        phi2_full[i] = phi2_reduced[k]
    for j, (v_ind, twist) in dependent.items():
        phi1_full[j] = twist * phi1_full[v_ind]
        phi2_full[j] = twist * phi2_full[v_ind]

    return phi1_full, phi2_full, eigenvalues


# ══════════════════════════════════════════════════════════════════════
#  Step 4: Extract x, compute y, embed in R^4
# ══════════════════════════════════════════════════════════════════════

def calibrate_and_embed(verts, phi1, phi2, tiling, tg):
    """Compute x = phi2/phi1 (up to Mobius), then y, then R^4."""

    ratio = np.where(np.abs(phi1) > 1e-15, phi2 / phi1, np.inf + 0j)

    # Find Weierstrass points and calibrate
    wp_disk = tiling.weierstrass_points_disk()
    wp_indices = []
    for w in wp_disk:
        dists = [abs(verts[i] - w) for i in range(len(verts))]
        wp_indices.append(np.argmin(dists))

    print("\n  Ratio phi2/phi1 at Weierstrass points:")
    ratios_at_wp = []
    for i, wi in enumerate(wp_indices):
        r = ratio[wi]
        print(f"    W{i} (vi={wi}): ratio = {r:.6f}, |phi1| = {abs(phi1[wi]):.6e}")
        ratios_at_wp.append(r)

    ratios_at_wp = np.array(ratios_at_wp)

    # Try all assignments and find best Mobius fit
    from itertools import permutations

    # Identify c-fixed Weierstrass points (closest and farthest from center)
    wp_dists = [(i, abs(wp_disk[i])) for i in range(6)]
    wp_dists.sort(key=lambda x: x[1])
    center_idx = wp_dists[0][0]
    far_idx = wp_dists[-1][0]
    non_fixed = [i for i in range(6) if i != center_idx and i != far_idx]

    best_error = np.inf
    best_mobius = None
    best_assignment = None

    for inf_choice in [center_idx, far_idx]:
        zero_choice = far_idx if inf_choice == center_idx else center_idx
        for perm in permutations([1, -1, 1j, -1j]):
            x_assign = {}
            x_assign[zero_choice] = 0
            x_assign[inf_choice] = np.inf
            for k, idx in enumerate(non_fixed):
                x_assign[idx] = perm[k]

            finite_r = []
            finite_x = []
            for idx in range(6):
                if np.isfinite(ratios_at_wp[idx]) and np.isfinite(x_assign[idx]):
                    finite_r.append(ratios_at_wp[idx])
                    finite_x.append(x_assign[idx])

            if len(finite_r) < 3:
                continue

            try:
                A, B, C, D = find_mobius_transform(finite_r[:3], finite_x[:3])
            except:
                continue

            error = sum(abs(apply_mobius(ratios_at_wp[idx], A, B, C, D) - x_assign[idx])
                        for idx in range(6)
                        if np.isfinite(ratios_at_wp[idx]) and np.isfinite(x_assign[idx]))

            if error < best_error:
                best_error = error
                best_mobius = (A, B, C, D)
                best_assignment = dict(x_assign)

    if best_mobius is None:
        print("  ERROR: Mobius calibration failed!")
        return None, None, None

    A, B, C, D = best_mobius
    print(f"\n  Best Mobius calibration (error = {best_error:.6e}):")
    for idx in range(6):
        print(f"    W{idx} -> x = {best_assignment[idx]}")

    # Apply Mobius to get x
    x_vals = apply_mobius(ratio, A, B, C, D)

    # y = sqrt(x^5 - x)
    f_x = x_vals**5 - x_vals
    y_vals = np.sqrt(f_x)

    # R^4 embedding with stereographic projection on x to avoid blowup
    # Map x to S^2: (2Re(x)/(1+|x|^2), 2Im(x)/(1+|x|^2), (|x|^2-1)/(|x|^2+1))
    x_abs2 = np.abs(x_vals)**2
    denom = 1 + x_abs2
    s0 = 2 * x_vals.real / denom  # stereographic x
    s1 = 2 * x_vals.imag / denom  # stereographic y
    s2 = (x_abs2 - 1) / denom     # stereographic z

    # For the y-coordinate, also apply bounded normalization
    y_bounded = np.arcsinh(y_vals / 2.0) * 2.0

    # R^4 coords: (s0, s1, s2, Re(y_bounded))  — 3 from stereographic x + 1 from y
    r4 = np.column_stack([s0, s1, s2, y_bounded.real])
    r4[~np.isfinite(r4)] = 0.0

    return x_vals, y_vals, r4


def find_mobius_transform(r_vals, x_vals):
    """Find Mobius x = (Ar+B)/(Cr+D) mapping r_vals to x_vals. See attempt 8."""
    r1, r2, r3 = r_vals[:3]
    x1, x2, x3 = x_vals[:3]
    M = np.array([
        [-r1, -1, r1*x1],
        [-r2, -1, r2*x2],
        [-r3, -1, r3*x3],
    ], dtype=complex)
    rhs = np.array([-x1, -x2, -x3], dtype=complex)
    ABC = np.linalg.solve(M, rhs)
    A, B, C = ABC
    D = 1.0 + 0j
    det = A*D - B*C
    s = np.sqrt(det)
    return A/s, B/s, C/s, D/s


def apply_mobius(r, A, B, C, D):
    return (A * r + B) / (C * r + D)


# ══════════════════════════════════════════════════════════════════════
#  Visualization
# ══════════════════════════════════════════════════════════════════════

def generate_html(r4, tris, x_vals, verts_disk):
    """Three.js visualization with stereographic embedding."""
    positions = r4.tolist()
    tri_list = [[int(a), int(b), int(c)] for a, b, c in tris]
    hues = (np.angle(x_vals) / (2 * np.pi) % 1.0).tolist()

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Bolza Surface — Attempt 9 (Discrete Hodge)</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{background:#0a0a12;overflow:hidden;font-family:'Courier New',monospace}}
#info{{position:absolute;top:12px;left:16px;color:#8899aa;font-size:13px;line-height:1.6;pointer-events:none;max-width:420px;text-shadow:0 0 8px rgba(0,0,0,0.8)}}
#info h2{{color:#ccd0d8;font-size:16px;margin-bottom:4px;font-weight:normal}}
#info .n{{color:#e8b84d}}
</style></head><body>
<div id="info">
<h2>Bolza Surface — Discrete Hodge (Attempt 9)</h2>
Holomorphic differentials from Laplacian kernel<br>
x = phi2/phi1 on y² = x⁵ − x<br>
<span class="n">{len(tri_list)}</span> triangles, <span class="n">{len(positions)}</span> vertices<br>
Stereographic projection on x ∈ CP¹ → S²
</div>
<script type="importmap">
{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}}}
</script>
<script type="module">
import*as THREE from'three';
import{{OrbitControls}}from'three/addons/controls/OrbitControls.js';

const R4={json.dumps(positions)};
const T={json.dumps(tri_list)};
const H={json.dumps(hues)};

const scene=new THREE.Scene();scene.background=new THREE.Color(0x0a0a12);
const cam=new THREE.PerspectiveCamera(45,innerWidth/innerHeight,0.01,200);
cam.position.set(0,3,6);
const ren=new THREE.WebGLRenderer({{antialias:true}});
ren.setSize(innerWidth,innerHeight);ren.setPixelRatio(Math.min(devicePixelRatio,2));
document.body.appendChild(ren.domElement);
const ctl=new OrbitControls(cam,ren.domElement);
ctl.enableDamping=true;ctl.dampingFactor=0.05;ctl.autoRotate=true;ctl.autoRotateSpeed=0.3;

scene.add(new THREE.AmbientLight(0x556677,0.6));
let l;
l=new THREE.DirectionalLight(0xffeedd,1.2);l.position.set(5,8,5);scene.add(l);
l=new THREE.DirectionalLight(0x99aacc,0.5);l.position.set(-4,2,-3);scene.add(l);

function hsl2rgb(h,s,l){{if(s===0)return[l,l,l];const q=l<0.5?l*(1+s):l+s-l*s,p=2*l-q;function f(t){{if(t<0)t+=1;if(t>1)t-=1;if(t<1/6)return p+(q-p)*6*t;if(t<1/2)return q;if(t<2/3)return p+(q-p)*(2/3-t)*6;return p}};return[f(h+1/3),f(h),f(h-1/3)]}}

const ps=[],ns=[],cs=[];
for(let i=0;i<T.length;i++){{
  const[a,b,c]=T[i];
  const pa=R4[a],pb=R4[b],pc=R4[c];
  const u=[pb[0]-pa[0],pb[1]-pa[1],pb[2]-pa[2]];
  const v=[pc[0]-pa[0],pc[1]-pa[1],pc[2]-pa[2]];
  const nx=u[1]*v[2]-u[2]*v[1],ny=u[2]*v[0]-u[0]*v[2],nz=u[0]*v[1]-u[1]*v[0];
  const nl=Math.sqrt(nx*nx+ny*ny+nz*nz)||1;
  for(const idx of[a,b,c]){{
    // Use first 3 of 4 R^4 coords as position
    ps.push(R4[idx][0],R4[idx][1],R4[idx][2]);
    ns.push(nx/nl,ny/nl,nz/nl);
    const[r,g,bl]=hsl2rgb(H[idx],0.7,0.55);
    cs.push(r,g,bl);
  }}
}}

const geo=new THREE.BufferGeometry();
geo.setAttribute('position',new THREE.Float32BufferAttribute(ps,3));
geo.setAttribute('normal',new THREE.Float32BufferAttribute(ns,3));
geo.setAttribute('color',new THREE.Float32BufferAttribute(cs,3));
scene.add(new THREE.Mesh(geo,new THREE.MeshPhongMaterial({{
  vertexColors:true,side:THREE.DoubleSide,shininess:30,specular:new THREE.Color(0x445566)
}})));

// Wireframe
const wp=[];
for(const[a,b,c]of T){{
  wp.push(R4[a][0],R4[a][1],R4[a][2],R4[b][0],R4[b][1],R4[b][2]);
  wp.push(R4[b][0],R4[b][1],R4[b][2],R4[c][0],R4[c][1],R4[c][2]);
  wp.push(R4[c][0],R4[c][1],R4[c][2],R4[a][0],R4[a][1],R4[a][2]);
}}
const wg=new THREE.BufferGeometry();
wg.setAttribute('position',new THREE.Float32BufferAttribute(wp,3));
scene.add(new THREE.LineSegments(wg,new THREE.LineBasicMaterial({{color:0x223344,transparent:true,opacity:0.12}})));

window.onresize=()=>{{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();ren.setSize(innerWidth,innerHeight)}};
(function a(){{requestAnimationFrame(a);ctl.update();ren.render(scene,cam)}})();
</script></body></html>"""


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  Bolza Surface — Attempt 9: Discrete Hodge Theory")
    print("=" * 60)

    # Build infrastructure
    print("\n[1/5] Building group infrastructure...")
    gl23 = GL23()
    tg = TriangleGroup()
    tiling = BolzaTiling(gl23, tg)
    print(f"  48 triangles, {len(tiling.vertices)} vertices")

    # Find side-pairing transformations
    print("\n[2/5] Finding side-pairing transformations...")
    side_pairings = find_side_pairings(tiling, tg)

    # Build fine mesh
    print("\n[3/5] Building fine mesh...")
    SUB_N = 12
    verts, tris = build_fine_mesh(tiling, sub_n=SUB_N)
    print(f"  {len(verts)} vertices, {len(tris)} triangles")

    # Classify boundary vertices
    vertex_side = classify_boundary(verts, tiling)
    print(f"  {len(vertex_side)} boundary vertices classified")
    side_counts = defaultdict(int)
    for _, (sname, _, _) in vertex_side.items():
        side_counts[sname] += 1
    for sname in sorted(side_counts):
        print(f"    {sname}: {side_counts[sname]} vertices")

    # Find boundary pairs
    boundary_pairs = find_boundary_pairs(verts, vertex_side, side_pairings, tg)
    print(f"  {len(boundary_pairs)} boundary pairs matched")

    # Solve for holomorphic differentials
    print("\n[4/5] Solving for holomorphic differentials...")
    phi1, phi2, eigenvalues = solve_holomorphic_differentials(
        verts, tris, boundary_pairs, tg)

    # Calibrate and embed
    print("\n[5/5] Computing R^4 embedding...")
    x_vals, y_vals, r4 = calibrate_and_embed(verts, phi1, phi2, tiling, tg)

    if r4 is None:
        print("FAILED")
        return

    # Validation
    print("\n  Validation:")
    finite_mask = np.isfinite(x_vals) & np.isfinite(y_vals)
    if finite_mask.any():
        residual = np.abs(y_vals[finite_mask]**2 - (x_vals[finite_mask]**5 - x_vals[finite_mask]))
        print(f"    |y^2 - (x^5-x)|: median={np.median(residual):.6e}, max={residual.max():.6e}")
    print(f"    Non-finite vertices: {np.sum(~np.isfinite(r4).all(axis=1))} / {len(r4)}")
    print(f"    R^4 ranges: {r4.min(axis=0)} to {r4.max(axis=0)}")

    # Generate HTML
    outpath = '/Users/robin/git/graphviz/bolza-surface/bolza_r4_v2.html'
    html = generate_html(r4, tris, x_vals, verts)
    with open(outpath, 'w') as f:
        f.write(html)
    print(f"\n  Written to {outpath}")
    print("Done!")


if __name__ == '__main__':
    main()
