#!/usr/bin/env python3
"""
Attempt 10: R^4 embedding using GL(2,3) labels as single source of truth.

Key insight: every triangle in the BFS carries a GL(2,3) label AND an SL(2,R)
matrix. Two triangles with the same label differ by a Gamma element. Side
pairings, Weierstrass points, and boundary identifications all follow from
this data — no geometric guessing.
"""

import numpy as np
import json
from collections import defaultdict, deque
from itertools import product
from scipy import sparse
from scipy.sparse.linalg import eigs


# ══════════════════════════════════════════════════════════════════════
#  GL(2,3) — inline to avoid import issues
# ══════════════════════════════════════════════════════════════════════

def build_gl23():
    mm3 = lambda A, B: (A @ B) % 3
    mk3 = lambda A: tuple((A % 3).flatten())
    md3 = lambda A: (A[0,0]*A[1,1] - A[0,1]*A[1,0]) % 3
    def mi3(A):
        d = md3(A); di = 1 if d == 1 else 2
        return (di * np.array([[A[1,1],(-A[0,1])%3],[(-A[1,0])%3,A[0,0]]])) % 3

    el = {}
    for a,b,c,d in product(range(3), repeat=4):
        A = np.array([[a,b],[c,d]])
        if md3(A) == 0: continue
        k = mk3(A)
        if k not in el: el[k] = A % 3

    mats = list(el.values())
    keys = list(el.keys())
    ki = {k: i for i, k in enumerate(keys)}

    mult = lambda i, j: ki[mk3(mm3(mats[i], mats[j]))]
    inv = lambda i: ki[mk3(mi3(mats[i]))]
    I = ki[mk3(np.eye(2, dtype=int))]

    orders = []
    for i in range(48):
        c = i
        for n in range(1, 49):
            if c == I: orders.append(n); break
            c = mult(c, i)

    # Find (2,3,8) generators
    ag = bg = cg = None
    for a in range(48):
        if orders[a] != 2: continue
        for b in range(48):
            if orders[b] != 3: continue
            if orders[mult(a, b)] == 8:
                ag, bg, cg = a, b, mult(a, b); break
        if ag: break

    bg_inv = inv(bg)
    return mats, mult, inv, I, orders, ag, bg, cg, bg_inv


# ══════════════════════════════════════════════════════════════════════
#  SL(2,R) generators
# ══════════════════════════════════════════════════════════════════════

def build_sl2r():
    a_real = np.array([[0, 1], [-1, 0]], dtype=float)
    cp8 = np.cos(np.pi/8); sp3 = np.sin(np.pi/3); cp3 = np.cos(np.pi/3)
    rhs = 2*cp8/sp3; fp = (rhs + np.sqrt(rhs**2 - 4))/2
    b_real = np.array([[cp3, -fp*sp3], [sp3/fp, cp3]], dtype=float)
    b_inv = np.linalg.inv(b_real)
    c_real = a_real @ b_real

    def fixed_pt(M):
        a,b,c,d = M[0,0],M[0,1],M[1,0],M[1,1]
        disc = 4-(a+d)**2
        z1 = ((a-d)+1j*np.sqrt(max(disc,0)))/(2*c)
        z2 = ((a-d)-1j*np.sqrt(max(disc,0)))/(2*c)
        return z1 if z1.imag > 0 else z2

    z_a, z_b, z_c = fixed_pt(a_real), fixed_pt(b_real), fixed_pt(c_real)

    def to_disk(z): return (z-1j)/(z+1j)
    def to_H(w): return 1j*(1+w)/(1-w)
    def mobius(M, z): return (M[0,0]*z+M[0,1])/(M[1,0]*z+M[1,1])

    w_c_raw = to_disk(z_c)
    def recenter(w): return (w - w_c_raw)/(1 - np.conj(w_c_raw)*w)
    def recenter_inv(w): return (w + w_c_raw)/(1 + np.conj(w_c_raw)*w)
    def disk_to_H(w): return to_H(recenter_inv(w))

    return (a_real, b_real, b_inv, c_real,
            z_a, z_b, z_c, to_disk, to_H, mobius,
            recenter, recenter_inv, disk_to_H)


# ══════════════════════════════════════════════════════════════════════
#  BFS with full data retention
# ══════════════════════════════════════════════════════════════════════

def full_bfs(sl2r, gl23, max_tri=5000):
    (a_real, b_real, b_inv, c_real,
     z_a, z_b, z_c, to_disk, to_H, mobius,
     recenter, recenter_inv, disk_to_H) = sl2r
    mats, mult, inv, I, orders, ag, bg, cg, bg_inv = gl23

    gen_sl = [('a', a_real), ('b', b_real), ('B', b_inv)]
    gen_gl = {'a': ag, 'b': bg, 'B': bg_inv}

    def mat_key(M, tol=1e-6):
        return tuple(round(x/tol)*tol for x in M.flatten())

    all_data = []  # (wa, wb, wc, M_sl, g_gl)
    visited = set()
    # Track triangles by disk position to deduplicate ±M (same geometry)
    geo_seen = set()
    queue = deque([(np.eye(2), I)])

    while queue and len(all_data) < max_tri:
        M, g = queue.popleft()
        mk = mat_key(M)
        if mk in visited: continue
        visited.add(mk)

        try:
            wa = recenter(to_disk(mobius(M, z_a)))
            wb = recenter(to_disk(mobius(M, z_b)))
            wc = recenter(to_disk(mobius(M, z_c)))
        except: continue

        if max(abs(wa), abs(wb), abs(wc)) > 0.997: continue
        area = abs((wb-wa).real*(wc-wa).imag - (wb-wa).imag*(wc-wa).real)
        if area < 1e-10: continue

        # Deduplicate: M and -M give the same triangle geometry
        geo_key = (round(wa.real,5), round(wa.imag,5),
                   round(wb.real,5), round(wb.imag,5))
        if geo_key in geo_seen: continue
        geo_seen.add(geo_key)

        all_data.append((wa, wb, wc, M.copy(), g))

        for gn, gmat in gen_sl:
            queue.append((M @ gmat, mult(g, gen_gl[gn])))

    return all_data


# ══════════════════════════════════════════════════════════════════════
#  Proper fundamental domain from GL(2,3) labels
# ══════════════════════════════════════════════════════════════════════

def build_fund_domain(all_data):
    """Select one triangle per GL(2,3) label (closest centroid to origin).

    The centroid-based selection gives the Dirichlet domain by construction.
    """
    by_label = defaultdict(list)
    for i, (wa, wb, wc, M, g) in enumerate(all_data):
        by_label[g].append(i)

    fund = {}  # gl_label -> index in all_data
    for g, indices in by_label.items():
        # Pick the triangle whose centroid is closest to origin (Dirichlet domain)
        best_idx = min(indices, key=lambda i: abs(
            (all_data[i][0] + all_data[i][1] + all_data[i][2]) / 3))
        fund[g] = best_idx

    assert len(fund) == 48, f"Expected 48 labels, got {len(fund)}"
    return fund


# ══════════════════════════════════════════════════════════════════════
#  Side pairings from boundary adjacency
# ══════════════════════════════════════════════════════════════════════

def find_side_pairings(all_data, fund, sl2r, gl23):
    """Find side-pairing Gamma elements ALGEBRAICALLY from GL(2,3) labels.

    For each fund domain triangle (label g, matrix M), its three neighbors
    in the tiling have labels g*a, g*b, g*B. If the neighbor's fund domain
    representative has a DIFFERENT SL(2,R) matrix than the expected M*gen,
    they're in different Gamma copies, and the ratio gives a Gamma generator.
    """
    (a_real, b_real, b_inv, c_real,
     z_a, z_b, z_c, to_disk, to_H, mobius,
     recenter, recenter_inv, disk_to_H) = sl2r
    mats, mult, inv, I, orders, ag, bg, cg, bg_inv = gl23

    gen_pairs = [('a', ag, a_real), ('b', bg, b_real), ('B', bg_inv, b_inv)]

    seen_gammas = {}

    for g_label, fi in fund.items():
        wa, wb, wc, M_fund, g = all_data[fi]

        for gen_name, gen_gl, gen_sl in gen_pairs:
            # The neighbor via this generator has GL(2,3) label:
            g_nbr = mult(g_label, gen_gl)

            # Its SL(2,R) matrix IF it were in the same Gamma copy:
            M_expected = M_fund @ gen_sl

            # The actual fund domain representative for label g_nbr:
            fi_nbr = fund[g_nbr]
            M_actual = all_data[fi_nbr][3]

            # The Gamma element relating them:
            gamma = M_actual @ np.linalg.inv(M_expected)
            tr = np.trace(gamma)
            det = np.linalg.det(gamma)

            # Skip trivial (±I): |tr| ≈ ±2
            if abs(abs(tr) - 2) < 0.01:
                continue

            # De-duplicate
            gkey = tuple(round(x, 2) for x in gamma.flatten())
            if gkey not in seen_gammas:
                seen_gammas[gkey] = gamma

    gammas = list(seen_gammas.values())
    gammas = [G for G in gammas if abs(np.linalg.det(G) - 1) < 0.01]
    print(f"  Found {len(gammas)} distinct side-pairing elements")
    for i, G in enumerate(gammas):
        det = np.linalg.det(G)
        tr = np.trace(G)
        print(f"    gamma_{i}: tr={tr:.4f}, det={det:.6f}")

    return gammas


# ══════════════════════════════════════════════════════════════════════
#  Weierstrass points from hyperelliptic involution
# ══════════════════════════════════════════════════════════════════════

def find_weierstrass_points(gammas, sl2r):
    """Find the 6 Weierstrass points as fixed points of gamma · c^4."""
    (a_real, b_real, b_inv, c_real,
     z_a, z_b, z_c, to_disk, to_H, mobius,
     recenter, recenter_inv, disk_to_H) = sl2r

    c4 = np.linalg.matrix_power(c_real, 4)  # hyperelliptic involution in SL(2,R)

    def fixed_pt_sl2r(M):
        a,b,c,d = M[0,0],M[0,1],M[1,0],M[1,1]
        disc = (a+d)**2 - 4
        if disc >= 0 and abs(c) < 1e-12: return None
        if abs(c) < 1e-12: return None
        sq = np.sqrt(complex(disc))
        z1 = ((a-d) + sq)/(2*c)
        z2 = ((a-d) - sq)/(2*c)
        # Return the one in the upper half-plane
        for z in [z1, z2]:
            if z.imag > 1e-10: return z
        return None

    # Identity gives the center
    wp_H = []
    z0 = fixed_pt_sl2r(c4)
    if z0: wp_H.append(z0)

    # Each gamma and combinations give more Weierstrass points
    candidates = [np.eye(2)]
    for g in gammas:
        candidates.append(g)
        candidates.append(np.linalg.inv(g))
    for i, g1 in enumerate(gammas):
        for g2 in gammas[i:]:
            candidates.append(g1 @ g2)
            candidates.append(g2 @ g1)
            candidates.append(np.linalg.inv(g1) @ g2)

    for g in candidates:
        M = g @ c4
        z = fixed_pt_sl2r(M)
        if z is None: continue
        # Check if z is in or near the fundamental domain
        w = recenter(to_disk(z))
        if abs(w) > 0.95: continue  # too far out
        # Check if it's a new point (not Gamma-equivalent to existing)
        is_new = True
        for zp in wp_H:
            if abs(z - zp) < 1e-4:
                is_new = False; break
        if is_new:
            wp_H.append(z)

    # Convert to disk coordinates
    wp_disk = [recenter(to_disk(z)) for z in wp_H]
    return wp_H, wp_disk


# ══════════════════════════════════════════════════════════════════════
#  Mesh, boundary pairs, Laplacian, solve
# ══════════════════════════════════════════════════════════════════════

def subdivide_fund_domain(all_data, fund, sub_n=10):
    """Subdivide fund domain triangles, returning mesh with vertex tracking."""
    fund_indices = sorted(fund.values())
    verts = []
    tris = []
    vert_map = {}
    tol = 1e-9

    def add_v(w):
        k = (round(w.real/tol)*tol, round(w.imag/tol)*tol)
        if k not in vert_map:
            vert_map[k] = len(verts)
            verts.append(w)
        return vert_map[k]

    for fi in fund_indices:
        wa, wb, wc, _, _ = all_data[fi]
        idx = {}
        for i in range(sub_n+1):
            for j in range(sub_n+1-i):
                k = sub_n - i - j
                w = (i*wa + j*wb + k*wc)/sub_n
                idx[(i,j)] = add_v(w)
        for i in range(sub_n):
            for j in range(sub_n-i):
                tris.append((idx[(i,j)], idx[(i+1,j)], idx[(i,j+1)]))
                if i+j+1 < sub_n:
                    tris.append((idx[(i+1,j)], idx[(i+1,j+1)], idx[(i,j+1)]))

    return verts, tris


def find_boundary_pairs_from_gammas(verts, gammas, sl2r, tol=0.01):
    """For each vertex, check if any gamma maps it near another vertex."""
    (a_real, b_real, b_inv, c_real,
     z_a, z_b, z_c, to_disk, to_H, mobius,
     recenter, recenter_inv, disk_to_H) = sl2r

    # Build spatial index
    from scipy.spatial import KDTree
    coords = np.array([(v.real, v.imag) for v in verts])
    tree = KDTree(coords)

    pairs = []  # (vi, vj, gamma_matrix)
    paired = set()

    for G in gammas:
        # For each vertex, apply gamma and find nearest mesh vertex
        for vi, w in enumerate(verts):
            if vi in paired: continue
            z = disk_to_H(w)
            gz = mobius(G, z)
            gw = recenter(to_disk(gz))
            if abs(gw) > 0.99: continue  # outside disk

            dist, vj = tree.query([gw.real, gw.imag])
            if dist < tol and vj != vi and vj not in paired:
                pairs.append((vi, vj, G))
                paired.add(vi)
                paired.add(vj)

    return pairs


def build_laplacian(verts, tris):
    """Cotangent-weight Laplacian."""
    N = len(verts)
    W = defaultdict(float)
    for i0, i1, i2 in tris:
        z = [verts[i0], verts[i1], verts[i2]]
        for a in range(3):
            b, c = (a+1)%3, (a+2)%3
            eb = z[b] - z[a]
            ec = z[c] - z[a]
            cross = eb.real*ec.imag - eb.imag*ec.real
            dot = eb.real*ec.real + eb.imag*ec.imag
            cot = dot/cross if abs(cross) > 1e-15 else 0
            W[(b if b < 3 else [i0,i1,i2][b], c if c < 3 else [i0,i1,i2][c])] += cot/2

    # Redo with actual indices
    W2 = defaultdict(float)
    for i0, i1, i2 in tris:
        z0, z1, z2 = verts[i0], verts[i1], verts[i2]
        e01, e02 = z1-z0, z2-z0
        e10, e12 = z0-z1, z2-z1
        e20, e21 = z0-z2, z1-z2

        cr0 = e01.real*e02.imag - e01.imag*e02.real
        ct0 = (e01.real*e02.real + e01.imag*e02.imag)/(cr0 if abs(cr0)>1e-15 else 1e-15)
        cr1 = e10.real*e12.imag - e10.imag*e12.real
        ct1 = (e10.real*e12.real + e10.imag*e12.imag)/(cr1 if abs(cr1)>1e-15 else 1e-15)
        cr2 = e20.real*e21.imag - e20.imag*e21.real
        ct2 = (e20.real*e21.real + e20.imag*e21.imag)/(cr2 if abs(cr2)>1e-15 else 1e-15)

        W2[(i1,i2)] += ct0/2; W2[(i2,i1)] += ct0/2
        W2[(i0,i2)] += ct1/2; W2[(i2,i0)] += ct1/2
        W2[(i0,i1)] += ct2/2; W2[(i1,i0)] += ct2/2

    rows, cols, vals = [], [], []
    for (i,j), w in W2.items():
        rows.append(i); cols.append(j); vals.append(w)
    L = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))
    diag = -np.array(L.sum(axis=1)).flatten()
    L = L + sparse.diags(diag)
    return L


def solve_holomorphic(verts, tris, boundary_pairs, sl2r):
    """Solve for x as a Gamma-invariant harmonic function.

    x satisfies: Laplacian(x) = 0 at interior vertices,
    x(vi) = x(vj) for identified boundary pairs.
    """
    (a_real, b_real, b_inv, c_real,
     z_a, z_b, z_c, to_disk, to_H, mobius,
     recenter, recenter_inv, disk_to_H) = sl2r

    N = len(verts)
    L = build_laplacian(verts, tris)

    # Boundary pairs: merge vj into vi (they must have same x value)
    dependent = {}  # vj -> vi
    for vi, vj, G in boundary_pairs:
        if vj not in dependent and vi not in dependent:
            dependent[vj] = vi

    independent = [i for i in range(N) if i not in dependent]
    N_ind = len(independent)
    ind_map = {v: k for k, v in enumerate(independent)}

    print(f"  {N_ind} independent, {len(dependent)} dependent (boundary merged)")

    # Build reduced Laplacian
    rows, cols, vals = [], [], []
    for i in independent:
        ki = ind_map[i]
        for ptr in range(L.indptr[i], L.indptr[i+1]):
            j = L.indices[ptr]
            w = L.data[ptr]
            # Map j to its independent representative
            j_ind = dependent.get(j, j)
            if j_ind in ind_map:
                kj = ind_map[j_ind]
                rows.append(ki); cols.append(kj); vals.append(float(w))

    L_red = sparse.csr_matrix((vals, (rows, cols)), shape=(N_ind, N_ind))

    print(f"  Reduced Laplacian: {N_ind}x{N_ind}, {L_red.nnz} nonzeros")

    # Find smallest eigenvalues (the kernel)
    n_eig = min(10, N_ind - 2)
    print(f"  Computing {n_eig} smallest eigenvalues...")
    eigenvalues, eigenvectors = eigs(L_red.tocsc(), k=n_eig, sigma=0, which='LM')
    order = np.argsort(np.abs(eigenvalues))
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    print("  Eigenvalues:")
    for i in range(min(6, len(eigenvalues))):
        print(f"    lambda_{i} = {eigenvalues[i]:.6e}")

    # Extract first two eigenvectors (hopefully the kernel)
    phi1_red = eigenvectors[:, 0].real  # real Laplacian -> real eigenvectors
    phi2_red = eigenvectors[:, 1].real

    # Expand to full vertex set
    phi1 = np.zeros(N)
    phi2 = np.zeros(N)
    for i in independent:
        k = ind_map[i]
        phi1[i] = phi1_red[k]
        phi2[i] = phi2_red[k]
    for j, vi in dependent.items():
        phi1[j] = phi1[vi]
        phi2[j] = phi2[vi]

    return phi1, phi2, eigenvalues


# ══════════════════════════════════════════════════════════════════════
#  Visualization (simplified — just render the mesh)
# ══════════════════════════════════════════════════════════════════════

def generate_html(verts, tris, phi1, phi2):
    """Simple Three.js visualization colored by phi1."""
    # Use phi1, phi2 as two coordinates + disk position as third
    # This gives a natural 3D embedding of the mesh
    coords = []
    for i, w in enumerate(verts):
        coords.append([float(phi1[i]), float(phi2[i]), float(w.real)])

    # Normalize
    arr = np.array(coords)
    arr -= arr.mean(axis=0)
    scale = np.percentile(np.linalg.norm(arr, axis=1), 95)
    if scale > 0: arr *= 2.0/scale
    coords = arr.tolist()

    tri_list = [[int(a),int(b),int(c)] for a,b,c in tris]

    # Color by phi1 value
    p1 = np.array([c[0] for c in coords])
    p1_norm = (p1 - p1.min()) / (p1.max() - p1.min() + 1e-15)
    hues = p1_norm.tolist()

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Bolza — Attempt 10 (GL(2,3) labels)</title>
<style>*{{margin:0;padding:0}}body{{background:#0a0a12;overflow:hidden;font-family:monospace}}
#info{{position:absolute;top:12px;left:16px;color:#8899aa;font-size:13px;line-height:1.6;pointer-events:none}}
#info h2{{color:#ccd0d8;font-size:16px;margin-bottom:4px;font-weight:normal}}
.n{{color:#e8b84d}}</style></head><body>
<div id="info"><h2>Bolza Surface — Attempt 10</h2>
Holomorphic differentials from Laplacian kernel<br>
Coords: (phi1, phi2, Re(disk))<br>
<span class="n">{len(tri_list)}</span> tris, <span class="n">{len(coords)}</span> verts</div>
<script type="importmap">{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}}}</script>
<script type="module">
import*as THREE from'three';
import{{OrbitControls}}from'three/addons/controls/OrbitControls.js';
const P={json.dumps(coords)};const T={json.dumps(tri_list)};const H={json.dumps(hues)};
const scene=new THREE.Scene();scene.background=new THREE.Color(0x0a0a12);
const cam=new THREE.PerspectiveCamera(45,innerWidth/innerHeight,0.01,100);
cam.position.set(0,4,7);
const ren=new THREE.WebGLRenderer({{antialias:true}});ren.setSize(innerWidth,innerHeight);
ren.setPixelRatio(Math.min(devicePixelRatio,2));document.body.appendChild(ren.domElement);
const ctl=new OrbitControls(cam,ren.domElement);ctl.enableDamping=true;ctl.autoRotate=true;ctl.autoRotateSpeed=0.4;
scene.add(new THREE.AmbientLight(0x556677,0.6));
let l=new THREE.DirectionalLight(0xffeedd,1.2);l.position.set(5,8,5);scene.add(l);
l=new THREE.DirectionalLight(0x99aacc,0.5);l.position.set(-4,2,-3);scene.add(l);
function hsl(h,s,l){{if(s===0)return[l,l,l];const q=l<.5?l*(1+s):l+s-l*s,p=2*l-q;
function f(t){{if(t<0)t+=1;if(t>1)t-=1;return t<1/6?p+(q-p)*6*t:t<.5?q:t<2/3?p+(q-p)*(2/3-t)*6:p}};return[f(h+1/3),f(h),f(h-1/3)]}}
const ps=[],ns=[],cs=[];
for(let i=0;i<T.length;i++){{const[a,b,c]=T[i];
const pa=P[a],pb=P[b],pc=P[c];
const u=[pb[0]-pa[0],pb[1]-pa[1],pb[2]-pa[2]];
const v=[pc[0]-pa[0],pc[1]-pa[1],pc[2]-pa[2]];
const nx=u[1]*v[2]-u[2]*v[1],ny=u[2]*v[0]-u[0]*v[2],nz=u[0]*v[1]-u[1]*v[0];
const nl=Math.sqrt(nx*nx+ny*ny+nz*nz)||1;
for(const idx of[a,b,c]){{ps.push(P[idx][0],P[idx][1],P[idx][2]);ns.push(nx/nl,ny/nl,nz/nl);
const[r,g,bl]=hsl(H[idx],0.7,0.55);cs.push(r,g,bl)}}}}
const geo=new THREE.BufferGeometry();
geo.setAttribute('position',new THREE.Float32BufferAttribute(ps,3));
geo.setAttribute('normal',new THREE.Float32BufferAttribute(ns,3));
geo.setAttribute('color',new THREE.Float32BufferAttribute(cs,3));
scene.add(new THREE.Mesh(geo,new THREE.MeshPhongMaterial({{vertexColors:true,side:THREE.DoubleSide,shininess:30}})));
const wp=[];for(const[a,b,c]of T){{wp.push(P[a][0],P[a][1],P[a][2],P[b][0],P[b][1],P[b][2]);
wp.push(P[b][0],P[b][1],P[b][2],P[c][0],P[c][1],P[c][2]);
wp.push(P[c][0],P[c][1],P[c][2],P[a][0],P[a][1],P[a][2])}}
const wg=new THREE.BufferGeometry();wg.setAttribute('position',new THREE.Float32BufferAttribute(wp,3));
scene.add(new THREE.LineSegments(wg,new THREE.LineBasicMaterial({{color:0x334455,transparent:true,opacity:0.15}})));
window.onresize=()=>{{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();ren.setSize(innerWidth,innerHeight)}};
(function a(){{requestAnimationFrame(a);ctl.update();ren.render(scene,cam)}})();
</script></body></html>"""


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  Bolza Surface — Attempt 10: GL(2,3) Labels")
    print("=" * 60)

    print("\n[1/6] Building GL(2,3) + SL(2,R)...")
    gl23 = build_gl23()
    sl2r = build_sl2r()
    mats, mult, inv, I, orders, ag, bg, cg, bg_inv = gl23

    print("\n[2/6] Full BFS with data retention...")
    all_data = full_bfs(sl2r, gl23, max_tri=5000)
    print(f"  {len(all_data)} triangles generated")

    print("\n[3/6] Proper fundamental domain from GL(2,3) labels...")
    fund = build_fund_domain(all_data)
    # Verify: count triangles per order-8 vertex
    wc_counts = defaultdict(int)
    for g, idx in fund.items():
        wc = all_data[idx][2]
        key = (round(wc.real, 4), round(wc.imag, 4))
        wc_counts[key] += 1
    print(f"  Order-8 vertex distribution:")
    for k, c in sorted(wc_counts.items(), key=lambda x: -x[1]):
        w = complex(*k)
        print(f"    {w:.4f} (|w|={abs(w):.4f}): {c} triangles")

    print("\n[4/6] Side pairings from boundary adjacency...")
    gammas = find_side_pairings(all_data, fund, sl2r, gl23)

    print("\n[5/6] Weierstrass points from hyperelliptic involution...")
    wp_H, wp_disk = find_weierstrass_points(gammas, sl2r)
    print(f"  Found {len(wp_disk)} Weierstrass points:")
    for i, (z, w) in enumerate(zip(wp_H, wp_disk)):
        print(f"    W{i}: disk={w:.6f} |w|={abs(w):.4f}")

    print("\n[6/6] Subdivide + Laplacian solve...")
    SUB_N = 10
    verts, tris = subdivide_fund_domain(all_data, fund, sub_n=SUB_N)
    print(f"  {len(verts)} vertices, {len(tris)} triangles")

    print("  Finding boundary pairs via Gamma elements...")
    pairs = find_boundary_pairs_from_gammas(verts, gammas, sl2r, tol=0.015)
    print(f"  {len(pairs)} boundary pairs matched")

    print("  Solving Laplacian eigenvalue problem...")
    phi1, phi2, eigenvalues = solve_holomorphic(verts, tris, pairs, sl2r)

    # Generate visualization
    outpath = '/Users/robin/git/graphviz/bolza-surface/bolza_r4_v3.html'
    html = generate_html(verts, tris, phi1, phi2)
    with open(outpath, 'w') as f:
        f.write(html)
    print(f"\n  Written to {outpath}")
    print("Done!")


if __name__ == '__main__':
    main()
