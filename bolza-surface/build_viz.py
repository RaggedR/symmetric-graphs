#!/usr/bin/env python3
"""
Bolza surface visualization: Poincaré disk tiling + 3D double torus.

Left: exact (2,3,8) hyperbolic tiling via GL(2,3) Möbius transformations
Right: smooth genus-2 surface (theta graph frame: 2 hubs + 3 tubes)
"""

import numpy as np
import json
from itertools import product
from collections import defaultdict, deque

# ══════════════════════════════════════════════════════════════════════
#  GL(2,3) and (2,3,8) generators
# ══════════════════════════════════════════════════════════════════════

print("Building GL(2,3)...")

def mm3(A, B): return (A @ B) % 3
def mk3(A): return tuple((A % 3).flatten())
def md3(A): return (A[0,0]*A[1,1] - A[0,1]*A[1,0]) % 3
def mi3(A):
    det = md3(A)
    det_inv = 1 if det == 1 else 2
    return (det_inv * np.array([[A[1,1], (-A[0,1])%3], [(-A[1,0])%3, A[0,0]]])) % 3

el = {}
for a, b, c, d in product(range(3), repeat=4):
    A = np.array([[a, b], [c, d]])
    if md3(A) == 0: continue
    k = mk3(A)
    if k not in el: el[k] = A % 3

pm = list(el.values()); pk = list(el.keys())
ki = {k: i for i, k in enumerate(pk)}
assert len(pm) == 48

def gm(i, j): return ki[mk3(mm3(pm[i], pm[j]))]
def ginv(i): return ki[mk3(mi3(pm[i]))]
I = ki[mk3(np.eye(2, dtype=int))]

def order(i):
    c = i
    for n in range(1, 49):
        if c == I: return n
        c = gm(c, i)
ords = [order(i) for i in range(48)]

# Find (2,3,8) generators
ag = bg = cg = None
for a in range(48):
    if ords[a] != 2: continue
    for b in range(48):
        if ords[b] != 3: continue
        if ords[gm(a, b)] == 8:
            ag, bg, cg = a, b, gm(a, b); break
    if ag: break

print(f"  |GL(2,3)| = 48, (2,3,8) generators found")


# ══════════════════════════════════════════════════════════════════════
#  Word decomposition (BFS on Cayley graph)
# ══════════════════════════════════════════════════════════════════════

print("Word decomposition...")

bg_inv = ginv(bg)
words = {I: []}
queue = deque([I])
gens = [(ag, 'a'), (bg, 'b'), (bg_inv, 'B')]

while queue:
    g = queue.popleft()
    for gi, gn in gens:
        h = gm(g, gi)
        if h not in words:
            words[h] = words[g] + [gn]
            queue.append(h)

assert len(words) == 48
print(f"  All 48 elements, max word length: {max(len(w) for w in words.values())}")


# ══════════════════════════════════════════════════════════════════════
#  SL(2,R) generators for (2,3,8) triangle group
# ══════════════════════════════════════════════════════════════════════

print("SL(2,R) generators...")

a_real = np.array([[0, 1], [-1, 0]], dtype=float)  # order 2, rotation by π at i

cos_pi8 = np.cos(np.pi / 8)
sin_pi3 = np.sin(np.pi / 3)
cos_pi3 = np.cos(np.pi / 3)

# b: tr(b) = 2cos(π/3) = 1, tr(ab) = 2cos(π/8)
# sin(π/3)(f + 1/f) = 2cos(π/8)
rhs = 2 * cos_pi8 / sin_pi3
f_param = (rhs + np.sqrt(rhs**2 - 4)) / 2

b_real = np.array([
    [cos_pi3, -f_param * sin_pi3],
    [sin_pi3 / f_param, cos_pi3]
], dtype=float)

b_inv_real = np.linalg.inv(b_real)
c_real = a_real @ b_real

# Verify
assert abs(np.trace(a_real)) < 1e-10
assert abs(np.trace(b_real) - 1) < 1e-10
assert abs(np.trace(c_real) - 2 * cos_pi8) < 1e-10

gen_real = {'a': a_real, 'b': b_real, 'B': b_inv_real}

# Fixed points of generators (fundamental triangle vertices)
def fixed_point(M):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    disc = 4 - (a+d)**2
    if disc < 0: return None
    z1 = ((a - d) + 1j * np.sqrt(disc)) / (2 * c)
    z2 = ((a - d) - 1j * np.sqrt(disc)) / (2 * c)
    return z1 if z1.imag > 0 else z2

z_a = fixed_point(a_real)  # π/2 vertex
z_b = fixed_point(b_real)  # π/3 vertex
z_c = fixed_point(c_real)  # π/8 vertex

print(f"  Fundamental triangle: z_a={z_a:.4f}, z_b={z_b:.4f}, z_c={z_c:.4f}")


# ══════════════════════════════════════════════════════════════════════
#  Poincaré disk tiling
#  Step 1: 48 fundamental domain triangles from GL(2,3) words
#  Step 2: BFS outward through Δ(2,3,8) to fill the disk
# ══════════════════════════════════════════════════════════════════════

print("Computing Poincaré disk tiling...")

def mobius(M, z):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    return (a * z + b) / (c * z + d)

def to_disk(z):
    return (z - 1j) / (z + 1j)

def mat_key(M, tol=1e-6):
    return tuple(round(x / tol) * tol for x in M.flatten())

# Recenter disk at the order-8 vertex (z_c) for canonical octagon
w_c_raw = to_disk(z_c)
print(f"  Recentering disk at order-8 vertex: w_c = {w_c_raw:.4f}")

def recenter(w):
    """Möbius disk automorphism mapping w_c to origin."""
    return (w - w_c_raw) / (1 - np.conj(w_c_raw) * w)

MAX_TRIANGLES = 4000
DISK_CUTOFF = 0.996

# Step 1: BFS to generate ALL triangles in the disk
print("  BFS through Δ(2,3,8) to fill disk...")
gen_mats = [('a', a_real), ('b', b_real), ('B', b_inv_real)]

all_matrices = []
all_triangles_raw = []
all_words = []  # track word for each triangle
visited = set()
queue = deque([(np.eye(2), "")])

while queue and len(all_triangles_raw) < MAX_TRIANGLES:
    M, word = queue.popleft()
    mk = mat_key(M)
    if mk in visited:
        continue
    visited.add(mk)

    try:
        wa = to_disk(mobius(M, z_a))
        wb = to_disk(mobius(M, z_b))
        wc = to_disk(mobius(M, z_c))
    except (ZeroDivisionError, FloatingPointError):
        continue

    if max(abs(wa), abs(wb), abs(wc)) > DISK_CUTOFF:
        continue

    area = abs((wb-wa).real * (wc-wa).imag - (wb-wa).imag * (wc-wa).real)
    if area < 1e-10:
        continue

    all_matrices.append(M)
    all_triangles_raw.append((recenter(wa), recenter(wb), recenter(wc)))
    all_words.append(word)

    for gn, gmat in gen_mats:
        queue.append((M @ gmat, word + gn))

print(f"  {len(all_triangles_raw)} total triangles")

# Step 2: Identify the Dirichlet fundamental domain
# For each triangle, compute hyperbolic distance of its center from origin.
# The "center" is the image of z_a (the order-2 fixed point, which maps to
# the disk origin for the identity element).
# Among all Γ-translates of a triangle, the Dirichlet domain keeps the one
# whose center is closest to the origin.
#
# Strategy: group triangles into Γ-cosets (48 triangles per coset),
# then from each coset pick the triangle closest to origin.

print("  Computing Dirichlet fundamental domain...")

def hyp_dist_from_origin(w):
    """Hyperbolic distance from origin in the Poincaré disk."""
    r = abs(w)
    if r >= 1: return 1e9
    return 2 * np.arctanh(r)

# For each triangle, compute the centroid in the recentered disk
tri_centers = []
for wa, wb, wc in all_triangles_raw:
    tri_centers.append((wa + wb + wc) / 3)

# Group into Γ-cosets: two matrices M1, M2 are in the same coset iff
# M1·M2⁻¹ maps to the identity in GL(2,3).
# Equivalently: they give the same image of the fundamental triangle
# modulo Γ-translation.
#
# Simpler approach: compute GL(2,3) element for each matrix by reducing
# the word modulo the GL(2,3) relations. But we don't have the words
# for BFS-generated matrices.
#
# Even simpler: just pick, for each small region of the disk, the
# Γ-translate closest to the origin. Since the Dirichlet domain is a
# Voronoi cell, each triangle belongs to the fundamental domain iff
# its center is closer to the origin than to any other Γ-translate of
# the origin.
#
# The Γ-translates of the origin are the centers of all copies of
# the fundamental domain. We can find them: they're the images of z_a
# under the 48 GL(2,3) coset representatives composed with each Γ element.
# But we don't know Γ explicitly.
#
# Practical approach: for each triangle, the center wa is in the disk.
# Check if wa is closer to the origin than to any other wa' that is a
# Γ-translate of wa. Two centers are Γ-translates if they differ by a
# Γ element — which means they're in the same Γ-orbit.
#
# All Γ-orbits of the origin are the set of ALL order-2 fixed points
# (images of z_a under Δ). We already have these: tri_centers.
# The Γ-orbit of the origin = {tri_centers[i] : matrix[i] is in Γ·identity}
# = {tri_centers[i] : GL(2,3) element of matrix[i] == identity}
#
# But without knowing which matrices map to which GL(2,3) elements,
# we can't directly identify the Γ-orbit.
#
# SIMPLEST PRACTICAL APPROACH: just check distance from origin.
# A triangle is in the Dirichlet domain iff ALL THREE of its vertices
# are closer to the origin than to any Γ-translate of the origin.
# Since Γ-translates of origin are dense, we approximate:
# a triangle is "fundamental" if its center is within a distance threshold.
#
# Better: the Dirichlet domain has a specific radius. For the Bolza octagon
# with 48 triangles, we know the octagon's circumradius in the Poincaré disk.

# Sort triangles by distance of center from origin
tri_dist = [(i, hyp_dist_from_origin(tri_centers[i])) for i in range(len(all_triangles_raw))]
tri_dist.sort(key=lambda x: x[1])

# The first 48 closest triangles form the Dirichlet fundamental domain
fund_indices = set(idx for idx, _ in tri_dist[:48])
N_FUND = 48

# Reorder: fundamental domain first, then the rest
triangles_disk = []
fund_mats = []
fund_words = []
for idx, _ in tri_dist[:48]:
    triangles_disk.append(all_triangles_raw[idx])
    fund_mats.append(all_matrices[idx])
    fund_words.append(all_words[idx])
for idx, _ in tri_dist[48:]:
    triangles_disk.append(all_triangles_raw[idx])

# Compute GL(2,3) element for each fundamental domain triangle
# by evaluating the word in the GL(2,3) generators
gen_gl = {'a': ag, 'b': bg, 'B': bg_inv}
fund_labels = []
for w in fund_words:
    g = I  # identity in GL(2,3)
    for letter in w:
        g = gm(g, gen_gl[letter])
    mat = pm[g]
    order_g = ords[g]
    # Label: matrix entries as compact string + order
    label = f"{mat[0,0]}{mat[0,1]};{mat[1,0]}{mat[1,1]}"
    fund_labels.append({
        'word': w if w else 'e',
        'matrix': label,
        'order': int(order_g),
    })

# Also compute triangle centers in disk coords for label placement
fund_centers = []
for wa, wb, wc in triangles_disk[:48]:
    cx = (wa + wb + wc) / 3
    fund_centers.append((float(cx.real), float(cx.imag)))

print(f"  Dirichlet domain: {N_FUND} triangles")
print(f"  Max center dist in fund. domain: {tri_dist[47][1]:.4f}")
print(f"  Min center dist outside: {tri_dist[48][1]:.4f}")
print(f"  {len(triangles_disk)} total triangles")
print(f"  Orders in fund. domain: {sorted(set(l['order'] for l in fund_labels))}")

# Index unique vertices
def ckey(z, tol=1e-8):
    return (round(z.real / tol) * tol, round(z.imag / tol) * tol)

vmap = {}; vlist = []
def vidx(z):
    k = ckey(z)
    if k not in vmap:
        vmap[k] = len(vlist)
        vlist.append(z)
    return vmap[k]

itris = []
for wa, wb, wc in triangles_disk:
    itris.append((vidx(wa), vidx(wb), vidx(wc)))

print(f"  {len(vlist)} unique vertices")

# Vertex types by valence
valence = defaultdict(int)
for a, b, c in itris:
    valence[a] += 1; valence[b] += 1; valence[c] += 1

types = defaultdict(int)
for v, val in valence.items():
    if val >= 14: types['V'] += 1      # order 8: 16 triangles meet
    elif val >= 5: types['F'] += 1     # order 3: 6 triangles meet
    else: types['E'] += 1             # order 2: 4 triangles meet
print(f"  Vertex types: {dict(types)}")

# Prepare colors by face (order-3 vertex = face center)
disk_verts = [(float(z.real), float(z.imag)) for z in vlist]
disk_tris = itris

def hsl(h, s, l):
    if s == 0: return (l, l, l)
    def f(p, q, t):
        if t < 0: t += 1
        if t > 1: t -= 1
        if t < 1/6: return p + (q-p)*6*t
        if t < 1/2: return q
        if t < 2/3: return p + (q-p)*(2/3-t)*6
        return p
    q = l*(1+s) if l<0.5 else l+s-l*s; p = 2*l-q
    return (f(p,q,h+1/3), f(p,q,h), f(p,q,h-1/3))

# Color triangles: first 48 are the fundamental domain (brighter),
# the rest are repeated copies (muted)
face_groups = defaultdict(list)
for ti, (ia, ib, ic) in enumerate(disk_tris):
    face_groups[ib].append(ti)

n_face_colors = 16
palette_bright = [hsl((i * 0.618034) % 1.0, 0.65, 0.55) for i in range(n_face_colors)]
palette_muted = [hsl((i * 0.618034) % 1.0, 0.25, 0.35) for i in range(n_face_colors)]

tri_colors = [(0.3, 0.3, 0.3)] * len(disk_tris)
for fi, (fv, tlist) in enumerate(face_groups.items()):
    for ti in tlist:
        is_fundamental = ti < 48  # first 48 triangles = fundamental domain
        palette = palette_bright if is_fundamental else palette_muted
        r, g, b = palette[fi % n_face_colors]
        tri_colors[ti] = (r, g, b)

# Find octagon boundary: edges shared between a first-48 triangle and an exterior triangle
print("  Computing octagon boundary...")
edge_tris = defaultdict(set)
for ti, (ia, ib, ic) in enumerate(disk_tris):
    for e in [(min(ia,ib),max(ia,ib)), (min(ib,ic),max(ib,ic)), (min(ia,ic),max(ia,ic))]:
        edge_tris[e].add(ti)

octagon_edges = []
for e, tris in edge_tris.items():
    has_fund = any(t < 48 for t in tris)
    has_ext = any(t >= 48 for t in tris)
    if has_fund and has_ext:
        octagon_edges.append(e)
print(f"  {len(octagon_edges)} boundary edges")


# ══════════════════════════════════════════════════════════════════════
#  3D genus-2 surface: theta graph (2 hubs + 3 tubes)
# ══════════════════════════════════════════════════════════════════════

print("\nBuilding 3D surface (theta graph frame)...")

# Two hubs connected by 3 tubes in a triangular arrangement
hub_L = np.array([-1.5, 0, 0], dtype=float)  # left hub
hub_R = np.array([1.5, 0, 0], dtype=float)   # right hub

TUBE_R = 0.38
HUB_R = TUBE_R * 1.1
N_ALONG = 28
N_AROUND = 20

# 3 tubes: top, bottom-left, bottom-right (in the YZ plane)
tube_offsets = [
    np.array([0, 0, 0.9]),       # top tube: arcs above
    np.array([0, 0.78, -0.45]),   # bottom-right
    np.array([0, -0.78, -0.45]),  # bottom-left
]

surf_verts = []; surf_normals = []; surf_tris = []

for tube_idx, offset in enumerate(tube_offsets):
    base = len(surf_verts)
    for k in range(N_ALONG + 1):
        t = k / N_ALONG
        # Curved path: arc from hub_L to hub_R through offset
        # Use quadratic Bezier: P(t) = (1-t)²·L + 2t(1-t)·M + t²·R
        mid = 0.5 * (hub_L + hub_R) + offset
        center = (1-t)**2 * hub_L + 2*t*(1-t) * mid + t**2 * hub_R

        # Tangent (derivative of Bezier)
        tangent = 2*(1-t)*(mid - hub_L) + 2*t*(hub_R - mid)
        tangent /= np.linalg.norm(tangent) + 1e-8

        # Perpendicular frame
        up = np.array([0, 0, 1.0])
        if abs(np.dot(tangent, up)) > 0.9:
            up = np.array([0, 1, 0.0])
        n1 = np.cross(tangent, up); n1 /= np.linalg.norm(n1) + 1e-8
        n2 = np.cross(tangent, n1)

        # Bulge at middle for smooth blending
        bulge = 1.0 + 0.12 * np.sin(np.pi * t)

        for m in range(N_AROUND):
            theta = 2 * np.pi * m / N_AROUND
            normal = np.cos(theta) * n1 + np.sin(theta) * n2
            surf_verts.append((center + TUBE_R * bulge * normal).tolist())
            surf_normals.append(normal.tolist())

    for k in range(N_ALONG):
        for m in range(N_AROUND):
            i0 = base + k*N_AROUND + m
            i1 = base + k*N_AROUND + (m+1)%N_AROUND
            i2 = base + (k+1)*N_AROUND + m
            i3 = base + (k+1)*N_AROUND + (m+1)%N_AROUND
            surf_tris.append((i0, i2, i1))
            surf_tris.append((i1, i2, i3))

# Hub spheres
N_SPH_LAT = 12; N_SPH_LON = 20
for hub in [hub_L, hub_R]:
    base = len(surf_verts)
    for lat in range(N_SPH_LAT + 1):
        phi = np.pi * lat / N_SPH_LAT
        for lon in range(N_SPH_LON):
            theta = 2 * np.pi * lon / N_SPH_LON
            normal = np.array([np.sin(phi)*np.cos(theta),
                              np.sin(phi)*np.sin(theta),
                              np.cos(phi)])
            surf_verts.append((hub + HUB_R * normal).tolist())
            surf_normals.append(normal.tolist())
    for lat in range(N_SPH_LAT):
        for lon in range(N_SPH_LON):
            i0 = base + lat*N_SPH_LON + lon
            i1 = base + lat*N_SPH_LON + (lon+1)%N_SPH_LON
            i2 = base + (lat+1)*N_SPH_LON + lon
            i3 = base + (lat+1)*N_SPH_LON + (lon+1)%N_SPH_LON
            surf_tris.append((i0, i2, i1))
            surf_tris.append((i1, i2, i3))

print(f"  {len(surf_verts)} vertices, {len(surf_tris)} triangles")


# ══════════════════════════════════════════════════════════════════════
#  Generate HTML
# ══════════════════════════════════════════════════════════════════════

print("Generating HTML...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8">
<title>Bolza Surface — Poincaré Disk & 3D Surface</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{background:#08080c;overflow:hidden;font-family:'Courier New',monospace;display:flex}}
.panel{{flex:1;position:relative;height:100vh}}
#diskPanel{{border-right:1px solid #1a1a2a}}
canvas{{display:block}}
.label{{position:absolute;top:12px;left:16px;color:#8899aa;font-size:12px;line-height:1.5;pointer-events:none}}
.label h2{{color:#ccd0d8;font-size:15px;margin-bottom:4px;font-weight:normal;letter-spacing:1px}}
.label .n{{color:#e8b84d}}
#info{{position:absolute;bottom:12px;left:16px;right:16px;color:#556677;font-size:11px;text-align:center;pointer-events:none}}
</style></head>
<body>
<div class="panel" id="diskPanel">
  <canvas id="diskCanvas"></canvas>
  <div class="label">
    <h2>Poincaré Disk</h2>
    <span class="n">48</span> fundamental domains<br>
    (2,3,8) hyperbolic tiling<br>
    <span class="n">6</span> vertices · <span class="n">24</span> edges · <span class="n">16</span> faces
  </div>
</div>
<div class="panel" id="threePanel">
  <div class="label">
    <h2>Genus-2 Surface</h2>
    &theta;-graph frame (2 hubs + 3 tubes)<br>
    Aut = GL(2,3), |G| = <span class="n">48</span>
  </div>
</div>
<div id="info">
  Bolza Surface · genus 2 · &chi; = &minus;2 ·
  H&sup1; = V&#x2082; &oplus; V&#x0305;&#x2082; ·
  dual graph: 16-vertex cubic (girth 6)
</div>

<script>
const DV={json.dumps(disk_verts)};
const DT={json.dumps(disk_tris)};
const DC={json.dumps(tri_colors)};
const OE={json.dumps(octagon_edges)};
const FL={json.dumps(fund_labels)};
const FC={json.dumps(fund_centers)};

const dc=document.getElementById('diskCanvas');
const ctx=dc.getContext('2d');

function resizeDisk(){{
  const p=document.getElementById('diskPanel');
  dc.width=p.clientWidth*devicePixelRatio;
  dc.height=p.clientHeight*devicePixelRatio;
  dc.style.width=p.clientWidth+'px';
  dc.style.height=p.clientHeight+'px';
  drawDisk();
}}

function drawDisk(){{
  const w=dc.width,h=dc.height,cx=w/2,cy=h/2,r=Math.min(w,h)*0.44;
  ctx.fillStyle='#08080c';ctx.fillRect(0,0,w,h);
  ctx.beginPath();ctx.arc(cx,cy,r,0,2*Math.PI);ctx.fillStyle='#0d0d14';ctx.fill();
  ctx.strokeStyle='#334455';ctx.lineWidth=1.5;ctx.stroke();

  for(let i=0;i<DT.length;i++){{
    const[ia,ib,ic]=DT[i];const[cr,cg,cb]=DC[i];
    const ax=cx+DV[ia][0]*r,ay=cy-DV[ia][1]*r;
    const bx=cx+DV[ib][0]*r,by=cy-DV[ib][1]*r;
    const ccx=cx+DV[ic][0]*r,ccy=cy-DV[ic][1]*r;
    const maxR=Math.max(DV[ia][0]**2+DV[ia][1]**2,DV[ib][0]**2+DV[ib][1]**2,DV[ic][0]**2+DV[ic][1]**2);
    if(maxR>0.97)continue;
    ctx.beginPath();ctx.moveTo(ax,ay);ctx.lineTo(bx,by);ctx.lineTo(ccx,ccy);ctx.closePath();
    ctx.fillStyle=`rgb(${{Math.round(cr*255)}},${{Math.round(cg*255)}},${{Math.round(cb*255)}})`;
    ctx.fill();ctx.strokeStyle='rgba(8,8,12,0.7)';ctx.lineWidth=0.8;ctx.stroke();
  }}
  // Hard outlines on all 48 fundamental domain triangles
  ctx.strokeStyle='#eedd88';ctx.lineWidth=1.8;ctx.lineCap='round';ctx.lineJoin='round';
  for(let i=0;i<Math.min(48,DT.length);i++){{
    const[ia,ib,ic]=DT[i];
    const ax=cx+DV[ia][0]*r,ay=cy-DV[ia][1]*r;
    const bx=cx+DV[ib][0]*r,by=cy-DV[ib][1]*r;
    const ccx=cx+DV[ic][0]*r,ccy=cy-DV[ic][1]*r;
    ctx.beginPath();ctx.moveTo(ax,ay);ctx.lineTo(bx,by);ctx.lineTo(ccx,ccy);ctx.closePath();ctx.stroke();
  }}

  // Labels on fundamental domain triangles
  ctx.textAlign='center';ctx.textBaseline='middle';
  const fontSize=Math.max(7, r*0.028);
  ctx.font=`${{fontSize}}px Courier New`;
  for(let i=0;i<FC.length;i++){{
    const x=cx+FC[i][0]*r, y=cy-FC[i][1]*r;
    // Only label if triangle is big enough to read
    const dist=Math.sqrt(FC[i][0]**2+FC[i][1]**2);
    if(dist>0.75)continue;
    ctx.fillStyle='rgba(255,255,240,0.85)';
    ctx.fillText(FL[i].matrix, x, y);
  }}

  // Disk boundary circle
  ctx.beginPath();ctx.arc(cx,cy,r,0,2*Math.PI);ctx.strokeStyle='#445566';ctx.lineWidth=2;ctx.stroke();
}}
window.addEventListener('resize',resizeDisk);resizeDisk();
</script>

<script type="importmap">
{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}}}
</script>
<script type="module">
import*as THREE from'three';
import{{OrbitControls}}from'three/addons/controls/OrbitControls.js';

const SV={json.dumps(surf_verts)};
const SN={json.dumps(surf_normals)};
const ST={json.dumps([(int(a),int(b),int(c)) for a,b,c in surf_tris])};

const container=document.getElementById('threePanel');
const scene=new THREE.Scene();scene.background=new THREE.Color(0x08080c);
const cam=new THREE.PerspectiveCamera(38,container.clientWidth/container.clientHeight,0.1,100);
cam.position.set(4,2.5,4);
const ren=new THREE.WebGLRenderer({{antialias:true}});
ren.setSize(container.clientWidth,container.clientHeight);
ren.setPixelRatio(Math.min(devicePixelRatio,2));
ren.toneMapping=THREE.ACESFilmicToneMapping;ren.toneMappingExposure=1.2;
container.appendChild(ren.domElement);
const ctl=new OrbitControls(cam,ren.domElement);
ctl.enableDamping=true;ctl.dampingFactor=0.05;ctl.autoRotate=true;ctl.autoRotateSpeed=0.5;

scene.add(new THREE.AmbientLight(0x445566,0.7));
let l;
l=new THREE.DirectionalLight(0xffeedd,1.3);l.position.set(5,8,5);scene.add(l);
l=new THREE.DirectionalLight(0x99aacc,0.5);l.position.set(-4,2,-3);scene.add(l);
l=new THREE.DirectionalLight(0xaabbdd,0.3);l.position.set(0,-4,5);scene.add(l);

const sp=[],sn=[];
for(const[a,b,c]of ST)for(const idx of[a,b,c]){{sp.push(...SV[idx]);sn.push(...SN[idx])}}
const sg=new THREE.BufferGeometry();
sg.setAttribute('position',new THREE.Float32BufferAttribute(sp,3));
sg.setAttribute('normal',new THREE.Float32BufferAttribute(sn,3));
scene.add(new THREE.Mesh(sg,new THREE.MeshPhongMaterial({{
  color:0xc49545,specular:0x554433,shininess:50,side:THREE.DoubleSide
}})));

window.addEventListener('resize',()=>{{
  cam.aspect=container.clientWidth/container.clientHeight;
  cam.updateProjectionMatrix();ren.setSize(container.clientWidth,container.clientHeight);
  resizeDisk();
}});
(function a(){{requestAnimationFrame(a);ctl.update();ren.render(scene,cam)}})();
</script>
</body></html>"""

with open('/Users/robin/git/graphviz/bolza-surface/bolza_3d.html', 'w') as f:
    f.write(html)
print("Done! → bolza_3d.html")
