#!/usr/bin/env python3
"""
Klein quartic visualization: Poincaré disk tiling + 3D tube surface.

Left panel: exact hyperbolic tiling via (2,3,7) Möbius transformations
Right panel: smooth genus-3 tube surface (rotating)
"""

import numpy as np
import json
from itertools import product, combinations
from collections import defaultdict, deque

# ══════════════════════════════════════════════════════════════════════
#  PART 1: PSL(2,7) as matrices over F_7
# ══════════════════════════════════════════════════════════════════════

print("Building PSL(2,7)...")

def mm7(A, B): return (A @ B) % 7
def mk(A):
    A = A % 7; n = (7 - A) % 7
    return min(tuple(A.flatten()), tuple(n.flatten()))

el = {}
for a, b, c, d in product(range(7), repeat=4):
    A = np.array([[a, b], [c, d]])
    if (A[0,0]*A[1,1] - A[0,1]*A[1,0]) % 7 != 1: continue
    k = mk(A)
    if k not in el: el[k] = A % 7

pm = list(el.values()); pk = list(el.keys())
ki = {k: i for i, k in enumerate(pk)}
assert len(pm) == 168

def gm(i, j): return ki[mk(mm7(pm[i], pm[j]))]
I = ki[mk(np.eye(2, dtype=int))]

def order(i):
    c = i
    for n in range(1, 169):
        if c == I: return n
        c = gm(c, i)
    return -1

ords = [order(i) for i in range(168)]

# Find (2,3,7) generators
ag = bg = cg = None
for a in range(168):
    if ords[a] != 2: continue
    for b in range(168):
        if ords[b] != 3: continue
        if ords[gm(a, b)] == 7:
            ag, bg, cg = a, b, gm(a, b); break
    if ag is not None: break

print(f"  (2,3,7) generators: a={ag}, b={bg}, c=ab={cg}")


# ══════════════════════════════════════════════════════════════════════
#  PART 2: Word decomposition (BFS on Cayley graph)
# ══════════════════════════════════════════════════════════════════════

print("Computing word decomposition...")

# b inverse in PSL(2,7)
bg_inv = None
for i in range(168):
    if gm(bg, i) == I:
        bg_inv = i; break

# BFS to find shortest word for each element
words = {I: []}  # element -> list of generator indices
queue = deque([I])
gen_list = [(ag, 'a'), (bg, 'b'), (bg_inv, 'B')]

while queue:
    g = queue.popleft()
    for gen_idx, gen_name in gen_list:
        h = gm(g, gen_idx)
        if h not in words:
            words[h] = words[g] + [gen_name]
            queue.append(h)

assert len(words) == 168
max_len = max(len(w) for w in words.values())
print(f"  All 168 elements reached, max word length: {max_len}")


# ══════════════════════════════════════════════════════════════════════
#  PART 3: Real SL(2,R) matrices for (2,3,7) generators
# ══════════════════════════════════════════════════════════════════════

print("Building SL(2,R) generators...")

# Generator a: rotation by π at i
a_real = np.array([[0, 1], [-1, 0]], dtype=float)

# Generator b: rotation by 2π/3
cos_pi7 = np.cos(np.pi / 7)
sin_pi3 = np.sin(np.pi / 3)
cos_pi3 = np.cos(np.pi / 3)
f_param = (2 * cos_pi7 / sin_pi3 + np.sqrt((2 * cos_pi7 / sin_pi3)**2 - 4)) / 2

b_real = np.array([
    [cos_pi3, -f_param * sin_pi3],
    [sin_pi3 / f_param, cos_pi3]
], dtype=float)

b_inv_real = np.linalg.inv(b_real)

# Verify
c_real = a_real @ b_real
assert abs(np.trace(a_real)) < 1e-10, f"tr(a) = {np.trace(a_real)}"
assert abs(np.trace(b_real) - 1) < 1e-10, f"tr(b) = {np.trace(b_real)}"
assert abs(np.trace(c_real) - 2 * cos_pi7) < 1e-10, f"tr(c) = {np.trace(c_real)}"

gen_real = {'a': a_real, 'b': b_real, 'B': b_inv_real}


# ══════════════════════════════════════════════════════════════════════
#  PART 4: Compute fundamental triangle and all 168 images
# ══════════════════════════════════════════════════════════════════════

print("Computing Poincaré disk tiling...")

def fixed_point_upper(M):
    """Fixed point of elliptic M in upper half-plane."""
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    tr = a + d
    disc = 4 - tr**2
    if disc < 0: return None
    # Two roots; pick the one in upper half-plane
    z1 = ((a - d) + 1j * np.sqrt(disc)) / (2 * c)
    z2 = ((a - d) - 1j * np.sqrt(disc)) / (2 * c)
    if z1.imag > 0: return z1
    if z2.imag > 0: return z2
    return z1  # fallback

# Fixed points of generators = vertices of fundamental triangle
z_a = fixed_point_upper(a_real)   # order-2 point (edge midpoint)
z_b = fixed_point_upper(b_real)   # order-3 point (face center)
z_c = fixed_point_upper(c_real)   # order-7 point (vertex)

print(f"  Fundamental triangle in ℍ:")
print(f"    z_a (π/2) = {z_a:.6f}")
print(f"    z_b (π/3) = {z_b:.6f}")
print(f"    z_c (π/7) = {z_c:.6f}")

# Möbius transformation
def mobius(M, z):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    return (a * z + b) / (c * z + d)

# Cayley transform: ℍ → 𝔻 (unit disk)
def to_disk(z):
    return (z - 1j) / (z + 1j)

# Map fundamental triangle to disk
w_a = to_disk(z_a)
w_b = to_disk(z_b)
w_c = to_disk(z_c)

print(f"  Fundamental triangle in 𝔻:")
print(f"    w_a = {w_a:.6f}  (|w|={abs(w_a):.4f})")
print(f"    w_b = {w_b:.6f}  (|w|={abs(w_b):.4f})")
print(f"    w_c = {w_c:.6f}  (|w|={abs(w_c):.4f})")

# Compute real matrix for each group element
def word_to_matrix(word_list):
    M = np.eye(2)
    for letter in word_list:
        M = M @ gen_real[letter]
    return M

# Apply all 168 elements to get 168 copies of the fundamental triangle
# Each copy has 3 vertices: image of (z_a, z_b, z_c)
triangles_disk = []  # list of (w_a', w_b', w_c') in disk coordinates

for g_idx in range(168):
    word = words[g_idx]
    M = word_to_matrix(word)
    za_img = mobius(M, z_a)
    zb_img = mobius(M, z_b)
    zc_img = mobius(M, z_c)
    wa = to_disk(za_img)
    wb = to_disk(zb_img)
    wc = to_disk(zc_img)
    triangles_disk.append((wa, wb, wc))

# Check: all vertices should be inside the unit disk
max_r = max(max(abs(wa), abs(wb), abs(wc)) for wa, wb, wc in triangles_disk)
print(f"  {len(triangles_disk)} triangles, max |w| = {max_r:.6f}")

# Collect unique vertices (with tolerance)
def complex_key(z, tol=1e-8):
    return (round(z.real / tol) * tol, round(z.imag / tol) * tol)

vertex_map = {}  # complex_key -> index
vertex_list = []  # list of complex numbers

def get_vertex_idx(z):
    k = complex_key(z)
    if k not in vertex_map:
        vertex_map[k] = len(vertex_list)
        vertex_list.append(z)
    return vertex_map[k]

indexed_triangles = []
for wa, wb, wc in triangles_disk:
    ia = get_vertex_idx(wa)
    ib = get_vertex_idx(wb)
    ic = get_vertex_idx(wc)
    indexed_triangles.append((ia, ib, ic))

print(f"  {len(vertex_list)} unique vertices")

# Determine vertex types by counting how many triangles meet
vertex_valence = defaultdict(int)
for ia, ib, ic in indexed_triangles:
    vertex_valence[ia] += 1
    vertex_valence[ib] += 1
    vertex_valence[ic] += 1

# In the (2,3,7) tiling:
# Order-7 vertices (V-type): 14 triangles meet (7 oriented pairs)
# Order-3 vertices (F-type): 6 triangles meet (3 oriented pairs)
# Order-2 vertices (E-type): 4 triangles meet (2 oriented pairs)
type_counts = defaultdict(int)
vertex_types = {}
for v, val in vertex_valence.items():
    if val >= 12:  # 14 = order 7
        vertex_types[v] = 'V'
    elif val >= 5:   # 6 = order 3
        vertex_types[v] = 'F'
    else:            # 4 = order 2
        vertex_types[v] = 'E'
    type_counts[vertex_types[v]] += 1

print(f"  Vertex types: {dict(type_counts)}")
# Should be V=24, E=84, F=56 but tiling in disk may show more
# (the 14-gon has boundary identifications we haven't enforced)


# ══════════════════════════════════════════════════════════════════════
#  PART 5: Prepare disk data for rendering
# ══════════════════════════════════════════════════════════════════════

# Convert to JSON-serializable format
disk_verts = [(float(z.real), float(z.imag)) for z in vertex_list]
disk_tris = [(int(a), int(b), int(c)) for a, b, c in indexed_triangles]

# Assign colors: group triangles by which face they belong to
# Two triangles sharing the F-vertex (order-3) are in the same {3,7} face
# Group triangles by their F-vertex (index ib in each triangle)
face_groups = defaultdict(list)
for tri_idx, (ia, ib, ic) in enumerate(indexed_triangles):
    # ib is the order-3 vertex (face center)
    face_groups[ib].append(tri_idx)

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

tri_colors = [(0.5, 0.5, 0.5)] * len(disk_tris)
for fi, (fv, tri_list) in enumerate(face_groups.items()):
    h = (fi * 0.618034) % 1.0
    r, g, b = hsl(h, 0.55, 0.5)
    for ti in tri_list:
        tri_colors[ti] = (r, g, b)


# ══════════════════════════════════════════════════════════════════════
#  PART 6: Build smooth 3D tube surface
# ══════════════════════════════════════════════════════════════════════

print("\nBuilding 3D tube surface...")

R_frame = 2.0
tet = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], dtype=float)
tet *= R_frame / np.linalg.norm(tet[0])

TUBE_R = 0.42
HUB_R = TUBE_R * 1.0
N_ALONG = 24
N_AROUND = 20

def tube_frame_fn(p0, p1):
    td = p1 - p0; td /= np.linalg.norm(td)
    up = np.array([0.,0.,1.])
    if abs(np.dot(td, up)) > 0.9: up = np.array([1.,0.,0.])
    n1 = np.cross(td, up); n1 /= np.linalg.norm(n1)
    n2 = np.cross(td, n1)
    return td, n1, n2

surf_verts = []; surf_normals = []; surf_tris = []
tube_edges_3d = [(i,j) for i in range(4) for j in range(i+1,4)]

for ti, tj in tube_edges_3d:
    p0, p1 = tet[ti], tet[tj]
    td, n1, n2 = tube_frame_fn(p0, p1)
    base = len(surf_verts)
    for k in range(N_ALONG + 1):
        t = k / N_ALONG
        bulge = 1.0 + 0.15 * np.sin(np.pi * t)
        center = (1-t)*p0 + t*p1
        for m in range(N_AROUND):
            theta = 2*np.pi*m/N_AROUND
            normal = np.cos(theta)*n1 + np.sin(theta)*n2
            surf_verts.append((center + TUBE_R*bulge*normal).tolist())
            surf_normals.append(normal.tolist())
    for k in range(N_ALONG):
        for m in range(N_AROUND):
            i0 = base + k*N_AROUND + m
            i1 = base + k*N_AROUND + (m+1)%N_AROUND
            i2 = base + (k+1)*N_AROUND + m
            i3 = base + (k+1)*N_AROUND + (m+1)%N_AROUND
            surf_tris.append((i0, i2, i1)); surf_tris.append((i1, i2, i3))

N_SPH_LAT = 10; N_SPH_LON = 20
for hi in range(4):
    base = len(surf_verts)
    center = tet[hi]
    for lat in range(N_SPH_LAT+1):
        phi = np.pi*lat/N_SPH_LAT
        for lon in range(N_SPH_LON):
            theta = 2*np.pi*lon/N_SPH_LON
            normal = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])
            surf_verts.append((center + HUB_R*normal).tolist())
            surf_normals.append(normal.tolist())
    for lat in range(N_SPH_LAT):
        for lon in range(N_SPH_LON):
            i0 = base + lat*N_SPH_LON + lon
            i1 = base + lat*N_SPH_LON + (lon+1)%N_SPH_LON
            i2 = base + (lat+1)*N_SPH_LON + lon
            i3 = base + (lat+1)*N_SPH_LON + (lon+1)%N_SPH_LON
            surf_tris.append((i0, i2, i1)); surf_tris.append((i1, i2, i3))

print(f"  {len(surf_verts)} vertices, {len(surf_tris)} triangles")


# ══════════════════════════════════════════════════════════════════════
#  PART 7: Generate HTML
# ══════════════════════════════════════════════════════════════════════

print("\nGenerating HTML...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8">
<title>Klein Quartic — Poincaré Disk & 3D Surface</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ background:#08080c; overflow:hidden; font-family:'Courier New',monospace; display:flex; }}
.panel {{ flex:1; position:relative; height:100vh; }}
#diskPanel {{ border-right: 1px solid #1a1a2a; }}
canvas {{ display:block; }}
.label {{
  position:absolute; top:12px; left:16px;
  color:#8899aa; font-size:12px; line-height:1.5;
  pointer-events:none;
}}
.label h2 {{ color:#ccd0d8; font-size:15px; margin-bottom:4px; font-weight:normal; letter-spacing:1px; }}
.label .n {{ color:#e8b84d; }}
#info {{
  position:absolute; bottom:12px; left:16px; right:16px;
  color:#556677; font-size:11px; text-align:center;
  pointer-events:none;
}}
</style></head>
<body>

<!-- Left: Poincaré disk -->
<div class="panel" id="diskPanel">
  <canvas id="diskCanvas"></canvas>
  <div class="label">
    <h2>Poincaré Disk</h2>
    <span class="n">168</span> fundamental domains<br>
    (2,3,7) hyperbolic tiling
  </div>
</div>

<!-- Right: 3D surface -->
<div class="panel" id="threePanel">
  <div class="label">
    <h2>Genus-3 Surface</h2>
    thickened K&#x2084; frame<br>
    Aut = PSL(2,7)
  </div>
</div>

<div id="info">
  Klein Quartic &middot; genus 3 &middot; &chi; = &minus;4 &middot;
  H&sup1; = V&#x2083; &oplus; V&#x0305;&#x2083; &middot;
  |Aut| = <span class="n" style="color:#e8b84d">168</span>
</div>

<script type="importmap">
{{
  "imports": {{
    "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"
  }}
}}
</script>

<script>
// ══════════════════════════════════════════════════════════════════
//  Poincaré Disk (2D Canvas)
// ══════════════════════════════════════════════════════════════════

const DV = {json.dumps(disk_verts)};
const DT = {json.dumps(disk_tris)};
const DC = {json.dumps(tri_colors)};

const diskCanvas = document.getElementById('diskCanvas');
const ctx = diskCanvas.getContext('2d');

function resizeDisk() {{
  const panel = document.getElementById('diskPanel');
  diskCanvas.width = panel.clientWidth * devicePixelRatio;
  diskCanvas.height = panel.clientHeight * devicePixelRatio;
  diskCanvas.style.width = panel.clientWidth + 'px';
  diskCanvas.style.height = panel.clientHeight + 'px';
  drawDisk();
}}

function drawDisk() {{
  const w = diskCanvas.width, h = diskCanvas.height;
  const cx = w/2, cy = h/2;
  const r = Math.min(w, h) * 0.44;

  ctx.fillStyle = '#08080c';
  ctx.fillRect(0, 0, w, h);

  // Disk boundary
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, 2*Math.PI);
  ctx.fillStyle = '#0d0d14';
  ctx.fill();
  ctx.strokeStyle = '#334455';
  ctx.lineWidth = 1.5;
  ctx.stroke();

  // Draw triangles
  for (let i = 0; i < DT.length; i++) {{
    const [ia, ib, ic] = DT[i];
    const [cr, cg, cb] = DC[i];

    const ax = cx + DV[ia][0]*r, ay = cy - DV[ia][1]*r;
    const bx = cx + DV[ib][0]*r, by = cy - DV[ib][1]*r;
    const ccx2 = cx + DV[ic][0]*r, ccy = cy - DV[ic][1]*r;

    // Skip triangles too close to disk edge (they get infinitesimally small)
    const maxR = Math.max(
      DV[ia][0]**2 + DV[ia][1]**2,
      DV[ib][0]**2 + DV[ib][1]**2,
      DV[ic][0]**2 + DV[ic][1]**2
    );
    if (maxR > 0.98) continue;

    ctx.beginPath();
    ctx.moveTo(ax, ay);
    ctx.lineTo(bx, by);
    ctx.lineTo(ccx2, ccy);
    ctx.closePath();

    const R = Math.round(cr * 255), G = Math.round(cg * 255), B = Math.round(cb * 255);
    ctx.fillStyle = `rgb(${{R}},${{G}},${{B}})`;
    ctx.fill();
    ctx.strokeStyle = 'rgba(8,8,12,0.6)';
    ctx.lineWidth = 0.8;
    ctx.stroke();
  }}

  // Redraw boundary circle on top
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, 2*Math.PI);
  ctx.strokeStyle = '#445566';
  ctx.lineWidth = 2;
  ctx.stroke();
}}

window.addEventListener('resize', resizeDisk);
resizeDisk();
</script>

<script type="module">
// ══════════════════════════════════════════════════════════════════
//  3D Tube Surface (Three.js)
// ══════════════════════════════════════════════════════════════════

import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';

const SV = {json.dumps(surf_verts)};
const SN = {json.dumps(surf_normals)};
const ST = {json.dumps([(int(a),int(b),int(c)) for a,b,c in surf_tris])};

const container = document.getElementById('threePanel');
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x08080c);

const cam = new THREE.PerspectiveCamera(38, container.clientWidth/container.clientHeight, 0.1, 100);
cam.position.set(5.5, 3.5, 5.5);

const ren = new THREE.WebGLRenderer({{ antialias: true }});
ren.setSize(container.clientWidth, container.clientHeight);
ren.setPixelRatio(Math.min(devicePixelRatio, 2));
ren.toneMapping = THREE.ACESFilmicToneMapping;
ren.toneMappingExposure = 1.2;
container.appendChild(ren.domElement);

const ctl = new OrbitControls(cam, ren.domElement);
ctl.enableDamping = true; ctl.dampingFactor = 0.05;
ctl.autoRotate = true; ctl.autoRotateSpeed = 0.5;

scene.add(new THREE.AmbientLight(0x445566, 0.7));
let l;
l = new THREE.DirectionalLight(0xffeedd, 1.3); l.position.set(5,8,5); scene.add(l);
l = new THREE.DirectionalLight(0x99aacc, 0.5); l.position.set(-4,2,-3); scene.add(l);
l = new THREE.DirectionalLight(0xaabbdd, 0.3); l.position.set(0,-4,5); scene.add(l);

const sp = [], sn = [];
for (const [a,b,c] of ST) {{
  for (const idx of [a,b,c]) {{ sp.push(...SV[idx]); sn.push(...SN[idx]); }}
}}
const sg = new THREE.BufferGeometry();
sg.setAttribute('position', new THREE.Float32BufferAttribute(sp, 3));
sg.setAttribute('normal', new THREE.Float32BufferAttribute(sn, 3));
scene.add(new THREE.Mesh(sg, new THREE.MeshPhongMaterial({{
  color: 0xc49545, specular: 0x554433, shininess: 50, side: THREE.DoubleSide
}})));

window.addEventListener('resize', () => {{
  cam.aspect = container.clientWidth / container.clientHeight;
  cam.updateProjectionMatrix();
  ren.setSize(container.clientWidth, container.clientHeight);
  resizeDisk();
}});

(function animate() {{
  requestAnimationFrame(animate);
  ctl.update();
  ren.render(scene, cam);
}})();
</script>
</body>
</html>"""

with open('/Users/robin/git/graphviz/klein-quartic/klein_quartic_3d.html', 'w') as f:
    f.write(html)
print("Done! → klein_quartic_3d.html")
