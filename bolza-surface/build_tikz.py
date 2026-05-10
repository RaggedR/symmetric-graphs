#!/usr/bin/env python3
"""
Generate TikZ Poincaré disk of the Bolza surface (2,3,8) tiling.
Outputs: bolza_disk.tex → bolza_disk.pdf → bolza_disk.jpg
"""

import numpy as np
from itertools import product
from collections import defaultdict, deque

# ══════════════════════════════════════════════════════════════════════
#  Build GL(2,3) and tiling (same as build_viz.py)
# ══════════════════════════════════════════════════════════════════════

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
def gm(i, j): return ki[mk3(mm3(pm[i], pm[j]))]
def ginv(i): return ki[mk3(mi3(pm[i]))]
I = ki[mk3(np.eye(2, dtype=int))]
def order(i):
    c = i
    for n in range(1, 49):
        if c == I: return n
        c = gm(c, i)
ords = [order(i) for i in range(48)]

ag = bg = cg = None
for a in range(48):
    if ords[a] != 2: continue
    for b in range(48):
        if ords[b] != 3: continue
        if ords[gm(a, b)] == 8:
            ag, bg, cg = a, b, gm(a, b); break
    if ag: break

bg_inv = ginv(bg)
words = {I: []}
queue = deque([I])
gens_gl = [(ag, 'a'), (bg, 'b'), (bg_inv, 'B')]
while queue:
    g = queue.popleft()
    for gi, gn in gens_gl:
        h = gm(g, gi)
        if h not in words:
            words[h] = words[g] + [gn]
            queue.append(h)

# SL(2,R) generators
a_real = np.array([[0, 1], [-1, 0]], dtype=float)
cos_pi8 = np.cos(np.pi / 8)
sin_pi3 = np.sin(np.pi / 3)
cos_pi3 = np.cos(np.pi / 3)
rhs = 2 * cos_pi8 / sin_pi3
f_param = (rhs + np.sqrt(rhs**2 - 4)) / 2
b_real = np.array([[cos_pi3, -f_param * sin_pi3], [sin_pi3 / f_param, cos_pi3]])
b_inv_real = np.linalg.inv(b_real)
c_real = a_real @ b_real
gen_real = {'a': a_real, 'b': b_real, 'B': b_inv_real}

def fixed_point(M):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    disc = 4 - (a+d)**2
    z1 = ((a - d) + 1j * np.sqrt(disc)) / (2 * c)
    z2 = ((a - d) - 1j * np.sqrt(disc)) / (2 * c)
    return z1 if z1.imag > 0 else z2

z_a = fixed_point(a_real)
z_b = fixed_point(b_real)
z_c = fixed_point(c_real)

def mobius(M, z):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    return (a * z + b) / (c * z + d)

def to_disk(z):
    return (z - 1j) / (z + 1j)

# Recenter at order-8 vertex
w_c_raw = to_disk(z_c)
def recenter(w):
    return (w - w_c_raw) / (1 - np.conj(w_c_raw) * w)

def mat_key(M, tol=1e-6):
    return tuple(round(x / tol) * tol for x in M.flatten())

def hyp_dist(w):
    r = abs(w)
    if r >= 1: return 1e9
    return 2 * np.arctanh(r)

# Generate all triangles
print("Generating tiling...")
MAX_TRIANGLES = 5000
DISK_CUTOFF = 0.997

all_tris = []
all_words_list = []
visited = set()
queue = deque([(np.eye(2), "")])

while queue and len(all_tris) < MAX_TRIANGLES:
    M, word = queue.popleft()
    mk = mat_key(M)
    if mk in visited: continue
    visited.add(mk)
    try:
        wa = recenter(to_disk(mobius(M, z_a)))
        wb = recenter(to_disk(mobius(M, z_b)))
        wc = recenter(to_disk(mobius(M, z_c)))
    except:
        continue
    if max(abs(wa), abs(wb), abs(wc)) > DISK_CUTOFF: continue
    area = abs((wb-wa).real * (wc-wa).imag - (wb-wa).imag * (wc-wa).real)
    if area < 1e-10: continue
    all_tris.append((wa, wb, wc))
    all_words_list.append(word)
    for gn, gmat in [('a', a_real), ('b', b_real), ('B', b_inv_real)]:
        queue.append((M @ gmat, word + gn))

print(f"  {len(all_tris)} triangles")

# Dirichlet domain: 48 closest to origin
tri_dist = [(i, hyp_dist((all_tris[i][0]+all_tris[i][1]+all_tris[i][2])/3))
            for i in range(len(all_tris))]
tri_dist.sort(key=lambda x: x[1])
fund_set = set(idx for idx, _ in tri_dist[:48])

# Compute GL(2,3) labels for fundamental triangles
gen_gl_map = {'a': ag, 'b': bg, 'B': bg_inv}
fund_labels = {}
for idx, _ in tri_dist[:48]:
    w = all_words_list[idx]
    g = I
    for letter in w:
        g = gm(g, gen_gl_map[letter])
    mat = pm[g]
    fund_labels[idx] = f"\\tiny ${mat[0,0]}{mat[0,1]}\\atop {mat[1,0]}{mat[1,1]}$"

# Boundary edges (shared between fund and non-fund)
def ckey(z, tol=1e-6):
    return (round(z.real/tol)*tol, round(z.imag/tol)*tol)
vmap = {}; vlist = []
def vidx(z):
    k = ckey(z)
    if k not in vmap:
        vmap[k] = len(vlist); vlist.append(z)
    return vmap[k]

itris = []
for wa, wb, wc in all_tris:
    itris.append((vidx(wa), vidx(wb), vidx(wc)))

edge_tris = defaultdict(set)
for ti, (ia, ib, ic) in enumerate(itris):
    for e in [(min(ia,ib),max(ia,ib)), (min(ib,ic),max(ib,ic)), (min(ia,ic),max(ia,ic))]:
        edge_tris[e].add(ti)

boundary_edges = []
for e, tris in edge_tris.items():
    has_f = any(t in fund_set for t in tris)
    has_x = any(t not in fund_set for t in tris)
    if has_f and has_x:
        boundary_edges.append(e)


# ══════════════════════════════════════════════════════════════════════
#  Generate TikZ
# ══════════════════════════════════════════════════════════════════════

print("Generating TikZ...")

SCALE = 6.0  # cm radius of the disk

def hsl_to_tikz(h, s, l):
    """Convert HSL to RGB for TikZ."""
    if s == 0: return f"{{{l:.2f},{l:.2f},{l:.2f}}}"
    def f(p, q, t):
        if t < 0: t += 1
        if t > 1: t -= 1
        if t < 1/6: return p + (q-p)*6*t
        if t < 1/2: return q
        if t < 2/3: return p + (q-p)*(2/3-t)*6
        return p
    q = l*(1+s) if l<0.5 else l+s-l*s; p = 2*l-q
    r, g, b = f(p,q,h+1/3), f(p,q,h), f(p,q,h-1/3)
    return f"{{{r:.3f},{g:.3f},{b:.3f}}}"

# Color definitions
n_colors = 16
color_defs = []
for i in range(n_colors):
    h = (i * 0.618034) % 1.0
    rgb_bright = hsl_to_tikz(h, 0.65, 0.55)
    rgb_muted = hsl_to_tikz(h, 0.20, 0.30)
    color_defs.append(f"\\definecolor{{fb{i}}}{{rgb}}{rgb_bright}")
    color_defs.append(f"\\definecolor{{fm{i}}}{{rgb}}{rgb_muted}")

# Group triangles by face center for coloring
face_groups = defaultdict(list)
for ti, (ia, ib, ic) in enumerate(itris):
    face_groups[ib].append(ti)

tri_color_names = ["fm0"] * len(itris)
for fi, (fv, tlist) in enumerate(face_groups.items()):
    ci = fi % n_colors
    for ti in tlist:
        if ti in fund_set:
            tri_color_names[ti] = f"fb{ci}"
        else:
            tri_color_names[ti] = f"fm{ci}"

# Build TikZ content
lines = []
lines.append(r"\documentclass[border=5mm]{standalone}")
lines.append(r"\usepackage{tikz}")
lines.append(r"\begin{document}")
lines.append("")
for cd in color_defs:
    lines.append(cd)
lines.append("")
lines.append(f"\\begin{{tikzpicture}}[scale=1]")
lines.append(f"  % Background disk")
lines.append(f"  \\fill[black!95] (0,0) circle ({SCALE}cm);")
lines.append(f"  \\draw[gray!40, line width=0.8pt] (0,0) circle ({SCALE}cm);")
lines.append("")

# Draw all triangles (muted first, then bright on top)
lines.append(f"  % Muted triangles (translates)")
for ti in range(len(itris)):
    if ti in fund_set: continue
    ia, ib, ic = itris[ti]
    wa, wb, wc = vlist[ia], vlist[ib], vlist[ic]
    ax, ay = wa.real * SCALE, wa.imag * SCALE
    bx, by = wb.real * SCALE, wb.imag * SCALE
    cx, cy = wc.real * SCALE, wc.imag * SCALE
    cn = tri_color_names[ti]
    lines.append(f"  \\fill[{cn}, draw=black!80, line width=0.1pt] "
                 f"({ax:.4f},{ay:.4f}) -- ({bx:.4f},{by:.4f}) -- ({cx:.4f},{cy:.4f}) -- cycle;")

lines.append("")
lines.append(f"  % Fundamental domain (48 bright triangles)")
for ti in range(len(itris)):
    if ti not in fund_set: continue
    ia, ib, ic = itris[ti]
    wa, wb, wc = vlist[ia], vlist[ib], vlist[ic]
    ax, ay = wa.real * SCALE, wa.imag * SCALE
    bx, by = wb.real * SCALE, wb.imag * SCALE
    cx, cy = wc.real * SCALE, wc.imag * SCALE
    cn = tri_color_names[ti]
    lines.append(f"  \\fill[{cn}, draw=black!90, line width=0.2pt] "
                 f"({ax:.4f},{ay:.4f}) -- ({bx:.4f},{by:.4f}) -- ({cx:.4f},{cy:.4f}) -- cycle;")

# Gold outlines on fundamental domain
lines.append("")
lines.append(f"  % Gold outlines on fundamental domain")
for ti in range(len(itris)):
    if ti not in fund_set: continue
    ia, ib, ic = itris[ti]
    wa, wb, wc = vlist[ia], vlist[ib], vlist[ic]
    ax, ay = wa.real * SCALE, wa.imag * SCALE
    bx, by = wb.real * SCALE, wb.imag * SCALE
    cx, cy = wc.real * SCALE, wc.imag * SCALE
    lines.append(f"  \\draw[yellow!60!orange, line width=0.5pt] "
                 f"({ax:.4f},{ay:.4f}) -- ({bx:.4f},{by:.4f}) -- ({cx:.4f},{cy:.4f}) -- cycle;")

# Labels
lines.append("")
lines.append(f"  % GL(2,3) matrix labels")
for ti in fund_set:
    if ti not in fund_labels: continue
    ia, ib, ic = itris[ti]
    wa, wb, wc = vlist[ia], vlist[ib], vlist[ic]
    cx = (wa.real + wb.real + wc.real) / 3 * SCALE
    cy = (wa.imag + wb.imag + wc.imag) / 3 * SCALE
    dist = np.sqrt((wa.real+wb.real+wc.real)**2 + (wa.imag+wb.imag+wc.imag)**2) / 3
    if dist > 0.7: continue
    label = fund_labels[ti]
    lines.append(f"  \\node[white, opacity=0.9, font=\\fontsize{{4}}{{4}}\\selectfont] at ({cx:.4f},{cy:.4f}) {{{label}}};")

# Disk border on top
lines.append("")
lines.append(f"  \\draw[gray!50, line width=1pt] (0,0) circle ({SCALE}cm);")

# Title
lines.append(f"  \\node[white, anchor=south, font=\\sffamily\\bfseries] at (0,{SCALE+0.4}) "
             f"{{Bolza Surface -- Poincar\\'e Disk}};")
lines.append(f"  \\node[gray!70, anchor=north, font=\\sffamily\\small] at (0,{-SCALE-0.3}) "
             f"{{(2,3,8) tiling $\\cdot$ GL(2,3) $\\cdot$ genus 2 $\\cdot$ $\\chi=-2$ $\\cdot$ "
             f"$H^1 = V_2 \\oplus \\bar V_2$}};")

lines.append(r"\end{tikzpicture}")
lines.append(r"\end{document}")

tex_path = '/Users/robin/git/graphviz/bolza-surface/bolza_disk.tex'
with open(tex_path, 'w') as f:
    f.write('\n'.join(lines))

print(f"  {len(lines)} lines of TikZ")
print(f"  → {tex_path}")
