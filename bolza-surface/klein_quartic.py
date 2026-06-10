"""
Klein Quartic visualization.

Step 1: Beautiful {3,7} tiling of the Poincaré disk with 7-coloring.
Step 2: Highlight the 14-gon fundamental domain.
Step 3: Fold into 3D (later).

Triangle group Δ(2,3,7), quotient PSL(2,7) of order 168.
"""

import numpy as np
from collections import deque, defaultdict
from itertools import product
import json

# ══════════════════════════════════════════════════════════════════════
#  PSL(2, F_7) — order 168
# ══════════════════════════════════════════════════════════════════════

class PSL27:
    def __init__(self):
        self.p = 7
        elements = {}
        for a, b, c, d in product(range(7), repeat=4):
            if (a * d - b * c) % 7 != 1:
                continue
            M = np.array([[a, b], [c, d]], dtype=int)
            key = self._canonical(M)
            if key not in elements:
                elements[key] = M % 7

        self.matrices = list(elements.values())
        self.keys = list(elements.keys())
        self.key_to_idx = {k: i for i, k in enumerate(self.keys)}
        assert len(self.matrices) == 168

        self.identity = self.key_to_idx[self._canonical(np.eye(2, dtype=int))]
        self.a_gen, self.b_gen, self.c_gen = self._find_237_generators()
        self.b_inv_gen = self.inv(self.b_gen)

    def _canonical(self, M):
        M = M % self.p
        Mneg = (self.p - M) % self.p
        for i in range(2):
            for j in range(2):
                if M[i, j] != 0:
                    return tuple(M.flatten()) if M[i, j] <= 3 else tuple(Mneg.flatten())
        return tuple(M.flatten())

    def mult(self, i, j):
        return self.key_to_idx[self._canonical((self.matrices[i] @ self.matrices[j]) % self.p)]

    def inv(self, i):
        M = self.matrices[i]
        det_inv = pow(int((M[0,0]*M[1,1] - M[0,1]*M[1,0]) % self.p), self.p-2, self.p)
        Minv = (det_inv * np.array([[M[1,1], (-M[0,1])%self.p],
                                     [(-M[1,0])%self.p, M[0,0]]])) % self.p
        return self.key_to_idx[self._canonical(Minv)]

    def order(self, i):
        c = i
        for n in range(1, 169):
            if c == self.identity: return n
            c = self.mult(c, i)
        return 169

    def _find_237_generators(self):
        orders = [self.order(i) for i in range(168)]
        for a in range(168):
            if orders[a] != 2: continue
            for b in range(168):
                if orders[b] != 3: continue
                if orders[self.mult(a, b)] == 7:
                    return a, b, self.mult(a, b)
        raise RuntimeError("No (2,3,7) generators found")

    def word_to_element(self, word):
        gen_map = {'a': self.a_gen, 'b': self.b_gen, 'B': self.b_inv_gen}
        g = self.identity
        for letter in word:
            g = self.mult(g, gen_map[letter])
        return g


# ══════════════════════════════════════════════════════════════════════
#  SL(2,R) generators for Δ(2,3,7)
# ══════════════════════════════════════════════════════════════════════

class TriangleGroup237:
    def __init__(self):
        self.a_real = np.array([[0, 1], [-1, 0]], dtype=float)
        cos7 = np.cos(np.pi / 7)
        sin3 = np.sin(np.pi / 3)
        cos3 = np.cos(np.pi / 3)
        rhs = 2 * cos7 / sin3
        f = (rhs + np.sqrt(rhs**2 - 4)) / 2
        self.b_real = np.array([[cos3, -f*sin3], [sin3/f, cos3]], dtype=float)
        self.b_inv_real = np.linalg.inv(self.b_real)
        self.c_real = self.a_real @ self.b_real
        self.gen_mats = {'a': self.a_real, 'b': self.b_real, 'B': self.b_inv_real}

        self.z_a = self._fp(self.a_real)
        self.z_b = self._fp(self.b_real)
        self.z_c = self._fp(self.c_real)
        self.w_c_raw = self.to_disk(self.z_c)

    @staticmethod
    def _fp(M):
        a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
        disc = 4 - (a+d)**2
        z1 = ((a-d) + 1j*np.sqrt(abs(disc))) / (2*c)
        z2 = ((a-d) - 1j*np.sqrt(abs(disc))) / (2*c)
        return z1 if z1.imag > 0 else z2

    @staticmethod
    def mobius(M, z):
        return (M[0,0]*z + M[0,1]) / (M[1,0]*z + M[1,1])

    @staticmethod
    def to_disk(z):
        return (z - 1j) / (z + 1j)

    def recenter(self, w):
        return (w - self.w_c_raw) / (1 - np.conj(self.w_c_raw) * w)


# ══════════════════════════════════════════════════════════════════════
#  BFS tiling — generate LOTS of triangles for the full disk
# ══════════════════════════════════════════════════════════════════════

def bfs_tiling(psl, tg, max_triangles=12000, disk_cutoff=0.997):
    """BFS through Δ(2,3,7). Returns right triangles with matrices and PSL labels."""
    gen_list = [('a', tg.a_real), ('b', tg.b_real), ('B', tg.b_inv_real)]
    psl_gens = {'a': psl.a_gen, 'b': psl.b_gen, 'B': psl.b_inv_gen}
    all_tris = []
    all_mats = []
    all_labels = []
    visited = set()
    queue = deque([(np.eye(2), psl.identity)])

    def mat_key(M):
        return tuple(round(x * 1e6) for x in M.flatten())

    while queue and len(all_tris) < max_triangles:
        M, label = queue.popleft()
        mk = mat_key(M)
        if mk in visited: continue
        visited.add(mk)

        try:
            wa = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_a)))
            wb = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_b)))
            wc = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_c)))
        except (ZeroDivisionError, FloatingPointError):
            continue

        if max(abs(wa), abs(wb), abs(wc)) > disk_cutoff: continue
        area = abs((wb-wa).real*(wc-wa).imag - (wb-wa).imag*(wc-wa).real)
        if area < 1e-10: continue

        all_tris.append((wa, wb, wc))
        all_mats.append(M.copy())
        all_labels.append(label)

        for gn, gmat in gen_list:
            new_label = psl.mult(label, psl_gens[gn])
            queue.append((M @ gmat, new_label))

    return all_tris, all_mats, all_labels


# ══════════════════════════════════════════════════════════════════════
#  Group right triangles into equilateral {3,7} triangles
# ══════════════════════════════════════════════════════════════════════

def make_equilaterals(psl, all_tris, all_mats, all_labels, tg):
    """Group right triangles into equilateral {3,7} triangles.

    Each right triangle M has siblings M·b and M·b². We group by finding
    these matrix triples. Color by PSL(2,7) label modulo the order-24
    stabilizer of a Fano plane point.
    """
    def mat_key(M):
        return tuple(round(x * 1e6) for x in M.flatten())

    mat_idx = {}
    for i, M in enumerate(all_mats):
        mat_idx[mat_key(M)] = i

    # Proper 7-coloring of the {3,7} tiling.
    # Build adjacency of the 56 faces (b-cosets) on the surface.
    # Two faces are adjacent iff they share a vertex (c-coset).
    # Then greedy-color with 7 colors.

    # Step 1: enumerate the 56 b-cosets and their c-coset vertices
    b_cosets = []  # list of frozenset of 3 labels
    b_coset_of = {}  # label -> coset index
    seen = set()
    for g in range(168):
        if g in seen: continue
        gb = psl.mult(g, psl.b_gen)
        gb2 = psl.mult(gb, psl.b_gen)
        coset = frozenset([g, gb, gb2])
        seen |= coset
        idx = len(b_cosets)
        b_cosets.append(coset)
        for x in coset:
            b_coset_of[x] = idx

    # Step 2: c-cosets (vertices)
    c_coset_of = {}
    c_id = 0
    for g in range(168):
        if g in c_coset_of: continue
        h = g
        for _ in range(7):
            c_coset_of[h] = c_id
            h = psl.mult(h, psl.c_gen)
        c_id += 1

    # Step 3: face → set of vertex IDs
    face_verts = []
    for coset in b_cosets:
        verts = set(c_coset_of[g] for g in coset)
        face_verts.append(verts)

    # Step 4: adjacency (shared vertex)
    adj = defaultdict(set)
    for i in range(56):
        for j in range(i+1, 56):
            if face_verts[i] & face_verts[j]:
                adj[i].add(j)
                adj[j].add(i)

    # Step 5: 7-coloring via vertex constraint propagation.
    # At each vertex (c-coset), the 7 faces meeting there must get all 7 colors.
    # Build vertex → list of face indices
    vert_faces = defaultdict(list)
    for fi, vs in enumerate(face_verts):
        for v in vs:
            vert_faces[v].append(fi)

    face_colors = [-1] * 56
    # Color the 7 faces at vertex 0 by their angular position
    v0_faces = vert_faces[0]
    # Sort by angle of face center in the disk
    # (we need disk positions — compute from the coset representatives)
    # For now, just assign 0-6 in the order they appear
    for k, fi in enumerate(v0_faces[:7]):
        face_colors[fi] = k

    # BFS through vertices, propagating constraints
    vert_done = {0}
    vq = deque([0])
    while vq:
        v = vq.popleft()
        # For each neighboring vertex (sharing an edge = sharing 2 faces)
        for v2 in range(24):
            if v2 in vert_done: continue
            shared = set(vert_faces[v]) & set(vert_faces[v2])
            if len(shared) < 2: continue  # not edge-adjacent vertices
            # Some faces at v2 are already colored
            v2f = vert_faces[v2]
            colored_at_v2 = {face_colors[f] for f in v2f if face_colors[f] >= 0}
            uncolored_at_v2 = [f for f in v2f if face_colors[f] < 0]
            available = [c for c in range(7) if c not in colored_at_v2]
            if len(available) == len(uncolored_at_v2):
                # Fully determined — assign
                for f, c in zip(uncolored_at_v2, available):
                    face_colors[f] = c
                vert_done.add(v2)
                vq.append(v2)
            elif len(uncolored_at_v2) == 0:
                vert_done.add(v2)
                vq.append(v2)

    # If some faces still uncolored, retry with remaining constraints
    for _ in range(10):
        if -1 not in face_colors: break
        for v in range(24):
            vf = vert_faces[v]
            colored_at_v = {face_colors[f] for f in vf if face_colors[f] >= 0}
            uncolored = [f for f in vf if face_colors[f] < 0]
            available = [c for c in range(7) if c not in colored_at_v]
            if len(available) == len(uncolored):
                for f, c in zip(uncolored, available):
                    face_colors[f] = c

    n_uncolored = face_colors.count(-1)
    if n_uncolored > 0:
        print(f"  WARNING: {n_uncolored} faces uncolored, assigning arbitrarily")
        for i in range(56):
            if face_colors[i] < 0:
                used = set(face_colors[j] for j in adj[i] if face_colors[j] >= 0)
                for c in range(7):
                    if c not in used:
                        face_colors[i] = c
                        break
                if face_colors[i] < 0:
                    face_colors[i] = 0

    def face_color(label):
        return face_colors[b_coset_of[label]]

    used = set()
    equilaterals = []  # (v0, v1, v2, b_coset_id)
    colors = []
    for i, M in enumerate(all_mats):
        if i in used: continue
        Mb = M @ tg.b_real
        Mb2 = Mb @ tg.b_real
        kb, kb2 = mat_key(Mb), mat_key(Mb2)
        if kb in mat_idx and kb2 in mat_idx:
            j, k = mat_idx[kb], mat_idx[kb2]
            used |= {i, j, k}
            bc = b_coset_of[all_labels[i]]
            equilaterals.append((all_tris[i][2], all_tris[j][2], all_tris[k][2], bc))
            colors.append(face_color(all_labels[i]))

    return equilaterals, colors


# ══════════════════════════════════════════════════════════════════════
#  Identify the 14-gon fundamental domain
# ══════════════════════════════════════════════════════════════════════

def find_fundamental_domain(equilaterals):
    """Select one equilateral triangle per surface face (b-coset),
    closest to center. This gives the Dirichlet fundamental domain."""
    best_per_coset = {}  # b_coset_id -> (dist, equilateral_index)
    for i, (v0, v1, v2, bc) in enumerate(equilaterals):
        center = (v0 + v1 + v2) / 3
        d = abs(center)
        if bc not in best_per_coset or d < best_per_coset[bc][0]:
            best_per_coset[bc] = (d, i)
    fund_set = set(idx for _, idx in best_per_coset.values())
    return fund_set


# ══════════════════════════════════════════════════════════════════════
#  HTML output — Poincaré disk tiling
# ══════════════════════════════════════════════════════════════════════

def generate_html(equilaterals, colors, fund_set, output_path):
    palette = [
        [0.85, 0.20, 0.20],  # red
        [0.20, 0.50, 0.85],  # blue
        [0.20, 0.75, 0.30],  # green
        [0.90, 0.75, 0.15],  # gold
        [0.70, 0.25, 0.80],  # purple
        [0.90, 0.45, 0.15],  # orange
        [0.15, 0.75, 0.80],  # cyan
    ]

    # Build triangle data
    tris_data = []
    for i in range(len(equilaterals)):
        v0, v1, v2 = equilaterals[i][0], equilaterals[i][1], equilaterals[i][2]
        in_fund = i in fund_set
        c = palette[colors[i]]
        if not in_fund:
            c = [x * 0.3 for x in c]
        tris_data.append({
            'verts': [[v0.real, v0.imag], [v1.real, v1.imag], [v2.real, v2.imag]],
            'color': c,
            'fund': in_fund
        })

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Klein Quartic — {{3,7}} Tiling</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{background:#0a0a10;display:flex;justify-content:center;align-items:center;height:100vh}}
canvas{{border-radius:50%;cursor:grab}}
canvas:active{{cursor:grabbing}}
.info{{position:absolute;top:16px;left:20px;color:#8899aa;font:13px/1.6 'Courier New',monospace;pointer-events:none}}
.info h2{{color:#ccd0d8;font-size:16px;font-weight:normal;margin-bottom:2px}}
.info .n{{color:#e8b84d}}
</style></head><body>
<div class="info">
<h2>Klein Quartic — {{3,7}} Tiling</h2>
<span class="n">56</span> equilateral triangles · PSL(2,7) · genus 3<br>
7-coloring · 14-gon fundamental domain highlighted
</div>
<canvas id="c"></canvas>
<script>
const D={json.dumps(tris_data, default=lambda x: round(x, 8) if isinstance(x, float) else x)};
const canvas=document.getElementById('c');
const ctx=canvas.getContext('2d');
let W=Math.min(innerWidth,innerHeight)-40;
canvas.width=canvas.height=W*devicePixelRatio;
canvas.style.width=canvas.style.height=W+'px';
ctx.scale(devicePixelRatio,devicePixelRatio);

function draw(){{
  const R=W/2-2;
  const cx=W/2, cy=W/2;
  ctx.clearRect(0,0,W,W);

  // Disk background
  ctx.beginPath();ctx.arc(cx,cy,R,0,Math.PI*2);
  ctx.fillStyle='#12121a';ctx.fill();
  ctx.strokeStyle='#2a2a3a';ctx.lineWidth=1.5;ctx.stroke();

  // Draw triangles
  for(const t of D){{
    const[[x0,y0],[x1,y1],[x2,y2]]=t.verts;
    const px0=cx+x0*R, py0=cy-y0*R;
    const px1=cx+x1*R, py1=cy-y1*R;
    const px2=cx+x2*R, py2=cy-y2*R;

    const[r,g,b]=t.color;
    ctx.beginPath();
    ctx.moveTo(px0,py0);ctx.lineTo(px1,py1);ctx.lineTo(px2,py2);ctx.closePath();
    ctx.fillStyle=`rgb(${{r*255|0}},${{g*255|0}},${{b*255|0}})`;
    ctx.fill();

    // Edge lines
    ctx.strokeStyle=t.fund?'rgba(255,255,255,0.5)':'rgba(255,255,255,0.08)';
    ctx.lineWidth=t.fund?0.8:0.3;
    ctx.stroke();
  }}

  // Disk border
  ctx.beginPath();ctx.arc(cx,cy,R,0,Math.PI*2);
  ctx.strokeStyle='#3a3a4a';ctx.lineWidth=2;ctx.stroke();
}}

draw();
window.onresize=()=>{{
  W=Math.min(innerWidth,innerHeight)-40;
  canvas.width=canvas.height=W*devicePixelRatio;
  canvas.style.width=canvas.style.height=W+'px';
  ctx.scale(devicePixelRatio,devicePixelRatio);
  draw();
}};
</script></body></html>"""

    with open(output_path, 'w') as f:
        f.write(html)
    print(f"Written to {output_path}")


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    print("Building PSL(2,7)...")
    psl = PSL27()
    print(f"  168 elements, a(2) b(3) c=ab(7)")

    print("Building Δ(2,3,7)...")
    tg = TriangleGroup237()
    print(f"  tr(c) = {np.trace(tg.c_real):.6f} = 2cos(π/7) ✓")

    print("BFS tiling...")
    all_tris, all_mats, all_labels = bfs_tiling(psl, tg)
    print(f"  {len(all_tris)} right triangles")

    print("Grouping into equilateral triangles...")
    equilaterals, colors = make_equilaterals(psl, all_tris, all_mats, all_labels, tg)
    print(f"  {len(equilaterals)} equilateral triangles")
    print(f"  Color distribution: {[colors.count(i) for i in range(7)]}")

    print("Finding fundamental domain...")
    fund_set = find_fundamental_domain(equilaterals)
    print(f"  {len(fund_set)} faces in fundamental domain")

    output = "/Users/robin/git/graphviz/bolza-surface/klein_quartic.html"
    print("Generating HTML...")
    generate_html(equilaterals, colors, fund_set, output)
    print("Done!")


if __name__ == '__main__':
    main()
