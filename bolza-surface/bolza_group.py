"""
Shared infrastructure for Bolza surface computations.

GL(2, F_3) of order 48, SL(2,R) generators for the (2,3,8) triangle group,
Poincare disk tiling, and Dirichlet fundamental domain.

Extracted from build_viz.py / fold_octagon.py / bolza_cohomology.py / build_tikz.py
to eliminate quadruple duplication.
"""

import numpy as np
from itertools import product
from collections import deque


# ══════════════════════════════════════════════════════════════════════
#  GL(2, F_3) — the automorphism group of the Bolza surface
# ══════════════════════════════════════════════════════════════════════

def _build_gl23():
    """Build GL(2, F_3) as a list of 48 2x2 matrices over F_3."""
    def mm3(A, B): return (A @ B) % 3
    def mk3(A): return tuple((A % 3).flatten())
    def md3(A): return (A[0,0]*A[1,1] - A[0,1]*A[1,0]) % 3
    def mi3(A):
        d = md3(A); di = 1 if d == 1 else 2
        return (di * np.array([[A[1,1], (-A[0,1])%3], [(-A[1,0])%3, A[0,0]]])) % 3

    el = {}
    for a, b, c, d in product(range(3), repeat=4):
        A = np.array([[a, b], [c, d]])
        if md3(A) == 0: continue
        k = mk3(A)
        if k not in el: el[k] = A % 3

    matrices = list(el.values())
    keys = list(el.keys())
    key_to_idx = {k: i for i, k in enumerate(keys)}
    assert len(matrices) == 48

    def mult(i, j):
        return key_to_idx[mk3(mm3(matrices[i], matrices[j]))]
    def inv(i):
        return key_to_idx[mk3(mi3(matrices[i]))]

    identity = key_to_idx[mk3(np.eye(2, dtype=int))]

    orders = []
    for i in range(48):
        c = i
        for n in range(1, 49):
            if c == identity:
                orders.append(n); break
            c = mult(c, i)

    return matrices, mult, inv, identity, orders


def _find_238_generators(mult, identity, orders):
    """Find generators a (order 2), b (order 3) with ab of order 8."""
    for a in range(48):
        if orders[a] != 2: continue
        for b in range(48):
            if orders[b] != 3: continue
            if orders[mult(a, b)] == 8:
                return a, b, mult(a, b)
    raise RuntimeError("No (2,3,8) generators found")


class GL23:
    """GL(2, F_3) of order 48, with (2,3,8) presentation."""

    def __init__(self):
        self.matrices, self._mult, self._inv, self.identity, self.orders = _build_gl23()
        self.a_gen, self.b_gen, self.c_gen = _find_238_generators(
            self._mult, self.identity, self.orders)
        self.b_inv_gen = self._inv(self.b_gen)

    def mult(self, i, j):
        return self._mult(i, j)

    def inv(self, i):
        return self._inv(i)

    def order(self, i):
        return self.orders[i]

    def word_to_element(self, word):
        """Evaluate a word in {a, b, B} to a GL(2,3) element index."""
        gen_map = {'a': self.a_gen, 'b': self.b_gen, 'B': self.b_inv_gen}
        g = self.identity
        for letter in word:
            g = self._mult(g, gen_map[letter])
        return g


# ══════════════════════════════════════════════════════════════════════
#  SL(2,R) generators for the (2,3,8) triangle group
# ══════════════════════════════════════════════════════════════════════

def _build_sl2r_generators():
    """Construct SL(2,R) matrices for the (2,3,8) triangle group.

    Returns a_real, b_real, c_real = a@b, b_inv_real, and the
    three fixed points z_a, z_b, z_c in the upper half-plane.
    """
    a_real = np.array([[0, 1], [-1, 0]], dtype=float)

    cos_pi8 = np.cos(np.pi / 8)
    sin_pi3 = np.sin(np.pi / 3)
    cos_pi3 = np.cos(np.pi / 3)

    # sin(pi/3)(f + 1/f) = 2*cos(pi/8)
    rhs = 2 * cos_pi8 / sin_pi3
    f_param = (rhs + np.sqrt(rhs**2 - 4)) / 2

    b_real = np.array([
        [cos_pi3, -f_param * sin_pi3],
        [sin_pi3 / f_param, cos_pi3]
    ], dtype=float)

    b_inv_real = np.linalg.inv(b_real)
    c_real = a_real @ b_real

    def fixed_point(M):
        a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
        disc = 4 - (a + d)**2
        if disc < 0: return None
        z1 = ((a - d) + 1j * np.sqrt(disc)) / (2 * c)
        z2 = ((a - d) - 1j * np.sqrt(disc)) / (2 * c)
        return z1 if z1.imag > 0 else z2

    z_a = fixed_point(a_real)
    z_b = fixed_point(b_real)
    z_c = fixed_point(c_real)

    return a_real, b_real, c_real, b_inv_real, z_a, z_b, z_c


class TriangleGroup:
    """Delta(2,3,8) as SL(2,R) matrices acting on H and the Poincare disk."""

    def __init__(self):
        (self.a_real, self.b_real, self.c_real, self.b_inv_real,
         self.z_a, self.z_b, self.z_c) = _build_sl2r_generators()

        self.gen_mats = {
            'a': self.a_real, 'b': self.b_real, 'B': self.b_inv_real
        }

        # Disk coordinates of the fixed points
        self.w_a = self.to_disk(self.z_a)
        self.w_b = self.to_disk(self.z_b)
        self.w_c = self.to_disk(self.z_c)  # order-8 vertex, recentered to origin

        # Recentering constant
        self._w_c_raw = self.w_c  # before recentering, this IS w_c

    @staticmethod
    def mobius(M, z):
        """Apply SL(2,R) matrix M as Mobius transformation on H."""
        a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
        return (a * z + b) / (c * z + d)

    @staticmethod
    def to_disk(z):
        """Cayley map: upper half-plane -> Poincare disk. Sends i -> 0."""
        return (z - 1j) / (z + 1j)

    @staticmethod
    def to_halfplane(w):
        """Inverse Cayley map: Poincare disk -> upper half-plane."""
        return 1j * (1 + w) / (1 - w)

    def recenter(self, w):
        """Mobius disk automorphism mapping the order-8 vertex to origin."""
        return (w - self._w_c_raw) / (1 - np.conj(self._w_c_raw) * w)

    def recenter_inv(self, w):
        """Inverse: origin -> order-8 vertex."""
        return (w + self._w_c_raw) / (1 + np.conj(self._w_c_raw) * w)

    def disk_to_halfplane(self, w_recentered):
        """Recentered disk coordinate -> upper half-plane coordinate."""
        w_raw = self.recenter_inv(w_recentered)
        return self.to_halfplane(w_raw)

    @staticmethod
    def hyp_dist(w):
        """Hyperbolic distance from origin in the Poincare disk."""
        r = abs(w)
        if r >= 1: return 1e9
        return 2 * np.arctanh(r)

    @staticmethod
    def mat_key(M, tol=1e-6):
        """Round matrix entries for de-duplication."""
        return tuple(round(x / tol) * tol for x in M.flatten())


# ══════════════════════════════════════════════════════════════════════
#  BFS tiling + Dirichlet fundamental domain
# ══════════════════════════════════════════════════════════════════════

class BolzaTiling:
    """BFS tiling of the Poincare disk, with 48-triangle Dirichlet domain."""

    def __init__(self, gl23: GL23, tg: TriangleGroup,
                 max_triangles=4000, disk_cutoff=0.996):
        self.gl23 = gl23
        self.tg = tg

        # BFS to generate triangles
        all_tris, all_mats, all_words = self._bfs(tg, max_triangles, disk_cutoff)

        # Select 48 closest as fundamental domain
        centers = [(t[0] + t[1] + t[2]) / 3 for t in all_tris]
        dists = [(i, tg.hyp_dist(centers[i])) for i in range(len(all_tris))]
        dists.sort(key=lambda x: x[1])

        self.fund_triangles = []  # 48 triangles as (wa, wb, wc) in recentered disk
        self.fund_mats = []       # SL(2,R) matrices that produced each triangle
        self.fund_words = []      # words in {a, b, B}
        self.fund_gl23 = []       # GL(2,3) element index for each triangle
        for idx, _ in dists[:48]:
            self.fund_triangles.append(all_tris[idx])
            self.fund_mats.append(all_mats[idx])
            self.fund_words.append(all_words[idx])
            self.fund_gl23.append(gl23.word_to_element(all_words[idx]))

        self.all_triangles = [all_tris[idx] for idx, _ in dists]

        # De-duplicate vertices
        self.vertices, self.tri_indices = self._index_vertices(self.fund_triangles)

    @staticmethod
    def _bfs(tg, max_triangles, disk_cutoff):
        gen_list = [('a', tg.a_real), ('b', tg.b_real), ('B', tg.b_inv_real)]
        all_tris = []
        all_mats = []
        all_words = []
        visited = set()
        queue = deque([(np.eye(2), "")])

        while queue and len(all_tris) < max_triangles:
            M, word = queue.popleft()
            mk = tg.mat_key(M)
            if mk in visited: continue
            visited.add(mk)

            try:
                wa = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_a)))
                wb = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_b)))
                wc = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_c)))
            except (ZeroDivisionError, FloatingPointError):
                continue

            if max(abs(wa), abs(wb), abs(wc)) > disk_cutoff:
                continue
            area = abs((wb-wa).real*(wc-wa).imag - (wb-wa).imag*(wc-wa).real)
            if area < 1e-10: continue

            all_tris.append((wa, wb, wc))
            all_mats.append(M.copy())
            all_words.append(word)

            for gn, gmat in gen_list:
                queue.append((M @ gmat, word + gn))

        return all_tris, all_mats, all_words

    @staticmethod
    def _index_vertices(triangles, tol=1e-8):
        """De-duplicate vertices, return (vertex_list, triangle_index_list)."""
        vmap = {}
        vlist = []
        def vidx(z):
            k = (round(z.real / tol) * tol, round(z.imag / tol) * tol)
            if k not in vmap:
                vmap[k] = len(vlist)
                vlist.append(z)
            return vmap[k]

        tri_indices = []
        for wa, wb, wc in triangles:
            tri_indices.append((vidx(wa), vidx(wb), vidx(wc)))

        return vlist, tri_indices

    def octagon_vertices(self):
        """Find the 8 octagon vertices (outermost points of the fundamental domain)."""
        R_oct = max(abs(v) for v in self.vertices)
        ov = sorted(
            [i for i, v in enumerate(self.vertices) if abs(abs(v) - R_oct) < 0.01],
            key=lambda i: np.angle(self.vertices[i])
        )
        return ov

    def weierstrass_points_disk(self):
        """Find the 6 Weierstrass points in the recentered Poincare disk.

        These are the 6 order-8 vertices of the {8,3} cell complex,
        i.e., the orbits of the order-8 fixed point under GL(2,3)/Stab.
        """
        gl23 = self.gl23
        tg = self.tg

        # The stabilizer of the order-8 vertex is <c> of order 8.
        # Coset representatives give the 6 = 48/8 Weierstrass points.
        c_subgroup = set()
        g = gl23.identity
        for _ in range(8):
            c_subgroup.add(g)
            g = gl23.mult(g, gl23.c_gen)

        covered = set()
        coset_reps = []
        for g in range(48):
            if g in covered: continue
            coset = frozenset(gl23.mult(g, h) for h in c_subgroup)
            coset_reps.append(g)
            covered |= coset

        # For each coset rep, find the SL(2,R) matrix and compute the
        # image of z_c (the order-8 fixed point)
        weierstrass = []
        for g in coset_reps:
            word = self._element_to_word(g)
            M = np.eye(2)
            for letter in word:
                M = M @ tg.gen_mats[letter]
            w = tg.recenter(tg.to_disk(tg.mobius(M, tg.z_c)))
            weierstrass.append(w)

        return weierstrass

    def _element_to_word(self, g):
        """Convert GL(2,3) element index to a word in {a, b, B} via BFS."""
        gl23 = self.gl23
        words = {gl23.identity: ""}
        queue = deque([gl23.identity])
        gens = [('a', gl23.a_gen), ('b', gl23.b_gen), ('B', gl23.b_inv_gen)]
        while queue:
            cur = queue.popleft()
            if cur == g:
                return words[g]
            for gn, gi in gens:
                h = gl23.mult(cur, gi)
                if h not in words:
                    words[h] = words[cur] + gn
                    queue.append(h)
        return words.get(g, "")
