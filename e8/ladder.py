"""
Shared machinery for the exceptional ladder E8 > E7 > E6 (minuscule descent).

This module makes the whole pipeline reproducible without the session
transcript: e8_mandala.py, e7_build.py and e6_build.py all import from here.
Nothing is cached here; the build scripts own the /tmp/*.npy files that the
drawers (e8_one.py, e7_one.py, e6_one.py, ...) read.

Everything is direct linear algebra in R^8:
  - the 240 E8 roots (length^2 = 2),
  - the Weyl-subsystem orthogonal to one/two fixed vectors (E7 = roots _|_ r0,
    E6 = roots _|_ {r0, w0}),
  - a Coxeter element of that subsystem and its Coxeter-plane eigenprojection.

The descent vectors are r0 = e1+e2 and w0 = e1+e3 (recovered from the original
session's caches; any pair in the same W(E8)-orbit gives an isomorphic ladder).
"""
import itertools
import numpy as np

R0 = np.array([1, 1, 0, 0, 0, 0, 0, 0], float)   # fixes E7 inside E8
W0 = np.array([1, 0, 1, 0, 0, 0, 0, 0], float)   # fixes E6 inside E7


def e8_roots():
    """The 240 roots of E8 in the standard R^8 coordinates, each length^2 = 2."""
    roots = []
    # 112 of type (+-1, +-1, 0^6)
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0]*8
                v[i] = si
                v[j] = sj
                roots.append(v)
    # 128 of type (+-1/2)^8 with an even number of minus signs
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            roots.append([s*0.5 for s in signs])
    R = np.array(roots, float)
    assert len(R) == 240 and np.allclose((R*R).sum(1), 2.0)
    return R


def simple_roots(S):
    """Simple roots of the root system whose rows are S.

    Pick a functional generic enough to split the roots into positive/negative
    with no ties, then keep the positive roots that are not a sum of two
    positive roots.  The irrational weights are what guarantee genericity
    (a rational functional can tie two roots and over-count)."""
    rank = np.linalg.matrix_rank(S)
    v = np.sqrt([2, 3, 5, 7, 11, 13, 17, 19.])
    c = S.T @ (S @ v)                        # a covector lying in span(S)
    height = np.round(S @ c, 9)
    pos = S[height > 0]
    posset = set(map(tuple, np.round(pos, 6)))
    simple = [r for r in pos
              if not any(tuple(np.round(r - s, 6)) in posset
                         for s in pos if not np.allclose(s, r))]
    simple = np.array(simple)
    assert len(simple) == rank, f"found {len(simple)} simple roots, rank {rank}"
    return simple


def coxeter_plane(simple, h):
    """Return (order, u1, u2): the order of the Coxeter element built from
    these simple reflections (must equal h) and an orthonormal basis of the
    Coxeter plane (the e^{2 pi i / h} eigenplane)."""
    def refl(a):
        return np.eye(8) - 2.0 * np.outer(a, a) / (a @ a)

    C = np.eye(8)
    for a in simple:
        C = C @ refl(a)
    order = next(k for k in range(1, 200)
                 if np.allclose(np.linalg.matrix_power(C, k), np.eye(8), atol=1e-7))
    assert order == h, f"Coxeter order {order} != {h}"
    w, V = np.linalg.eig(C)
    idx = np.argmin(np.abs(np.angle(w) - 2*np.pi/h))
    v = V[:, idx]
    u1 = np.real(v); u1 /= np.linalg.norm(u1)
    u2 = np.imag(v); u2 -= (u2 @ u1) * u1; u2 /= np.linalg.norm(u2)
    return order, u1, u2


def project(vectors, u1, u2):
    """Coxeter-plane coordinates of each vector."""
    return np.column_stack([vectors @ u1, vectors @ u2])


def ring_histogram(P):
    """radius -> count, the concentric-ring signature of a Coxeter projection."""
    from collections import Counter
    radii = np.round(np.hypot(P[:, 0], P[:, 1]), 4)
    return dict(sorted(Counter(radii).items()))
