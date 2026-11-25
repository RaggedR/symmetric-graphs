"""
The exceptional ladder as a tower of ANTIPODAL DOUBLE COVERS.
Resolves open thread #1 of miriam-handover.md (and corrects its conjecture).

The conjecture was: the ladder folds into a tower of named strongly-regular
graphs, with E7 folding to SRG(28,12,6,4).  That is FALSE.  What is true, and
cleaner, is:

  * each +1-graph is an antipodal distance-regular graph of diameter 3;
  * its antipodal quotient (fold) is the sigma-even half of its spectrum;
  * E8 folds to a famous SRG, but E7 folds only to the complete graph K28,
    and E6 (odd order 27) has no antipodal pairs at all -- it is terminal.

  rung | +1-graph (deg)      | vertices | antipodal quotient
  -----+---------------------+----------+------------------------------------
  E8   | 4_21 / Gosset (56)  |   240    | SRG(120,56,28,24)  [O+(8,2):2]
  E7   | Gosset graph  (27)  |    56    | K28  (complete -- NOT SRG(28,12,6,4))
  E6   | Schlafli SRG  (10)  |    27    | terminal (odd order, no fold)

Spectral punchline: folding keeps the sigma-even eigenspaces and kills the
sigma-odd ones; the rank-n reflection representation is exactly the sigma-odd
part, sitting at one eigenvalue with multiplicity 8 -> 7 -> 6 (the rank descent
E8 > E7 > E6).

Reads /tmp/{e8_R,e7_W,e6_W}.npy (run e8_mandala.py, e7_build.py, e6_build.py).
Pure linear algebra; nothing is drawn.
"""
import collections
import numpy as np


def spectrum(A):
    w = np.linalg.eigvalsh(A.astype(float))
    return dict(sorted(collections.Counter(np.round(w, 4)).items(), reverse=True))


def srg_params(A):
    """(n,k,lambda,mu) if A is a strongly regular graph, else None."""
    n = A.shape[0]
    if len(set(A.sum(1))) != 1:
        return None
    A2 = A @ A
    lam = {int(A2[i, j]) for i in range(n) for j in range(n) if i != j and A[i, j]}
    mu = {int(A2[i, j]) for i in range(n) for j in range(n) if i != j and not A[i, j]}
    if len(lam) == 1 and len(mu) == 1:
        return (n, int(A.sum(1)[0]), lam.pop(), mu.pop())
    return None


def bfs_distances(A):
    n = A.shape[0]
    D = np.full((n, n), -1)
    for s in range(n):
        D[s, s] = 0
        q = collections.deque([s])
        while q:
            u = q.popleft()
            for v in np.nonzero(A[u])[0]:
                if D[s, v] < 0:
                    D[s, v] = D[s, u] + 1
                    q.append(v)
    return D


def antipode_split(A, sigma):
    """Spectra of A restricted to the sigma-even and sigma-odd subspaces."""
    n = A.shape[0]
    Psig = np.zeros((n, n))
    Psig[np.arange(n), sigma] = 1.0
    val, vec = np.linalg.eigh(Psig)
    even, odd = vec[:, val > 0.5], vec[:, val < -0.5]
    return spectrum(even.T @ A @ even), spectrum(odd.T @ A @ odd)


def analyse(name, V, adj_ip, group_label):
    G = np.round(V @ V.T).astype(int)
    n = len(V)
    A = (G == adj_ip).astype(int)
    D = bfs_distances(A)
    diam = int(D.max())
    dist = dict(sorted(collections.Counter(D[0]).items()))
    print(f"\n{'='*66}\n{name}  --  {n} vertices, +1-graph degree {A.sum(1)[0]}\n{'='*66}")
    print(f"  diameter {diam}, distance distribution from a vertex: {dist}")
    print(f"  +1-graph spectrum: {spectrum(A)}")

    far = (D == diam).astype(int)
    if far.sum(1)[0] != 1:
        print("  no antipodal matching (far-degree != 1) -> TERMINAL rung.")
        print(f"  it is already {('SRG'+str(srg_params(A))) if srg_params(A) else 'not SRG'}.")
        return
    sigma = far.argmax(1)
    m = n // 2
    print(f"  antipodal: distance-{diam} relation is a perfect matching -> "
          f"double cover of a {m}-vertex graph")

    # the fold (antipodal quotient): merge each antipodal pair
    reps = [i for i in range(n) if i < sigma[i]]
    Q = np.zeros((m, m), int)
    for a, i in enumerate(reps):
        for b, j in enumerate(reps):
            if a != b and (A[i, j] or A[i, sigma[j]]):
                Q[a, b] = 1
    p = srg_params(Q)
    if p:
        print(f"  fold = SRG{p}   [{group_label}]   spectrum {spectrum(Q)}")
    elif np.array_equal(Q, 1 - np.eye(m, dtype=int)):
        print(f"  fold = K{m}  (complete graph -- trivial quotient)   spectrum {spectrum(Q)}")
    else:
        print(f"  fold = regular, not SRG; spectrum {spectrum(Q)}")

    ev, od = antipode_split(A.astype(float), sigma)
    print(f"  sigma-EVEN spectrum (= the fold's): {ev}")
    print(f"  sigma-ODD  spectrum (killed by fold): {od}")
    odd_mult = min(od.values())
    print(f"  -> reflection representation = sigma-odd part, "
          f"smallest odd multiplicity {odd_mult} = rank")


if __name__ == "__main__":
    analyse("E8  (4_21 / Gosset edge graph)", np.load("/tmp/e8_R.npy"), 1,
            "E8 root graph, Aut = O+(8,2):2")
    analyse("E7  (Gosset graph, 56 minuscule weights)", np.load("/tmp/e7_W.npy"), 1,
            "")
    analyse("E6  (Schlafli graph, 27 lines)", np.load("/tmp/e6_W.npy"), 0,
            "PSp(4,3)")
