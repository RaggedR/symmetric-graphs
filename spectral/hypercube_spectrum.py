#!/usr/bin/env python3
"""
Spectral eigenvalues of hypercubes Q_n for various dimensions.

The n-dimensional hypercube Q_n has 2^n vertices, each of degree n.
Its adjacency spectrum has a beautiful closed form:
    eigenvalue (n - 2k) with multiplicity C(n, k)  for k = 0, 1, ..., n

This follows from Q_n = K_2 □ K_2 □ ... □ K_2 (n-fold Cartesian product)
and the Cartesian product spectral theorem.

We verify the closed form against direct numerical computation of the
adjacency matrix for small dimensions.
"""

import numpy as np
from math import comb
from itertools import product as cartesian_product


def hypercube_adjacency(n):
    """Build the 2^n × 2^n adjacency matrix of Q_n.

    Two binary strings are adjacent iff they differ in exactly one bit.
    """
    N = 1 << n
    A = np.zeros((N, N), dtype=float)
    for v in range(N):
        for bit in range(n):
            w = v ^ (1 << bit)
            A[v, w] = 1.0
    return A


def closed_form_spectrum(n):
    """Return (eigenvalue, multiplicity) pairs from the closed form."""
    return [(n - 2 * k, comb(n, k)) for k in range(n + 1)]


def numerical_spectrum(n):
    """Compute adjacency eigenvalues numerically and round to integers."""
    A = hypercube_adjacency(n)
    eigs = np.linalg.eigvalsh(A)
    # Round to nearest integer (exact for hypercubes)
    return np.sort(np.round(eigs).astype(int))[::-1]


def spectral_summary(n):
    """Print a full spectral summary for Q_n."""
    spec = closed_form_spectrum(n)
    N = 1 << n

    print(f"\n{'='*60}")
    print(f"  Q_{n}  —  {N} vertices, {n}-regular")
    print(f"{'='*60}")

    print(f"\n  Adjacency spectrum (closed form):")
    print(f"  {'Eigenvalue':>12}  {'Multiplicity':>12}")
    print(f"  {'─'*12}  {'─'*12}")
    for val, mult in spec:
        print(f"  {val:>12}  {mult:>12}")

    # Key spectral quantities
    lambda_max = spec[0][0]       # = n
    lambda_2 = spec[1][0]         # = n - 2
    lambda_min = spec[-1][0]      # = -n
    spectral_gap = lambda_max - lambda_2  # always 2 for hypercubes

    print(f"\n  Key quantities:")
    print(f"    λ_max          = {lambda_max}")
    print(f"    λ_2            = {lambda_2}")
    print(f"    λ_min          = {lambda_min}")
    print(f"    Spectral gap   = {spectral_gap}")

    # Laplacian spectrum: for n-regular, Laplacian eigenvalues = n - λ_adj
    print(f"\n  Laplacian spectrum:")
    print(f"  {'Eigenvalue':>12}  {'Multiplicity':>12}")
    print(f"  {'─'*12}  {'─'*12}")
    for val, mult in spec:
        lap_val = n - val
        print(f"  {lap_val:>12}  {mult:>12}")

    fiedler = n - lambda_2  # = 2 always
    print(f"\n    Algebraic connectivity (Fiedler) = {fiedler}")

    # Ramanujan check: λ_2 ≤ 2√(n-1)
    ramanujan_bound = 2 * np.sqrt(n - 1) if n > 1 else 0
    is_ramanujan = lambda_2 <= ramanujan_bound + 1e-10
    print(f"    Ramanujan bound 2√(n−1) = {ramanujan_bound:.4f}")
    print(f"    Is Ramanujan?  {'YES' if is_ramanujan else 'NO'} (λ₂ = {lambda_2})")

    # Verify against numerical computation for small n
    if n <= 12:
        num_eigs = numerical_spectrum(n)
        closed_eigs = []
        for val, mult in spec:
            closed_eigs.extend([val] * mult)
        closed_eigs.sort(reverse=True)
        match = np.allclose(num_eigs, closed_eigs)
        print(f"\n    Numerical verification: {'PASS ✓' if match else 'FAIL ✗'}")


if __name__ == "__main__":
    for n in range(1, 11):
        spectral_summary(n)

    # Summary table
    print(f"\n\n{'='*60}")
    print(f"  Summary Table")
    print(f"{'='*60}")
    print(f"  {'Q_n':>4}  {'|V|':>6}  {'deg':>4}  {'λ₂':>4}  {'gap':>4}  {'Fiedler':>8}  {'Ramanujan':>10}")
    print(f"  {'─'*4}  {'─'*6}  {'─'*4}  {'─'*4}  {'─'*4}  {'─'*8}  {'─'*10}")
    for n in range(1, 16):
        N = 1 << n
        spec = closed_form_spectrum(n)
        lam2 = spec[1][0]
        gap = 2
        fiedler = 2
        ram = 2 * np.sqrt(n - 1) if n > 1 else 0
        is_ram = "YES" if lam2 <= ram + 1e-10 else "NO"
        print(f"  Q_{n:<2}  {N:>6}  {n:>4}  {lam2:>4}  {gap:>4}  {fiedler:>8}  {is_ram:>10}")
