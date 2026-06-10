#!/usr/bin/env python3
"""
Attempt 11: Period matrix of the Bolza curve y^2 = x^5 - x.

Computes the period matrix Omega by numerically integrating the
holomorphic differentials dx/(2y) and x*dx/(2y) along a canonical
homology basis. No Fuchsian group or Poincare disk needed.
"""

import numpy as np
from scipy.integrate import quad


# ══════════════════════════════════════════════════════════════════════
#  The curve y^2 = x^5 - x = x(x-1)(x+1)(x-i)(x+i)
# ══════════════════════════════════════════════════════════════════════

BRANCH_POINTS = [0, 1, -1, 1j, -1j]  # 5 finite branch points + infinity


def f_curve(x):
    """f(x) = x^5 - x. The curve is y^2 = f(x)."""
    return complex(x)**5 - complex(x)


def y_branch(x, sheet=0):
    """Compute y = sqrt(x^5 - x) on a given sheet.

    sheet=0: principal branch (arg in (-pi, pi])
    sheet=1: other sheet (negated)
    """
    val = f_curve(x)
    y = np.sqrt(val)  # numpy's principal sqrt
    if sheet == 1:
        y = -y
    return y


# ══════════════════════════════════════════════════════════════════════
#  Period computation via path integration
# ══════════════════════════════════════════════════════════════════════

def integrate_along_path(path_func, path_deriv, t_range, omega_idx, n_points=2000):
    """Integrate omega along a parameterized path in the x-plane.

    path_func(t) -> x(t)  (complex)
    path_deriv(t) -> dx/dt (complex)
    omega_idx: 0 for dx/(2y), 1 for x*dx/(2y)

    Uses composite trapezoidal rule with continuous branch tracking.
    """
    t0, t1 = t_range
    ts = np.linspace(t0, t1, n_points)
    dt = ts[1] - ts[0]

    xs = np.array([complex(path_func(t)) for t in ts], dtype=complex)
    dxs = np.array([complex(path_deriv(t)) for t in ts], dtype=complex)

    # Compute y with continuous branch tracking
    ys = np.zeros(n_points, dtype=complex)
    val0 = f_curve(xs[0])
    ys[0] = np.sqrt(complex(val0))
    for k in range(1, n_points):
        y_candidate = np.sqrt(complex(f_curve(xs[k])))
        # Choose the branch closer to the previous value (continuity)
        if abs(y_candidate - ys[k-1]) > abs(-y_candidate - ys[k-1]):
            y_candidate = -y_candidate
        ys[k] = y_candidate

    # Integrand: f_k(x) * dx / (2y)
    if omega_idx == 0:
        integrand = dxs / (2 * ys)
    else:
        integrand = xs * dxs / (2 * ys)

    # Trapezoidal rule
    return np.trapz(integrand, ts)


def circle_path(center, radius):
    """Path around a circle: t in [0, 2*pi]."""
    def path(t): return center + radius * np.exp(1j * t)
    def deriv(t): return 1j * radius * np.exp(1j * t)
    return path, deriv, (0, 2 * np.pi)


def line_path(z0, z1):
    """Straight line from z0 to z1: t in [0, 1]."""
    def path(t): return z0 + t * (z1 - z0)
    def deriv(t): return z1 - z0
    return path, deriv, (0, 1)


def arc_path(center, radius, theta0, theta1):
    """Arc from angle theta0 to theta1."""
    def path(t): return center + radius * np.exp(1j * (theta0 + t * (theta1 - theta0)))
    def deriv(t): return 1j * (theta1 - theta0) * radius * np.exp(1j * (theta0 + t * (theta1 - theta0)))
    return path, deriv, (0, 1)


# ══════════════════════════════════════════════════════════════════════
#  Homology cycles for y^2 = x^5 - x
# ══════════════════════════════════════════════════════════════════════

def compute_a_period(b1, b2, omega_idx, radius=0.05, n_points=4000):
    """Compute A-period: integral around a small loop encircling the cut b1--b2.

    The loop goes around the midpoint of [b1, b2] at distance radius,
    starting and ending on the upper sheet.
    """
    mid = (b1 + b2) / 2
    r = abs(b2 - b1) / 2 + radius  # just outside the cut endpoints

    path, deriv, t_range = circle_path(mid, r)
    return integrate_along_path(path, deriv, t_range, omega_idx, n_points)


def compute_b_period(cut_a, cut_b, omega_idx, n_points=4000):
    """Compute B-period: integral along a path connecting two cuts.

    Goes from the midpoint of cut_a to the midpoint of cut_b,
    avoiding other branch points.
    """
    mid_a = (cut_a[0] + cut_a[1]) / 2
    mid_b = (cut_b[0] + cut_b[1]) / 2

    # Route via a point away from all branch points
    waypoint = mid_a + 2 * (mid_b - mid_a) * 1j  # perpendicular detour

    # Path: mid_a -> waypoint -> mid_b, then back on the other sheet
    # Actually for B-cycles, the path goes from one cut to another on the upper sheet
    # and returns on the lower sheet. The net integral is the B-period.
    #
    # Simpler: just integrate along the upper sheet from mid_a to mid_b.
    # The B-period is this integral (times 2 if we need the full cycle).

    path, deriv, t_range = line_path(mid_a, mid_b)
    return integrate_along_path(path, deriv, t_range, omega_idx, n_points)


# ══════════════════════════════════════════════════════════════════════
#  Alternative: direct integration between branch points
# ══════════════════════════════════════════════════════════════════════

def period_between_branch_points(b1, b2, omega_idx, n_points=4000):
    """Integrate omega on the upper sheet along a straight line from b1 to b2.

    This is half the A-period for cuts along this segment.
    Handles the sqrt singularities at the endpoints by approaching from offset.
    """
    eps = 1e-8
    # Slightly inside to avoid singularities at endpoints
    z0 = b1 + eps * (b2 - b1)
    z1 = b2 - eps * (b2 - b1)
    path, deriv, t_range = line_path(z0, z1)
    return integrate_along_path(path, deriv, t_range, omega_idx, n_points)


# ══════════════════════════════════════════════════════════════════════
#  Compute the full period matrix
# ══════════════════════════════════════════════════════════════════════

def compute_period_matrix():
    """Compute the period matrix of y^2 = x^5 - x.

    Uses Tretkoff-Tretkoff style direct integration:
    - A-period = integral across the cut on the upper sheet (times a factor)
    - B-period = 2 * integral from one cut to another on the upper sheet

    Branch point pairing:
    Cut 1: 0 to 1 (A1 cycle encircles this)
    Cut 2: -1 to -i (A2 cycle encircles this)
    Cut 3: i to infinity (not directly used for A-cycles)

    Returns Omega (2x2 complex), A_matrix, B_matrix.
    """
    print("Computing period matrix of y^2 = x^5 - x...")
    print()

    # A-periods: integral across each cut on the upper sheet.
    # The full A-cycle integral = integral from b1 to b2 of dx/y (upper sheet).
    # This equals 2 * integral of omega = 2 * integral of dx/(2y).
    #
    # Parameterize: for cut from b1 to b2, use x = b1 + t*(b2-b1), t in [eps, 1-eps]
    print("  A-periods (integrals across cuts):")
    A = np.zeros((2, 2), dtype=complex)

    cuts = [(0, 1), (-1, -1j)]
    for k, (b1, b2) in enumerate(cuts):
        for j in range(2):
            # Integrate across the cut on the upper sheet
            # The A-period = 2 * integral_{b1}^{b2} omega_j (upper sheet)
            val = period_between_branch_points(b1, b2, j, n_points=8000)
            A[j, k] = 2 * val
        print(f"    A[0,{k}] = {A[0,k]:.8f}")
        print(f"    A[1,{k}] = {A[1,k]:.8f}")

    print()

    # B-periods: path from cut_k to cut_3 on the upper sheet, times 2.
    # The full B-cycle goes upper-sheet forward + lower-sheet backward = 2x upper.
    #
    # B1: from midpoint of cut1 (0.5) to near cut3 (i)
    # B2: from midpoint of cut2 (-0.5 - 0.5i) to near cut3 (i)
    #
    # Path must avoid crossing other cuts!
    print("  B-periods (paths between cuts, x2):")
    B = np.zeros((2, 2), dtype=complex)

    # B1: 0.5 -> i, avoiding cut2 (from -1 to -i, in Q3)
    # Straight line from 0.5 to i stays in Q1/Q2, doesn't cross cut2 ✓
    for j in range(2):
        val = period_between_branch_points(0.5, 1j, j, n_points=8000)
        B[j, 0] = 2 * val
    print(f"    B[0,0] = {B[0,0]:.8f}")
    print(f"    B[1,0] = {B[1,0]:.8f}")

    # B2: midpoint of cut2 -> i
    # Cut2 midpoint: (-1 + (-i))/2 = -0.5 - 0.5i
    # Straight line to i: goes through the origin area, must not cross cut1 (0 to 1)
    # Line from -0.5-0.5i to i: parameterize and check...
    # At t=0: -0.5 - 0.5i. At t=1: i. The real part goes from -0.5 to 0.
    # The line doesn't cross the positive real axis [0,1], so it doesn't cross cut1 ✓
    mid2 = (-1 + (-1j)) / 2
    for j in range(2):
        val = period_between_branch_points(mid2, 1j, j, n_points=8000)
        B[j, 1] = 2 * val
    print(f"    B[0,1] = {B[0,1]:.8f}")
    print(f"    B[1,1] = {B[1,1]:.8f}")

    print()
    print(f"  A matrix:\n{A}")
    print(f"  B matrix:\n{B}")

    # Period matrix: Omega = A^{-1} @ B
    try:
        Omega = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("  ERROR: A matrix is singular!")
        return None, A, B

    print(f"\n  Period matrix Omega:")
    print(f"    [{Omega[0,0]:.8f}  {Omega[0,1]:.8f}]")
    print(f"    [{Omega[1,0]:.8f}  {Omega[1,1]:.8f}]")

    # Validation
    print("\n  Validation:")
    sym_err = np.abs(Omega - Omega.T).max()
    print(f"    Symmetric? |Omega - Omega^T| = {sym_err:.6e}")
    Im_Omega = Omega.imag
    eigvals = np.linalg.eigvalsh(Im_Omega)
    print(f"    Im(Omega) eigenvalues: {eigvals}")
    print(f"    Im(Omega) > 0? {np.all(eigvals > 0)}")

    # For Bolza surface, check D4 symmetry
    print(f"    Omega[0,0] = {Omega[0,0]:.8f}")
    print(f"    Omega[1,1] = {Omega[1,1]:.8f}")
    print(f"    Omega[0,1] = {Omega[0,1]:.8f}")
    print(f"    Omega[1,0] = {Omega[1,0]:.8f}")

    return Omega, A, B


# ══════════════════════════════════════════════════════════════════════
#  Riemann theta function
# ══════════════════════════════════════════════════════════════════════

def riemann_theta(z, Omega, N_max=6):
    """Compute the Riemann theta function Theta(z; Omega).

    Theta(z; Omega) = sum_{n in Z^g} exp(pi*i * n^T Omega n + 2*pi*i * n^T z)

    Args:
        z: complex vector of length g (= 2)
        Omega: g x g period matrix
        N_max: truncation radius (sum over |n| <= N_max)

    Returns: complex number
    """
    z = np.asarray(z, dtype=complex)
    g = len(z)
    total = 0j

    # Iterate over all n in Z^g with |n| <= N_max
    for n1 in range(-N_max, N_max + 1):
        for n2 in range(-N_max, N_max + 1):
            n = np.array([n1, n2], dtype=complex)
            quad_form = np.dot(n, Omega @ n)
            lin_form = np.dot(n, z)
            total += np.exp(np.pi * 1j * quad_form + 2 * np.pi * 1j * lin_form)

    return total


# ══════════════════════════════════════════════════════════════════════
#  Abel-Jacobi map on the curve
# ══════════════════════════════════════════════════════════════════════

def abel_jacobi(x_target, y_target, base_x=0.5, n_points=2000):
    """Compute the Abel-Jacobi map from base_point to (x_target, y_target).

    u = integral_{base}^{target} (dx/(2y), x*dx/(2y))

    Integrates along a straight line in the x-plane from base_x to x_target.
    """
    path, deriv, t_range = line_path(base_x, x_target)
    u0 = integrate_along_path(path, deriv, t_range, 0, n_points)
    u1 = integrate_along_path(path, deriv, t_range, 1, n_points)
    return np.array([u0, u1])


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  Bolza Surface — Attempt 11: Period Matrix")
    print("=" * 60)
    print()

    Omega, A, B = compute_period_matrix()
    if Omega is None:
        return

    print("\n" + "=" * 60)
    print("  Theta function test")
    print("=" * 60)

    theta_0 = riemann_theta([0, 0], Omega)
    print(f"  Theta(0; Omega) = {theta_0:.8f}")
    print(f"  |Theta(0)| = {abs(theta_0):.8f}")

    # Test at a few points
    for z in [[0.1, 0.2], [0.5, 0.5], [0.3, 0.7j]]:
        z = np.array(z, dtype=complex)
        th = riemann_theta(z, Omega)
        print(f"  Theta({z}) = {th:.6f}, |Theta| = {abs(th):.6f}")

    print("\n" + "=" * 60)
    print("  Abel-Jacobi at branch points")
    print("=" * 60)

    # Compute Abel-Jacobi images of the Weierstrass points
    # (relative to base point x = 0.5)
    for bp in [0, 1, -1, 1j, -1j]:
        y_bp = 0  # Weierstrass points have y = 0
        u = abel_jacobi(bp, y_bp, base_x=0.5 + 0.01j)  # slight offset to avoid real axis
        print(f"  AJ({bp}) = [{u[0]:.6f}, {u[1]:.6f}]")

    print("\nDone!")


if __name__ == '__main__':
    main()
