"""
Compute the 27 lines on the Clebsch diagonal surface.

The Clebsch surface in P^4: sum(x_i) = 0, sum(x_i^3) = 0.
Affine chart: x3=1, x4=-(x0+x1+x2+1), so f(x,y,z) = x^3+y^3+z^3+1-(x+y+z+1)^3 = 0.

The 27 lines split under S5 into:
  - 15 lines (Type I): for each partition {i,j} | {m} | {k,l} of {0..4},
    the line span(e_i - e_j, e_k - e_l) in the hyperplane sum=0.
  - 12 lines (Type II): for each unordered 5-cycle {sigma, sigma^{-1}} in S5,
    the real line from the DFT eigenvectors of the permutation matrix.
"""

import numpy as np
from itertools import combinations, permutations

def proj_to_affine(p):
    """Project from P^4 (with sum=0) to affine R^3 via x3=1."""
    # Point in R^5 with sum=0. Dehomogenize by x3.
    if abs(p[3]) < 1e-12:
        return None  # point at infinity in this chart
    s = p[3]
    return np.array([p[0]/s, p[1]/s, p[2]/s])

def line_to_affine_segment(A, B, clip=2.0):
    """Convert a projective line span(A,B) to an affine segment clipped to a sphere.

    A, B are points in R^5 with sum=0. The line is {sA + tB : [s:t] in P^1}.
    Affine chart: divide by x3 coordinate.
    """
    # The affine point for parameter t is: P(t) = A + t*B, then (P[0]/P[3], P[1]/P[3], P[2]/P[3])
    # This is a rational map R -> R^3 with a pole at t = -A[3]/B[3] (if B[3]!=0).
    #
    # Strategy: sample densely, including around the pole, and collect clipped points.

    # Handle the case where the line is entirely at infinity (x3=0 for all points)
    if abs(A[3]) < 1e-12 and abs(B[3]) < 1e-12:
        return None

    # If B[3] ~ 0, reparametrize: swap A and B or use a different basis
    if abs(B[3]) < abs(A[3]):
        A, B = B, A

    # Now B[3] != 0. The pole is at t_pole = -A[3]/B[3].
    # Sample away from the pole on both sides.
    t_pole = -A[3] / B[3] if abs(B[3]) > 1e-12 else 1e12

    # Sample on both sides of the pole
    segments = []
    for t_range in [np.linspace(t_pole + 0.01, t_pole + 100, 2000),
                    np.linspace(t_pole - 100, t_pole - 0.01, 2000)]:
        points = []
        for t in t_range:
            P = A + t * B
            if abs(P[3]) < 1e-10:
                continue
            aff = np.array([P[0]/P[3], P[1]/P[3], P[2]/P[3]])
            if np.linalg.norm(aff) < clip:
                points.append(aff)
        if len(points) >= 2:
            segments.append((np.array(points[0]), np.array(points[-1])))

    if not segments:
        return None

    # Return the longest segment
    best = max(segments, key=lambda s: np.linalg.norm(s[1] - s[0]))
    return best


def compute_15_lines():
    """The 15 Type I lines: for each partition {i,j}|{m}|{k,l} of {0..4}.

    A line is determined by span(e_i-e_j, e_k-e_l), which equals span(e_k-e_l, e_i-e_j).
    So the partition {i,j}|{m}|{k,l} is the same line as {k,l}|{m}|{i,j}.
    We just need: singleton m, and one (unordered) pair from the remaining 4.
    That gives 5 * 3 = 15 lines.
    """
    lines = []
    seen = set()

    for m in range(5):
        rest = sorted(i for i in range(5) if i != m)
        # All ways to split rest into two unordered pairs (3 per singleton)
        for i, j in combinations(rest, 2):
            k, l = sorted(x for x in rest if x != i and x != j)
            # Only take one ordering: require min(i,j) < min(k,l)
            if min(i,j) > min(k,l):
                continue

            A = np.zeros(5); A[i] = 1; A[j] = -1
            B = np.zeros(5); B[k] = 1; B[l] = -1
            lines.append((A, B, f"L({i}{j}|{m}|{k}{l})"))

    return lines


def compute_12_lines():
    """The 12 Type II lines from unordered 5-cycles.

    For a 5-cycle sigma, the line is span(u, w) where:
      u_k = cos(2*pi*k/5) applied to sigma-permuted coordinates
      w_k = sin(2*pi*k/5) applied to sigma-permuted coordinates
    """
    lines = []

    # Generate all 5-cycles in S5 as tuples (a,b,c,d,e) meaning a->b->c->d->e->a
    # Two 5-cycles are "the same unordered" if one is the inverse of the other
    # (or a cyclic rotation, but we fix the starting element)

    five_cycles = []
    for perm in permutations(range(5)):
        # Check if this is a 5-cycle
        # A permutation is a 5-cycle if its cycle decomposition is a single 5-cycle
        visited = [False]*5
        cycle_len = 0
        pos = 0
        while not visited[pos]:
            visited[pos] = True
            pos = perm[pos]
            cycle_len += 1
        if cycle_len == 5 and all(visited):
            # Normalize: start with 0
            cycle = [0]
            pos = perm[0]
            while pos != 0:
                cycle.append(pos)
                pos = perm[pos]
            five_cycles.append(tuple(cycle))

    # Remove duplicates (cyclic rotations all start with 0 now, so just unique tuples)
    five_cycles = list(set(five_cycles))  # Should be 24/5 = ... wait
    # Starting with 0: there are 4! = 24 permutations of remaining 4, but only
    # 24/5*5 = 24 five-cycles total, and fixing start=0 gives 24/5... hmm
    # Actually: 5-cycles starting with 0 are just (0, a, b, c, d) for permutations of {1,2,3,4}
    # that form a valid cycle. There are 4! = 24 of them, but we mod out by cyclic rotation.
    # Since we fixed starting point = 0, there's no more cyclic freedom. So we have 24 cycles.
    # But wait: |5-cycles in S5| = 5!/5 = 24. All start with 0 after normalization. ✓

    # Now group into unordered pairs {sigma, sigma^{-1}}
    # sigma = (0,a,b,c,d) => sigma^{-1} = (0,d,c,b,a)
    seen = set()
    unordered_cycles = []
    for cycle in five_cycles:
        inv = (cycle[0],) + tuple(reversed(cycle[1:]))
        key = tuple(sorted([cycle, inv]))
        if key not in seen:
            seen.add(key)
            unordered_cycles.append(cycle)

    print(f"Found {len(unordered_cycles)} unordered 5-cycles")

    for cycle in unordered_cycles:
        # DFT eigenvectors
        u = np.zeros(5)
        w = np.zeros(5)
        for k in range(5):
            angle = 2 * np.pi * k / 5
            coord = cycle[k]
            u[coord] = np.cos(angle)
            w[coord] = np.sin(angle)

        # Verify sum = 0
        assert abs(np.sum(u)) < 1e-10, f"u sum = {np.sum(u)}"
        assert abs(np.sum(w)) < 1e-10, f"w sum = {np.sum(w)}"

        lines.append((u, w, f"L5({','.join(map(str,cycle))})"))

    return lines


def verify_line_on_surface(A, B, name=""):
    """Verify that span(A,B) lies on sum=0, sum_cubes=0."""
    # Check sum=0
    assert abs(np.sum(A)) < 1e-10, f"{name}: sum(A)={np.sum(A)}"
    assert abs(np.sum(B)) < 1e-10, f"{name}: sum(B)={np.sum(B)}"

    # Check: for all s,t, sum((sA+tB)_i^3) = 0
    # Coefficients of s^3, s^2*t, s*t^2, t^3:
    c3 = np.sum(A**3)
    c2 = np.sum(3 * A**2 * B)  # 3*sum(a_i^2 * b_i)
    c1 = np.sum(3 * A * B**2)  # 3*sum(a_i * b_i^2)
    c0 = np.sum(B**3)

    ok = all(abs(c) < 1e-8 for c in [c3, c2, c1, c0])
    if not ok:
        print(f"  FAILED {name}: coeffs = {c3:.6f}, {c2:.6f}, {c1:.6f}, {c0:.6f}")
    return ok


def main():
    print("Computing 27 lines on the Clebsch diagonal surface...")
    print("=" * 60)

    lines_15 = compute_15_lines()
    print(f"\nType I (15 lines from pair partitions): {len(lines_15)} lines")

    ok_count = 0
    for A, B, name in lines_15:
        if verify_line_on_surface(A, B, name):
            ok_count += 1
    print(f"  Verified: {ok_count}/15")

    lines_12 = compute_12_lines()
    print(f"\nType II (12 lines from 5-cycles): {len(lines_12)} lines")

    ok_count = 0
    for A, B, name in lines_12:
        if verify_line_on_surface(A, B, name):
            ok_count += 1
    print(f"  Verified: {ok_count}/12")

    all_lines = lines_15 + lines_12
    print(f"\nTotal: {len(all_lines)} lines")

    # Project to affine coordinates and output as JSON for the visualization
    print("\n" + "=" * 60)
    print("Affine segments (x3=1 chart, clipped to [-2,2]^3):")
    print("=" * 60)

    clip = 8.0
    segments = []
    for A, B, name in all_lines:
        seg = line_to_affine_segment(A, B, clip=clip)
        if seg is not None:
            p1, p2 = seg
            segments.append((p1.tolist(), p2.tolist(), name))
            print(f"  {name}: ({p1[0]:.4f},{p1[1]:.4f},{p1[2]:.4f}) -> ({p2[0]:.4f},{p2[1]:.4f},{p2[2]:.4f})")
        else:
            print(f"  {name}: not visible in this chart")

    print(f"\n{len(segments)} lines visible in affine chart")

    # Output as JavaScript array for the HTML visualization
    print("\n// JavaScript data for clebsch.html:")
    print("const LINES_27 = [")
    for p1, p2, name in segments:
        print(f"  [{p1[0]:.6f},{p1[1]:.6f},{p1[2]:.6f}, {p2[0]:.6f},{p2[1]:.6f},{p2[2]:.6f}], // {name}")
    print("];")


if __name__ == "__main__":
    main()
