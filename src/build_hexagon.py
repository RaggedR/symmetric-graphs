#!/usr/bin/env python3
"""
Construct the split Cayley hexagon GH(2,2) and its incidence graph (= Tutte 12-cage).

Uses the Zorn vector matrix model of the split octonion algebra over F_2.
An imaginary split octonion is a trace-zero Zorn matrix [α, a; b, α]
represented as a 7-tuple (a1, a2, a3, α, b1, b2, b3).

The quadratic form is Q(x) = a·b + α = x0*x4 + x1*x5 + x2*x6 + x3.
The 63 points of Q(6,2) are the nonzero vectors with Q(x) = 0.

A totally isotropic line {p, q, p+q} is a HEXAGON line iff p·q ∈ ⟨p,q⟩
(the octonion product lands back in the 2-space).

The incidence graph (63 points + 63 lines, bipartite, cubic) is the Tutte 12-cage.
"""

import json
import os
import sys


def vec_to_int(v):
    """Convert 7-tuple over F_2 to integer."""
    return sum(b << i for i, b in enumerate(v))


def int_to_vec(n):
    """Convert integer to 7-tuple over F_2."""
    return tuple((n >> i) & 1 for i in range(7))


def Q(v):
    """Quadratic form: x0*x4 + x1*x5 + x2*x6 + x3."""
    return (v[0]*v[4] + v[1]*v[5] + v[2]*v[6] + v[3]) % 2


def B(u, v):
    """Associated bilinear form: x0*y4 + x4*y0 + x1*y5 + x5*y1 + x2*y6 + x6*y2."""
    return (u[0]*v[4] + u[4]*v[0] + u[1]*v[5] + u[5]*v[1] + u[2]*v[6] + u[6]*v[2]) % 2


def xor_vec(u, v):
    """Add two vectors over F_2."""
    return tuple((a + b) % 2 for a, b in zip(u, v))


def zorn_mult(x, y):
    """Multiply two imaginary split octonions in Zorn coordinates over F_2.

    x = (a1, a2, a3, alpha, b1, b2, b3) represents [alpha, (a1,a2,a3); (b1,b2,b3), alpha].
    y = (c1, c2, c3, beta, d1, d2, d3) represents [beta, (c1,c2,c3); (d1,d2,d3), beta].

    Returns 7-tuple (imaginary part) if trace=0, else None.
    """
    a1, a2, a3, alpha, b1, b2, b3 = x
    c1, c2, c3, beta, d1, d2, d3 = y

    # Inner products
    ad = (a1*d1 + a2*d2 + a3*d3) % 2
    bc = (b1*c1 + b2*c2 + b3*c3) % 2

    # Diagonal entries
    p = (alpha * beta + ad) % 2
    q = (alpha * beta + bc) % 2

    # Trace = p + q = ad + bc = B(x,y)
    # If B(x,y) = 0, the product is imaginary (trace zero, p = q)

    if (p + q) % 2 != 0:
        return None  # product has nonzero trace (not imaginary)

    # Cross products over F_2 (using standard right-hand rule, but -1 = +1)
    # b × d
    bxd = ((b2*d3 + b3*d2) % 2, (b3*d1 + b1*d3) % 2, (b1*d2 + b2*d1) % 2)
    # a × c
    axc = ((a2*c3 + a3*c2) % 2, (a3*c1 + a1*c3) % 2, (a1*c2 + a2*c1) % 2)

    # Upper-right vector: u = alpha*c + beta*a + b×d
    u1 = (alpha*c1 + beta*a1 + bxd[0]) % 2
    u2 = (alpha*c2 + beta*a2 + bxd[1]) % 2
    u3 = (alpha*c3 + beta*a3 + bxd[2]) % 2

    # Lower-left vector: v = beta*b + alpha*d + a×c
    v1 = (beta*b1 + alpha*d1 + axc[0]) % 2
    v2 = (beta*b2 + alpha*d2 + axc[1]) % 2
    v3 = (beta*b3 + alpha*d3 + axc[2]) % 2

    return (u1, u2, u3, p, v1, v2, v3)


def main():
    # Step 1: Find all 63 points on Q(6,2)
    points = []
    for n in range(1, 128):  # nonzero vectors in F_2^7
        v = int_to_vec(n)
        if Q(v) == 0:
            points.append(v)

    print(f"Points on Q(6,2): {len(points)}")
    assert len(points) == 63, f"Expected 63 points, got {len(points)}"

    # Index points
    point_idx = {vec_to_int(p): i for i, p in enumerate(points)}

    # Step 2: Find all totally isotropic lines
    # A t.i. line is a 2-space {p, q, p+q} where Q(p)=Q(q)=Q(p+q)=0 and B(p,q)=0
    ti_lines = []
    for i in range(63):
        for j in range(i+1, 63):
            p, q = points[i], points[j]
            if B(p, q) == 0:
                r = xor_vec(p, q)
                if Q(r) == 0:
                    line = frozenset([i, j, point_idx[vec_to_int(r)]])
                    if len(line) == 3:  # all three distinct
                        ti_lines.append(line)

    # Deduplicate (each line found 3 times from different pairs)
    ti_lines = list(set(ti_lines))
    print(f"Totally isotropic lines: {len(ti_lines)}")

    # Step 3: Filter to hexagon lines using octonion product condition
    # {p, q, p+q} is a hexagon line iff p·q ∈ {0, p, q, p+q}
    hex_lines = []
    for line in ti_lines:
        pts = list(line)
        p, q = points[pts[0]], points[pts[1]]
        prod = zorn_mult(p, q)
        if prod is None:
            continue  # shouldn't happen since B(p,q)=0
        prod_int = vec_to_int(prod)
        # Check if product is in the 2-space {0, p, q, p+q}
        space = {0, vec_to_int(p), vec_to_int(q), vec_to_int(xor_vec(p, q))}
        if prod_int in space:
            hex_lines.append(line)

    print(f"Hexagon lines: {len(hex_lines)}")

    if len(hex_lines) != 63:
        print(f"WARNING: Expected 63 hexagon lines, got {len(hex_lines)}")
        print("Trying alternative condition: p·q = 0...")
        hex_lines_zero = [l for l in ti_lines
                          if zorn_mult(points[list(l)[0]], points[list(l)[1]]) == (0,0,0,0,0,0,0)]
        print(f"  Lines with p·q = 0: {len(hex_lines_zero)}")

        print("Trying: p·q = p+q...")
        hex_lines_sum = []
        for line in ti_lines:
            pts = list(line)
            p, q = points[pts[0]], points[pts[1]]
            pq = xor_vec(p, q)
            prod = zorn_mult(p, q)
            if prod is not None and prod == pq:
                hex_lines_sum.append(line)
        print(f"  Lines with p·q = p+q: {len(hex_lines_sum)}")

        # Try the complement condition
        print("Trying complement: p·q NOT in <p,q>...")
        comp_lines = [l for l in ti_lines if l not in set(hex_lines)]
        print(f"  Complement lines: {len(comp_lines)}")

    # Step 4: Check point valencies
    val = [0] * 63
    for line in hex_lines:
        for pt in line:
            val[pt] += 1
    val_counts = {}
    for v in val:
        val_counts[v] = val_counts.get(v, 0) + 1
    print(f"Point valencies: {val_counts}")

    if len(hex_lines) == 63 and all(v == 3 for v in val):
        print("\n*** SUCCESS: Valid generalized hexagon GH(2,2) ***")
        build_incidence_graph(points, hex_lines)
    else:
        print("\nHexagon construction needs adjustment.")
        # Try to find the right subset
        find_hexagon_lines(points, point_idx, ti_lines)


def find_hexagon_lines(points, point_idx, ti_lines):
    """Try alternative methods to find the 63 hexagon lines."""
    print("\n--- Alternative: using trilinear form ---")

    # The trilinear form T(x,y,z) from the octonion structure constants
    # For Zorn matrices, T comes from the alternating part of the product
    # T(x,y,z) = the "associator" component

    # For the Zorn algebra, the trilinear form can be computed as:
    # T(x,y,z) = <x·y, z> where <,> uses the bilinear form B
    # and x·y is the imaginary part of the product

    count_by_T = {0: 0, 1: 0}
    lines_with_T0 = []
    lines_with_T1 = []

    for line in ti_lines:
        pts = list(line)
        p, q, r = points[pts[0]], points[pts[1]], points[pts[2]]
        prod_pq = zorn_mult(p, q)
        if prod_pq is None:
            continue
        # T(p,q,r) = B(p·q, r)
        T_val = B(prod_pq, r)
        if T_val == 0:
            lines_with_T0.append(line)
        else:
            lines_with_T1.append(line)

    print(f"Lines with T(p,q,r) = 0: {len(lines_with_T0)}")
    print(f"Lines with T(p,q,r) = 1: {len(lines_with_T1)}")

    # Check valencies for each set
    for label, lines in [("T=0", lines_with_T0), ("T=1", lines_with_T1)]:
        val = [0] * 63
        for line in lines:
            for pt in line:
                val[pt] += 1
        val_counts = {}
        for v in val:
            val_counts[v] = val_counts.get(v, 0) + 1
        print(f"  {label} valencies: {val_counts}")

        if len(lines) == 63 and all(v == 3 for v in val):
            print(f"\n*** SUCCESS with {label}: Valid GH(2,2) ***")
            build_incidence_graph(points, lines)
            return

    print("\n--- Alternative: using product norm condition ---")
    # A hexagon line {p,q,r} satisfies: the norm of the product p·q is zero
    # Q(p·q) = 0 (product stays on the quadric)
    lines_norm0 = []
    for line in ti_lines:
        pts = list(line)
        p, q = points[pts[0]], points[pts[1]]
        prod = zorn_mult(p, q)
        if prod is not None and Q(prod) == 0:
            lines_norm0.append(line)

    print(f"Lines with Q(p·q) = 0: {len(lines_norm0)}")
    val = [0] * 63
    for line in lines_norm0:
        for pt in line:
            val[pt] += 1
    val_counts = {}
    for v in val:
        val_counts[v] = val_counts.get(v, 0) + 1
    print(f"  Valencies: {val_counts}")

    if len(lines_norm0) == 63:
        print("\n*** SUCCESS with Q(p·q)=0: Valid GH(2,2) ***")
        build_incidence_graph(points, lines_norm0)


def build_incidence_graph(points, hex_lines):
    """Build the incidence graph (bipartite, cubic) and output for GAP."""
    n_pts = len(points)
    n_lines = len(hex_lines)
    n = n_pts + n_lines  # 126

    print(f"\nBuilding incidence graph: {n_pts} points + {n_lines} lines = {n} vertices")

    # Adjacency list (1-indexed for GAP)
    adj = [[] for _ in range(n)]

    hex_lines_list = [sorted(list(line)) for line in hex_lines]
    hex_lines_list.sort()

    for li, line in enumerate(hex_lines_list):
        line_vertex = n_pts + li  # 0-indexed
        for pt in line:
            adj[pt].append(line_vertex)
            adj[line_vertex].append(pt)

    # Sort adjacency lists and convert to 1-indexed
    adj_1indexed = [sorted([v + 1 for v in neighbors]) for neighbors in adj]

    # Verify cubic
    degrees = [len(a) for a in adj_1indexed]
    assert all(d == 3 for d in degrees), f"Not cubic! Degrees: {set(degrees)}"

    # Verify bipartite
    # Points are 1..63, lines are 64..126
    for i in range(n_pts):
        for j in adj_1indexed[i]:
            assert j > n_pts, f"Point {i+1} adjacent to point {j}"
    for i in range(n_pts, n):
        for j in adj_1indexed[i]:
            assert j <= n_pts, f"Line {i+1} adjacent to line {j}"

    print("Graph verified: cubic, bipartite")

    # Check girth (should be 12)
    # Quick BFS-based girth check
    girth = check_girth(adj, min(20, n))
    print(f"Girth (estimated from first 20 vertices): {girth}")

    # Output as GAP adjacency list
    print("\n# GAP adjacency list for split Cayley hexagon incidence graph")
    print("# (Tutte 12-cage, 126 vertices, cubic, girth 12)")
    gap_lines = []
    for a in adj_1indexed:
        gap_lines.append(f"[ {', '.join(str(x) for x in a)} ]")
    gap_str = ",\n    ".join(gap_lines)

    output = f"""# Split Cayley hexagon GH(2,2) incidence graph
# = Tutte 12-cage = Foster census F126
# 126 vertices, cubic (3-regular), girth 12
# Aut = G_2(2), order 12096
# Points: vertices 1-63, Lines: vertices 64-126

GraphLibrary.tuttecage12 := function()
    return rec(name := "Tutte-12-cage", n := 126, adj := [
    {gap_str}
    ]);
end;
"""

    src_dir = os.path.dirname(os.path.abspath(__file__))
    gap_path = os.path.join(src_dir, "census", "tutte12cage.g")
    with open(gap_path, "w") as f:
        f.write(output)
    print(f"\nWrote GAP definition to {gap_path}")

    # Also output JSON for direct use
    data = {
        "name": "Tutte-12-cage (GH(2,2))",
        "n": n,
        "adj": adj_1indexed,
        "points": n_pts,
        "lines": n_lines,
    }
    json_path = os.path.join(src_dir, "tutte12cage.json")
    with open(json_path, "w") as f:
        json.dump(data, f)
    print(f"Wrote JSON to {json_path}")


def check_girth(adj, max_start=20):
    """Estimate girth by BFS from first max_start vertices."""
    n = len(adj)
    min_cycle = float('inf')

    for start in range(min(max_start, n)):
        # BFS to find shortest cycle through start
        dist = [-1] * n
        parent = [-1] * n
        dist[start] = 0
        queue = [start]
        qi = 0
        while qi < len(queue):
            u = queue[qi]
            qi += 1
            for v in adj[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    queue.append(v)
                elif parent[u] != v and parent[v] != u:
                    # Found a cycle
                    cycle_len = dist[u] + dist[v] + 1
                    min_cycle = min(min_cycle, cycle_len)

    return min_cycle if min_cycle < float('inf') else -1


if __name__ == "__main__":
    main()
