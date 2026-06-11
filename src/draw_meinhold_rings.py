#!/usr/bin/env python3
"""Meinhold-72 multi-ring layout via cyclic symmetry.

Tries three cyclic ring structures (orbits exported by meinhold_ring_export.g
to /tmp/meinhold_orbits.json):

  C7:  10 rings of 7 + 2 fixed points at centre (7-fold symmetry)
  C8:  6 rings of 8 + 6 rings of 4 (8-fold)
  C12: 4 rings of 12 + 4 rings of 6 (12-fold)

Orbits of size s < k sit at angular positions offset + j*(k/s), preserving
the k-fold rotational symmetry. Offset/reflection search minimises edge
length, then local refinement minimises true crossings. Best structure by
crossing count wins.

Output: individual_graphs/meinhold72-rings.{tex,pdf}
"""

import json
import math
import pathlib
import random
import re
import subprocess

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
N = 72


def load_adjacency():
    text = (ROOT / "src" / "symmetric-layout.g").read_text()
    m = re.search(
        r'GraphLibrary\.meinhold72 := function\(\)\s*return rec\('
        r'name := "Meinhold-72", n := 72, adj := (\[.*?\])\);',
        text, re.S)
    adj = json.loads(m.group(1).replace(" ", ""))
    assert len(adj) == N and all(len(row) == 7 for row in adj)
    return adj


ADJ = load_adjacency()
EDGES = np.array([(v, u - 1) for v in range(N) for u in ADJ[v] if v < u - 1],
                 dtype=np.int32)
assert len(EDGES) == 252


# ============================================================
# Coordinates: k-fold layout with mixed orbit sizes s | k
# ============================================================

def compute_radii(n_rings):
    R_max, R_min = 7.5, 2.2
    if n_rings == 1:
        return [R_max]
    return [R_max - i * (R_max - R_min) / (n_rings - 1) for i in range(n_rings)]


def make_coords(rings, offsets, radii, k):
    """rings: list of generator-ordered orbits (1-indexed vertices).
    Size-1 orbits go to the centre."""
    coords = np.zeros((N, 2))
    centre = [r for r in rings if len(r) == 1]
    for ri, orb in enumerate(rings):
        s = len(orb)
        if s == 1:
            # spread fixed points slightly around centre
            ci = centre.index(orb)
            ang = 2 * math.pi * ci / max(1, len(centre))
            coords[orb[0] - 1] = (0.55 * math.cos(ang), 0.55 * math.sin(ang))
            continue
        step = k // s
        for j, v in enumerate(orb):
            theta = 2 * math.pi * ((offsets[ri] + j * step) % k) / k + math.pi / 2
            coords[v - 1] = (radii[ri] * math.cos(theta),
                             radii[ri] * math.sin(theta))
    return coords


def edge_length_sq(coords):
    diff = coords[EDGES[:, 0]] - coords[EDGES[:, 1]]
    return float((diff ** 2).sum())


def count_crossings(coords):
    crossings = 0
    p1 = coords[EDGES[:, 0]]
    p2 = coords[EDGES[:, 1]]
    n_edges = len(EDGES)
    for i in range(n_edges):
        ei = EDGES[i]
        rest = EDGES[i + 1:]
        mask = ((rest[:, 0] != ei[0]) & (rest[:, 0] != ei[1]) &
                (rest[:, 1] != ei[0]) & (rest[:, 1] != ei[1]))
        if not mask.any():
            continue
        ax, ay = p1[i]
        bx, by = p2[i]
        cx, cy = p1[i + 1:][mask].T
        dx, dy = p2[i + 1:][mask].T
        d1 = (dx - cx) * (ay - cy) - (dy - cy) * (ax - cx)
        d2 = (dx - cx) * (by - cy) - (dy - cy) * (bx - cx)
        d3 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        d4 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax)
        denom1, denom2 = d1 - d2, d3 - d4
        valid = (np.abs(denom1) > 1e-12) & (np.abs(denom2) > 1e-12)
        if not valid.any():
            continue
        t = np.where(valid, d1 / np.where(valid, denom1, 1.0), 0.0)
        s = np.where(valid, d3 / np.where(valid, denom2, 1.0), 0.0)
        EPS = 1e-9
        crossings += int((valid & (t > EPS) & (t < 1 - EPS)
                          & (s > EPS) & (s < 1 - EPS)).sum())
    return crossings


# ============================================================
# Search: reflections (sampled) + greedy offsets + local refinement
# ============================================================

def search_layout(orbits, k, n_mirror_samples=128):
    # big rings outermost: sort by size desc so dense rings get more room
    rings0 = sorted(orbits, key=len, reverse=True)
    n_rings = len(rings0)
    radii = compute_radii(n_rings)

    rotatable = [i for i, r in enumerate(rings0) if len(r) > 1]
    total_mirrors = 1 << max(0, len(rotatable) - 1)
    if n_mirror_samples >= total_mirrors:
        mirror_list = range(total_mirrors)
    else:
        mirror_list = [0] + [random.randint(1, total_mirrors - 1)
                             for _ in range(n_mirror_samples - 1)]

    best = (float("inf"), None, None)
    for mirror_bits in mirror_list:
        rings = [list(r) for r in rings0]
        for bit, ri in enumerate(rotatable[1:]):
            if mirror_bits & (1 << bit):
                rings[ri] = rings[ri][::-1]
        offsets = [0] * n_rings
        for ri in rotatable[1:]:
            best_off, best_score = 0, float("inf")
            for off in range(k):
                offsets[ri] = off
                score = edge_length_sq(make_coords(rings, offsets, radii, k))
                if score < best_score:
                    best_score, best_off = score, off
            offsets[ri] = best_off
        proxy = edge_length_sq(make_coords(rings, offsets, radii, k))
        if proxy < best[0]:
            best = (proxy, [list(r) for r in rings], list(offsets))

    _, rings, offsets = best
    best_cross = count_crossings(make_coords(rings, offsets, radii, k))
    improved = True
    while improved:
        improved = False
        for ri in rotatable:
            for off in range(k):
                if off == offsets[ri]:
                    continue
                trial = list(offsets)
                trial[ri] = off
                c = count_crossings(make_coords(rings, trial, radii, k))
                if c < best_cross:
                    best_cross, offsets = c, trial
                    improved = True
    coords = make_coords(rings, offsets, radii, k)
    return coords, rings, best_cross


# ============================================================
# TikZ (house palette, blue outer -> red inner)
# ============================================================

RING_PALETTE = [
    "0033CC", "0066FF", "0099FF", "00CCCC", "00CC66",
    "33CC00", "99CC00", "CCCC00", "FF9900", "FF6600",
    "FF3300", "CC0000",
]


def generate_tikz(coords, rings, subtitle, crossings):
    n_rings = len(rings)
    ring_of = {}
    for ri, orb in enumerate(rings):
        for v in orb:
            ring_of[v - 1] = ri

    def color_idx(ri):
        return min(11, int(ri * 12 / n_rings))

    lines = [
        r"\documentclass[border=12pt]{standalone}",
        r"\usepackage{tikz}",
        r"\usepackage{amsmath, amssymb}",
    ]
    for idx, hex_c in enumerate(RING_PALETTE):
        lines.append(f"\\definecolor{{ring{idx}}}{{HTML}}{{{hex_c}}}")
    lines += [r"\begin{document}", r"\begin{tikzpicture}"]

    y_title = coords[:, 1].max() + 1.6
    lines.append(f"\\node[font=\\Large\\bfseries] at (0, {y_title:.1f}) "
                 r"{Meinhold-72};")
    lines.append(f"\\node[font=\\small, text=gray] at (0, {y_title - 0.7:.1f}) "
                 f"{{{subtitle}}};")
    lines.append(f"\\node[font=\\footnotesize, text=gray!70] at "
                 f"(0, {y_title - 1.2:.1f}) {{{crossings} crossings}};")

    # Bipartite graph: every edge is inter-ring. Colour each edge by its
    # outer endpoint's ring so edge bundles read as coloured fans.
    for v, u in EDGES:
        rv, ru = ring_of[int(v)], ring_of[int(u)]
        outer = min(rv, ru)  # ring 0 is outermost
        x1, y1 = coords[v]
        x2, y2 = coords[u]
        lines.append(f"\\draw[line width=0.35pt, ring{color_idx(outer)}, "
                     f"opacity=0.40] ({x1:.3f}, {y1:.3f}) -- "
                     f"({x2:.3f}, {y2:.3f});")

    for v in range(N):
        x, y = coords[v]
        c = f"ring{color_idx(ring_of[v])}"
        lines.append(f"\\fill[{c}] ({x:.3f}, {y:.3f}) circle (2.6pt);")

    lines += [r"\end{tikzpicture}", r"\end{document}"]
    return "\n".join(lines)


def compile_pdf(tex_path):
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", tex_path.name],
        capture_output=True, text=True, cwd=tex_path.parent)
    if result.returncode != 0:
        print(result.stdout[-600:])
        raise SystemExit(f"pdflatex failed for {tex_path}")
    for ext in (".aux", ".log"):
        p = tex_path.with_suffix(ext)
        if p.exists():
            p.unlink()


if __name__ == "__main__":
    random.seed(42)
    orbit_sets = json.load(open("/tmp/meinhold_orbits.json"))

    results = {}
    for key, k in [("c7", 7), ("c8", 8), ("c12", 12)]:
        print(f"=== {key}: {len(orbit_sets[key])} rings, {k}-fold ===")
        coords, rings, crossings = search_layout(orbit_sets[key], k)
        print(f"  crossings: {crossings}")
        results[key] = (coords, rings, crossings, k)

    best_key = min(results, key=lambda kk: results[kk][2])
    coords, rings, crossings, k = results[best_key]
    print(f"\nBest: {best_key} with {crossings} crossings")

    n_big = sum(1 for r in rings if len(r) > 1)
    subtitle = (r"$\mathrm{Aut} \cong \mathrm{PSU}(3,3) \rtimes C_2$, "
                r"$|\mathrm{Aut}| = 12096$, 72 vertices, 7-regular --- "
                f"$C_{{{k}}}$ rings")
    tex = generate_tikz(coords, rings, subtitle, crossings)
    out = ROOT / "individual_graphs" / "meinhold72-rings.tex"
    out.write_text(tex)
    compile_pdf(out)
    print(f"written {out} (+pdf)")
