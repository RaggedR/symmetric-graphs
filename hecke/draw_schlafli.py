"""
Draw the Schläfli graph — SRG(27, 10, 1, 5)
Layout: 2₂₁ polytope projection — 2 concentric rings of 12 + 3 at center.
Uses an order-12 automorphism from Aut(Schläfli) = W(E₆).
"""

import numpy as np
import math
import subprocess, shutil

# Adjacency list (1-indexed)
ADJ_1 = [
    [2, 3, 4, 7, 10, 13, 16, 20, 23, 26],
    [1, 3, 5, 8, 11, 14, 17, 21, 24, 27],
    [1, 2, 6, 9, 12, 15, 18, 19, 22, 25],
    [1, 5, 6, 7, 12, 15, 18, 21, 24, 27],
    [2, 4, 6, 8, 10, 13, 16, 19, 22, 25],
    [3, 4, 5, 9, 11, 14, 17, 20, 23, 26],
    [1, 4, 8, 9, 11, 14, 17, 19, 22, 25],
    [2, 5, 7, 9, 12, 15, 18, 20, 23, 26],
    [3, 6, 7, 8, 10, 13, 16, 21, 24, 27],
    [1, 5, 9, 11, 12, 14, 18, 19, 23, 27],
    [2, 6, 7, 10, 12, 15, 16, 20, 24, 25],
    [3, 4, 8, 10, 11, 13, 17, 21, 22, 26],
    [1, 5, 9, 12, 14, 15, 17, 20, 24, 25],
    [2, 6, 7, 10, 13, 15, 18, 21, 22, 26],
    [3, 4, 8, 11, 13, 14, 16, 19, 23, 27],
    [1, 5, 9, 11, 15, 17, 18, 21, 22, 26],
    [2, 6, 7, 12, 13, 16, 18, 19, 23, 27],
    [3, 4, 8, 10, 14, 16, 17, 20, 24, 25],
    [3, 5, 7, 10, 15, 17, 20, 21, 24, 26],
    [1, 6, 8, 11, 13, 18, 19, 21, 22, 27],
    [2, 4, 9, 12, 14, 16, 19, 20, 23, 25],
    [3, 5, 7, 12, 14, 16, 20, 23, 24, 27],
    [1, 6, 8, 10, 15, 17, 21, 22, 24, 25],
    [2, 4, 9, 11, 13, 18, 19, 22, 23, 26],
    [3, 5, 7, 11, 13, 18, 21, 23, 26, 27],
    [1, 6, 8, 12, 14, 16, 19, 24, 25, 27],
    [2, 4, 9, 10, 15, 17, 20, 22, 25, 26],
]

N = 27
adj = [[v - 1 for v in row] for row in ADJ_1]

# Order-12 automorphism orbits (from GAP, converted to 0-indexed):
# Outer ring (12): [0, 2, 8, 26, 3, 17, 15, 21, 6, 24, 20, 19]
# Inner ring (12): [1, 5, 9, 14, 23, 16, 4, 13, 10, 22, 18, 12]
# Center (3):      [7, 25, 11]

outer = [0, 2, 8, 26, 3, 17, 15, 21, 6, 24, 20, 19]  # 1-indexed: 1,3,9,27,4,18,16,22,7,25,21,20
inner = [1, 5, 9, 14, 23, 16, 4, 13, 10, 22, 18, 12]  # 1-indexed: 2,6,10,15,24,17,5,14,11,23,19,13
center = [7, 25, 11]                                    # 1-indexed: 8,26,12

# Positions
W, H_svg = 900, 900
cx, cy = W / 2, H_svg / 2
r_outer = 320
r_inner = 190
r_center = 0  # all at center, but slightly spread

pos = np.zeros((N, 2))

# Try different angular offsets between rings to minimise crossings
best_score = float('inf')
best_offset = 0

for offset_deg in range(30):
    offset = offset_deg * math.pi / 180
    test_pos = np.zeros((N, 2))
    for k, v in enumerate(outer):
        angle = 2 * math.pi * k / 12 - math.pi / 2
        test_pos[v] = [cx + r_outer * math.cos(angle), cy + r_outer * math.sin(angle)]
    for k, v in enumerate(inner):
        angle = 2 * math.pi * k / 12 - math.pi / 2 + offset
        test_pos[v] = [cx + r_inner * math.cos(angle), cy + r_inner * math.sin(angle)]
    for k, v in enumerate(center):
        angle = 2 * math.pi * k / 3 - math.pi / 2
        test_pos[v] = [cx + 35 * math.cos(angle), cy + 35 * math.sin(angle)]

    # Score: sum of inverse distances (penalize short edges = crossings)
    score = 0
    for i in range(N):
        for j in adj[i]:
            if j > i:
                dx = test_pos[i][0] - test_pos[j][0]
                dy = test_pos[i][1] - test_pos[j][1]
                d = math.sqrt(dx*dx + dy*dy)
                if d < 50:
                    score += 1000
                score += 1.0 / max(d, 1)

    if score < best_score:
        best_score = score
        best_offset = offset

# Apply best offset
for k, v in enumerate(outer):
    angle = 2 * math.pi * k / 12 - math.pi / 2
    pos[v] = [cx + r_outer * math.cos(angle), cy + r_outer * math.sin(angle)]
for k, v in enumerate(inner):
    angle = 2 * math.pi * k / 12 - math.pi / 2 + best_offset
    pos[v] = [cx + r_inner * math.cos(angle), cy + r_inner * math.sin(angle)]
for k, v in enumerate(center):
    angle = 2 * math.pi * k / 3 - math.pi / 2
    pos[v] = [cx + 35 * math.cos(angle), cy + 35 * math.sin(angle)]

px = pos[:, 0]
py = pos[:, 1]

# --- SVG ---
vr = 8
lines = []
lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H_svg}">')

# Background
lines.append(f'  <rect width="{W}" height="{H_svg}" fill="#0f0f1a"/>')

# Title
lines.append(f'  <text x="{W/2}" y="35" text-anchor="middle" '
             f'font-family="Georgia, serif" font-size="22" fill="#e8e8e8" '
             f'letter-spacing="1.5">The Schläfli Graph</text>')
lines.append(f'  <text x="{W/2}" y="56" text-anchor="middle" '
             f'font-family="Georgia, serif" font-size="13" fill="#7799aa">'
             f'Sab(PSp(4,3), (C₂⁴):A₅, HaH)  —  27 lines on a cubic surface</text>')

# Classify edges by ring membership for layered drawing
ring_of = {}
for v in outer:
    ring_of[v] = 'outer'
for v in inner:
    ring_of[v] = 'inner'
for v in center:
    ring_of[v] = 'center'

# Draw edges
lines.append('  <g>')
for i in range(N):
    for j in adj[i]:
        if j <= i:
            continue
        ri, rj = ring_of[i], ring_of[j]

        if ri == rj == 'outer':
            color, opacity, width = "#5599bb", 0.7, 1.0
        elif ri == rj == 'inner':
            color, opacity, width = "#66bb99", 0.7, 1.0
        elif ri == rj == 'center':
            color, opacity, width = "#dddddd", 0.8, 1.2
        elif 'center' in (ri, rj):
            color, opacity, width = "#9988aa", 0.6, 0.9
        else:
            # outer-inner
            color, opacity, width = "#7799aa", 0.55, 0.8

        lines.append(
            f'    <line x1="{px[i]:.1f}" y1="{py[i]:.1f}" '
            f'x2="{px[j]:.1f}" y2="{py[j]:.1f}" '
            f'stroke="{color}" stroke-width="{width}" opacity="{opacity}"/>'
        )
lines.append('  </g>')

# Draw vertices
for i in range(N):
    r = ring_of[i]
    if r == 'outer':
        fill = "#4ea8de"
    elif r == 'inner':
        fill = "#48c9a0"
    else:
        fill = "#e8e8e8"

    lines.append(
        f'  <circle cx="{px[i]:.1f}" cy="{py[i]:.1f}" r="{vr}" '
        f'fill="{fill}" stroke="#ffffff" stroke-width="1.0" opacity="0.92"/>'
    )
    lines.append(
        f'  <text x="{px[i]:.1f}" y="{py[i] + 3.2:.1f}" text-anchor="middle" '
        f'font-family="Helvetica, sans-serif" font-size="8" fill="#0f0f1a" '
        f'font-weight="bold">{i + 1}</text>'
    )

# Legend
ly = H_svg - 55
lines.append(f'  <text x="18" y="{ly + 14}" font-family="Georgia, serif" '
             f'font-size="11" fill="#7799aa">27 vertices · 135 edges · 10-regular · λ = 1, μ = 5</text>')
lines.append(f'  <text x="18" y="{ly + 32}" font-family="Georgia, serif" '
             f'font-size="11" fill="#7799aa">Layout: 2₂₁ polytope projection (order-12 automorphism)</text>')

lines.append(f'  <text x="{W - 15}" y="{H_svg - 12}" text-anchor="end" '
             f'font-family="Georgia, serif" font-size="11" fill="#445566">'
             f'Aut = W(E₆) = O(5,3):C₂, order 51840</text>')

lines.append('</svg>')

with open('schlafli.svg', 'w') as f:
    f.write('\n'.join(lines))

# Convert to PNG then JPEG
try:
    import cairosvg
    cairosvg.svg2png(url='schlafli.svg', write_to='schlafli.png',
                     output_width=1800, output_height=1800)
    subprocess.run(['sips', '-s', 'format', 'jpeg', 'schlafli.png',
                    '--out', 'schlafli.jpg'], capture_output=True)
    print("Written schlafli.svg, schlafli.jpg (1800×1800)")
except ImportError:
    print("Written schlafli.svg (install cairosvg for JPEG)")
