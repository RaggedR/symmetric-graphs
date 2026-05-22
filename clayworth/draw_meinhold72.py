#!/usr/bin/env python3
"""Draw Meinhold-72 with multi-ring layout based on GAP orbit decomposition."""

import subprocess
import math
import os
import tempfile

N = 72
adj = [
    [4,24,32,38,43,52,59],[10,15,17,20,33,44,69],[10,11,20,29,34,37,47],
    [1,14,18,28,48,65,71],[11,22,25,35,37,56,66],[7,15,22,33,38,59,66],
    [6,8,23,39,41,54,57],[7,9,22,24,25,32,53],[8,12,23,26,58,61,64],
    [2,3,31,39,40,51,57],[3,5,14,36,40,65,72],[9,25,30,44,56,60,69],
    [26,45,49,50,61,67,70],[4,11,16,20,43,56,69],[2,6,18,45,57,70,71],
    [14,18,26,36,45,46,64],[2,18,19,23,39,48,64],[4,15,16,17,21,38,63],
    [17,21,30,34,44,47,68],[2,3,14,45,48,49,62],[18,19,27,28,46,50,70],
    [5,6,8,36,45,61,62],[7,9,17,30,33,55,63],[1,8,31,39,48,58,62],
    [5,8,12,31,57,65,71],[9,13,16,32,43,60,63],[21,35,38,42,52,66,68],
    [4,21,30,55,56,59,66],[3,41,49,51,54,67,72],[12,19,23,28,54,65,72],
    [10,24,25,34,35,44,52],[1,8,26,36,54,65,67],[2,6,23,40,49,61,72],
    [3,19,31,48,50,65,67],[5,27,31,40,50,58,61],[11,16,22,32,38,47,68],
    [3,5,41,45,46,50,57],[1,6,18,27,36,39,40],[7,10,17,24,38,42,47],
    [10,11,33,35,38,43,63],[7,29,37,42,53,55,66],[27,39,41,46,51,58,64],
    [1,14,26,40,49,51,58],[2,12,19,31,61,62,70],[13,15,16,20,22,37,53],
    [16,21,37,42,47,56,60],[3,19,36,39,46,54,62],[4,17,20,24,34,53,55],
    [13,20,29,33,43,55,59],[13,21,34,35,37,55,63],[10,29,42,43,52,60,69],
    [1,27,31,51,67,70,71],[8,41,45,48,64,67,71],[7,29,30,32,47,59,60],
    [23,28,41,48,49,50,58],[5,12,14,28,46,58,62],[7,10,15,25,37,60,63],
    [9,24,35,42,43,55,56],[1,6,28,49,54,62,70],[12,26,46,51,54,57,70],
    [9,13,22,33,35,44,68],[20,22,24,44,47,56,59],[18,23,26,40,50,57,65],
    [9,16,17,42,53,68,69],[4,11,25,30,32,34,63],[5,6,27,28,41,71,72],
    [13,29,32,34,52,53,68],[19,27,36,61,64,67,72],[2,12,14,51,64,71,72],
    [13,15,21,44,52,59,60],[4,15,25,52,53,66,69],[11,29,30,33,66,68,69],
]

# Write GAP script to find orbit decomposition
gap_script = """LoadPackage("grape");;
adj := [
""" + ',\n'.join(str(a) for a in adj) + """
];;
gamma := Graph(Group(()), [1..72],
    function(x,g) return x; end,
    function(x,y) return y in adj[x]; end, true);;
aut := AutGroupGraph(gamma);;

# Find subgroup giving 6-9 equal-sized orbits
best := fail;; bestScore := 100;;
for cc in ConjugacyClassesSubgroups(aut) do
    rep := Representative(cc);;
    orbs := Orbits(rep, [1..72]);;
    sizes := SortedList(List(orbs, Length));;
    nOrbs := Length(orbs);;
    if nOrbs >= 5 and nOrbs <= 12 then
        maxSz := Maximum(sizes);; minSz := Minimum(sizes);;
        score := (maxSz - minSz) + AbsInt(nOrbs - 8);;
        if score < bestScore then
            bestScore := score; best := rec(H := rep, orbs := orbs);
        fi;
    fi;
od;

if best <> fail then
    orbs := best.orbs;;
    SortBy(orbs, function(o) return -Length(o); end);;
    Print("ORBITS:", List(orbs, function(o) return o; end), "\\n");;
    Print("SIZES:", List(orbs, Length), "\\n");;
    Print("H:", StructureDescription(best.H), "\\n");;
fi;
QuitGap(0);
"""

with tempfile.NamedTemporaryFile(mode='w', suffix='.g', delete=False) as f:
    f.write(gap_script)
    gap_file = f.name

print("Finding orbit decomposition via GAP...")
result = subprocess.run(['gap', '-q', gap_file], capture_output=True, text=True, timeout=120)
os.unlink(gap_file)

# Parse orbits
orbits = []
for line in result.stdout.split('\n'):
    if line.startswith('ORBITS:'):
        # Parse nested GAP list
        s = line[7:].strip()
        # Manual parse: split on ], [
        depth = 0
        current = []
        for ch in s:
            if ch == '[':
                depth += 1
                if depth == 2:
                    current = []
            elif ch == ']':
                if depth == 2 and current:
                    orbits.append(current)
                    current = []
                depth -= 1
            elif ch.isdigit():
                if not current or not isinstance(current[-1], str) or not current[-1].isdigit():
                    current.append(ch)
                else:
                    current[-1] += ch
            elif ch == ',' and depth == 2:
                pass  # separator
            elif ch == ' ':
                if current and isinstance(current[-1], str) and current[-1].isdigit():
                    current[-1] = int(current[-1])
                    current.append('')
        # Convert strings to ints
        orbits = [[int(x) if isinstance(x, str) and x else x for x in o if x] for o in orbits]
        orbits = [[x for x in o if isinstance(x, int)] for o in orbits if o]
    elif line.startswith('H:'):
        print(f"  Subgroup: {line[2:]}")
    elif line.startswith('SIZES:'):
        print(f"  Orbit sizes: {line[6:]}")

if not orbits:
    # Fallback: BFS-based rings from vertex 1
    print("  Falling back to BFS layout")
    dist = [-1] * (N + 1)
    dist[1] = 0
    queue = [1]
    while queue:
        v = queue.pop(0)
        for w in adj[v-1]:
            if dist[w] < 0:
                dist[w] = dist[v] + 1
                queue.append(w)
    max_d = max(dist[1:])
    orbits = [[] for _ in range(max_d + 1)]
    for v in range(1, N + 1):
        orbits[dist[v]].append(v)

print(f"  {len(orbits)} rings: {[len(o) for o in orbits]}")

# Place vertices on concentric rings
positions = {}
n_rings = len(orbits)
for ring_idx, orbit in enumerate(orbits):
    # Radius: outermost = 5.0, innermost = 1.0
    if n_rings > 1:
        r = 5.0 - ring_idx * 4.0 / (n_rings - 1)
    else:
        r = 3.0
    n_in_ring = len(orbit)
    for j, v in enumerate(orbit):
        angle = 2 * math.pi * j / n_in_ring - math.pi / 2
        positions[v] = (r * math.cos(angle), r * math.sin(angle))

# Generate TikZ
lines = []
lines.append(r"""\documentclass[border=10pt]{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[
    vertex/.style={circle, fill=black, inner sep=1pt, minimum size=3pt},
    edge/.style={draw=gray!40, line width=0.2pt}
]""")

lines.append(r"\node[above, font=\Large\bfseries] at (0,6.2) {Meinhold-72};")
lines.append(r"\node[above, font=\small] at (0,5.7) {72 vertices, 7-regular, girth 4, diameter 5};")
lines.append(r"\node[above, font=\small] at (0,5.3) {$\mathrm{Aut} = G_2(2)$, quotient of Clayworth via $K = \mathrm{PSL}(3,2)$};")

# Edges
for v in range(1, N + 1):
    x1, y1 = positions[v]
    for w in adj[v - 1]:
        if w > v:
            x2, y2 = positions[w]
            lines.append(f"  \\draw[edge] ({x1:.4f},{y1:.4f}) -- ({x2:.4f},{y2:.4f});")

# Vertices
for v in range(1, N + 1):
    x, y = positions[v]
    lines.append(f"  \\node[vertex] at ({x:.4f},{y:.4f}) {{}};")

lines.append(r"\end{tikzpicture}")
lines.append(r"\end{document}")

tex_path = "named_graphs/meinhold72.tex"
os.makedirs("named_graphs", exist_ok=True)
with open(tex_path, 'w') as f:
    f.write('\n'.join(lines))

print(f"Compiling {tex_path}...")
result = subprocess.run(
    ['pdflatex', '-interaction=nonstopmode', '-output-directory=named_graphs', tex_path],
    capture_output=True, text=True, timeout=30
)

pdf_path = "named_graphs/meinhold72.pdf"
if os.path.exists(pdf_path):
    print(f"Created {pdf_path}")
    subprocess.run(['open', pdf_path])
else:
    print("pdflatex failed")
    for line in result.stdout.split('\n'):
        if '!' in line or 'Error' in line:
            print(f"  {line}")
