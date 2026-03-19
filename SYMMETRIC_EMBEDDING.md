# Symmetric Graph Embedding Tool

## What This Is

A pipeline for drawing symmetric (arc-transitive) graphs with visible symmetry.
GAP computes the group theory, Python computes coordinates, output is TikZ PDF or interactive D3.js HTML.

Robin's honours thesis was on "Symmetric Graphs and Their Quotients" — this tool automates the drawing.

## Quick Start

```bash
cd ~/git/graphviz

# Single graph → PDF (auto-detects best layout)
python3 draw_graph.py petersen
python3 draw_graph.py dodecahedron --layout rings
python3 draw_graph.py heawood --layout cyclic

# Single graph → interactive HTML (draggable vertices)
python3 draw_graph_d3.py pappus

# All graphs → single gallery PDF
python3 draw_all.py
```

## Files

| File | Purpose |
|------|---------|
| `symmetric-layout.g` | **GAP script**: computes Aut(Γ), enumerates subgroups, finds orbit quotients, outputs JSON |
| `draw_graph.py` | **Python → PDF**: reads GAP JSON, auto-detects cyclic vs cylindrical, outputs TikZ, compiles |
| `draw_graph_d3.py` | **Python → HTML**: same pipeline but outputs interactive D3.js with draggable vertices |
| `draw_all.py` | **Gallery**: draws all graphs in one PDF (calls draw_graph.py internally) |
| `spectral_layout.py` | Earlier attempt: eigenvector embedding. Abandoned — eigenvalue multiplicity kills it |

### Supporting files
| File | Purpose |
|------|---------|
| `dodecahedron-quotient.g` | GAP proof: dodecahedron antipodal quotient ≅ Petersen |
| `z5-orbits.g` | Detailed Z₅ orbit analysis of dodecahedron |
| `z5-dodecahedron.tex` | Hand-crafted Z₅ concentric pentagon layout (reference quality) |
| `dodec-to-petersen.tex` | Side-by-side dodecahedron → Petersen quotient diagram |
| `nftmarket-arch.dot` | Unrelated: NFTmarket architecture diagram |

## Pipeline Architecture

```
symmetric-layout.g (GAP)          draw_graph.py (Python)
┌─────────────────────┐          ┌──────────────────────┐
│ 1. Build graph      │          │ 1. Call GAP           │
│ 2. Compute Aut(Γ)   │  JSON   │ 2. Parse JSON         │
│ 3. Enumerate subs   │ ──────> │ 3. Select best quot   │
│ 4. Compute quotients│          │ 4. Detect layout type │
│ 5. Analyse offsets  │          │ 5. Compute coords     │
│ 6. Output JSON      │          │ 6. Generate TikZ      │
└─────────────────────┘          │ 7. pdflatex → open    │
                                 └──────────────────────┘
```

### How GAP is called
Python runs `gap -q` with stdin:
```
graph_name := "petersen";
max_deg := 4;
Read("/Users/robin/git/graphviz/symmetric-layout.g");
```
GAP outputs JSON between `=== JSON OUTPUT ===` and `=== END JSON ===` markers.
Python parses with `eval()` (GAP list syntax ≈ Python list syntax).

## Graph Library (defined in `symmetric-layout.g`)

All adjacency lists verified against known |Aut| values.

| Graph | n | |Aut| | Best Layout | Detected As |
|-------|---|------|-------------|-------------|
| petersen | 10 | 120 (S₅) | 2 rings | cylindrical (GP(5,2)) |
| heawood | 14 | 336 (PSL(3,2):C₂) | single circle | cyclic |
| pappus | 18 | 216 ((C₃×C₃):C₃):D₈ | 3 hexagonal rings | cylindrical |
| dodecahedron | 20 | 120 (A₅×Z₂) | 4 pentagonal rings | cylindrical |
| desargues | 20 | 240 (C₂×S₅) | 2 rings | cylindrical (GP(10,3)) |
| moebiuskantor | 16 | 96 (GL(2,3):C₂) | 2 rings | cylindrical (GP(8,3)) |
| cube | 8 | 48 (C₂×S₄) | varies | cylindrical |
| icosahedron | 12 | 120 (A₅×Z₂) | varies | cylindrical |
| k33 | 6 | 72 | varies | — |

### Adding a new graph
Add to `GraphLibrary` in `symmetric-layout.g`:
```gap
GraphLibrary.newgraph := function()
    return rec(
        name := "New Graph",
        n := 10,
        adj := [ [2,5,6], [1,3,7], ... ]  # 1-indexed adjacency lists
    );
end;
```
**Always verify** |Aut| matches the known value. A single wrong edge kills symmetry.

### Building from LCF notation
Many cubic symmetric graphs have LCF notation. To convert:
```gap
lcf := [5,7,-7,7,-7,-5, ...];  # repeated pattern
adj := List([1..n], i -> []);
for i in [1..n] do
    j := (i mod n) + 1;              # Hamiltonian cycle edge
    AddSet(adj[i], j); AddSet(adj[j], i);
    k := ((i-1+lcf[i]) mod n) + 1;   # LCF chord
    AddSet(adj[i], k); AddSet(adj[k], i);
od;
```

## Theory: Eades & Hong (Graph Drawing Handbook, Ch. 3)

### Theorem 3.2 — Which automorphism groups are "geometric"?

An automorphism group A can be displayed as a 2D symmetry iff:

- **(a) Reflection**: |A|=2 and fix_A induces disjoint paths
- **(b) Rotation**: A=⟨ρ⟩ cyclic, |fix_A|≤1, and if A fixes an edge then |fix_A|=0
- **(c) Dihedral**: A=⟨α,ρ⟩ dihedral, |fix_A|≤1, fix_α induces disjoint paths

### The Rotation Algorithm (Theorem 3.2b)
Given cyclic A=⟨ρ⟩ of order k, with m=n/k orbits O₁,...,Oₘ:
1. Choose uᵢ from each orbit Oᵢ
2. Place ρʲ(uᵢ) at angle **2π(i + j·m)/n**
3. ALL vertices on ONE circle, orbits interleaved
4. Rotation by 2π/k is a geometric symmetry of the drawing

### NP-completeness
Finding the **maximum** geometric automorphism group is NP-complete (Thm 3.3).
But **given** a specific automorphism, drawing it is O(n) (Corollary 3.2).

### Reference
- Eades & Hong, "Symmetric Graph Drawing", Graph Drawing Handbook Ch. 3
  (PDF: cs.brown.edu/people/rtamassi/gdhandbook/chapters/symmetry.pdf)
- Carr & Kocay, "An Algorithm for Drawing a Graph Symmetrically" (1999)

## Two Layout Modes

### Cyclic (Eades-Hong rotation) — `--layout cyclic`
- All vertices on ONE circle
- Orbits interleaved: orbit i at positions i, i+m, i+2m, ...
- Shows Z_k rotation as geometric rotation by 2π/k
- Good for: Heawood, any graph where interleaving looks clean

### Cylindrical (concentric rings) — `--layout rings`
- Each orbit on its own concentric circle
- Angular offsets between rings computed from cross-edge analysis
- Good for: generalized Petersen graphs (2 rings), dodecahedron (4 rings), Pappus (3 rings)

### Auto-detection (`detect_symmetry_type` in draw_graph.py)
- 2 equal orbits → cylindrical
- Quotient is a path → cylindrical
- ≤4 equal orbits, max degree ≤2 → cylindrical
- Otherwise → cyclic

## Scoring Function (how the best quotient is selected)

In `draw_graph.py`, `score_quotient()`:
- Quotient is a cycle: +100
- 2 equal orbits with internal cycle (GP pattern): **+150** (highest priority)
- Equal orbit sizes: +30
- Low max degree: -10 per degree
- 3-10 orbits: +20
- Internal cycles per orbit: +5 each

## Key Design Decisions & Lessons

1. **Spectral embedding fails** for highly symmetric graphs because eigenvalue multiplicity
   makes the basis choice arbitrary. More symmetry = worse spectral layout.

2. **Concentric rings are NOT the Eades-Hong rotation layout.** Rotation puts ALL orbits
   on one circle, interleaved. Concentric rings are for the dihedral case. We use concentric
   rings anyway for GP graphs because it matches the folk drawing convention.

3. **The angular offset between concentric rings** is computed by GAP's `AnalyseOffsets`:
   for each pair of adjacent orbits, it checks which positions are connected by cross-edges
   and finds the modal offset. A half-step offset (e.g., 36° for pentagons) means each
   outer vertex connects to the inner vertex at its position AND the adjacent position.

4. **GAP generator finding**: must pick the **highest-order** orbit-preserving element,
   not just any element. Low-order elements (reflections) don't trace full orbit cycles,
   causing `has_internal_cycle` to report false negatives.

5. **Always verify |Aut|** when adding a graph. The Pappus graph had two wrong edges
   (|Aut| was 1 instead of 216). Use LCF notation for cubic symmetric graphs.

6. **The "right" drawing depends on which symmetry you want to show.** The dodecahedron
   can be drawn with Z₅ (4 concentric pentagons), Z₂ (Petersen quotient), or Z₁₀ (2 rings).
   Each is valid but shows different structure. No algorithm can universally pick the "best" one.

7. **Buchheim & Hong's crossing minimization** is the key to getting angular offsets right
   for concentric ring layouts. But for graphs where all inter-orbit edges are "parallel spokes"
   (like GP graphs), the crossing count is always 0 — the spoke deviation metric is what
   actually determines the right offset. This is our addition to their algorithm.

8. **The Pappus graph** looks best with Eades-Hong cyclic layout (single circle, interleaved),
   NOT concentric rings. The auto-detection picks cylindrical because it has 3 equal orbits
   with a cycle quotient, but the orbit ordering from GAP doesn't produce clean hexagons.
   Override with `--layout cyclic`.

### Reference: Buchheim & Hong
"Crossing Minimization for Symmetries", Christoph Buchheim & Seok-Hee Hong, 2002.
PDF: kups.ub.uni-koeln.de/54867/1/zaik2002-440.pdf
Key results: Theorem 4 (O(m log m) for SCM+2), Corollary 1 (path orbit graph).

## Buchheim & Hong Crossing Minimization — IMPLEMENTED

Reference: "Crossing Minimization for Symmetries" (2002), Buchheim & Hong.

For concentric ring layouts, the key parameter is the **angular offset** (exponent)
between adjacent rings. The algorithm:

1. For each pair of adjacent orbit rings, try all k possible offsets (0 to k-1 steps)
2. **Primary metric**: minimize inter-orbit edge crossings (Theorem 4, O(m log m))
3. **Secondary metric** (our addition): minimize spoke deviation from radial —
   sum of squared angular distances between each spoke's endpoints
4. Pick the offset with (min crossings, min deviation)

This solved the Petersen inner-ring rotation problem: offset=3 steps (216°) makes
spokes nearly radial, giving the classic pentagon+pentagram drawing.

Results after implementing:
- Petersen: offset=3 (216°) — correct!
- Desargues GP(10,3): offset=6 (216°) — correct!
- Dodecahedron Z₅: rings 1-2 offset=0, ring 3 offset=4 (288°)
- Möbius-Kantor: offset=0

## TODO / Known Issues

1. **Pappus PDF** uses `--layout cyclic` (Eades-Hong single circle) because the
   cylindrical 3-ring layout has wrong orbit ordering from GAP. The interactive
   D3 Pappus looks great with cyclic layout. Need to make `draw_graph.py` auto-detect
   that Pappus should use cyclic not cylindrical.

2. **D3.js version** (`draw_graph_d3.py`) is out of sync with `draw_graph.py` —
   doesn't have the Buchheim & Hong crossing minimization or the GP orbit swap.
   Should be unified.

3. **Gallery PDF** (`draw_all.py`) works but layout/spacing may need tuning.

4. **No dihedral layout implemented** — Eades-Hong Theorem 3.2(c) describes
   circular grid for dihedral groups. We approximate with concentric rings.

5. Could add: Coxeter graph (28 vertices), or let user provide adjacency list via file.

## Dependencies

- **GAP 4.15** with GRAPE and Digraphs: `brew install gap-system/gap/gap`
- **Python 3** (no pip dependencies for draw_graph.py)
- **pdflatex** with tikz package
- **Graphviz** (dot, neato) for .dot files only
