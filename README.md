# Symmetric Graph Drawing

Quotient-driven layouts for symmetric graphs, using automorphism group orbits to determine vertex placement on concentric rings with rotational symmetry.

## Key files

| File | Description |
|------|-------------|
| `gallery.pdf` | All 7 graphs on one page |
| `draw_graph.py` | Layout engine: GAP interface, crossing minimisation, TikZ output |
| `draw_all.py` | Generates the gallery |
| `symmetric-layout.g` | GAP script: computes Aut(G), orbit decompositions, quotient graphs |
| `SYMMETRIC_EMBEDDING.md` | Theory and implementation notes |
| `thesis.pdf` | Reference: Eades & Hong, *Symmetric Graph Drawing* |

## How it works

1. **GAP** computes Aut(G), enumerates subgroups, finds orbit quotients with low degree
2. **Python** selects the best quotient (by symmetry score), detects cyclic vs cylindrical layout, and places vertices on concentric rings
3. **Crossing minimisation**: path reordering of orbits (Schlegel-like), global angular offset search with rotations and reflections
4. **TikZ** renders the final PDF

## Graphs

All seven graphs are **cubic** (3-regular) — every vertex has degree 3. Cubic graphs are the simplest non-trivial regular graphs, and most of the classic symmetric graphs happen to be cubic: the Petersen graph, the generalized Petersen graphs GP(n,k), the Levi graphs of point-line configurations, and the polyhedral skeletons of the Platonic solids. They are a natural testbed for symmetry-aware drawing because their automorphism groups are large relative to their size.

Given a graph Γ, we compute G = Aut(Γ) and select a subgroup H ≤ G whose orbits give a clean quotient. The layout places each H-orbit on a concentric ring.

| Graph | n | G = Aut(Γ) | \|G\| | H | \|H\| | Orbits | Crossings |
|-------|:-:|------------|:-----:|---|:-----:|--------|:---------:|
| Petersen | 10 | S₅ | 120 | D₁₀ | 10 | [5, 5] | 5 |
| Heawood | 14 | PGL(2,7) | 336 | C₂ × C₂ | 4 | [2, 2, 2, 4, 4] | 14 |
| Pappus | 18 | (C₃ × C₃) ⋊ S₃ | 216 | S₃ | 6 | [6, 6, 6] | 7 |
| Dodecahedron | 20 | A₅ × C₂ | 120 | C₅ | 5 | [5, 5, 5, 5] | 0 |
| Desargues GP(10,3) | 20 | S₅ × C₂ | 240 | C₁₀ | 10 | [10, 10] | 20 |
| Mobius-Kantor GP(8,3) | 16 | GL(2,3) ⋊ C₂ | 96 | C₈ | 8 | [8, 8] | 16 |
| Cube | 8 | S₄ × C₂ | 48 | C₄ | 4 | [4, 4] | 0 |

## Usage

```bash
# Single graph
python3 draw_graph.py dodecahedron

# Full gallery
python3 draw_all.py
```

Requires [GAP](https://www.gap-system.org/) with the GRAPE package, Python 3, and pdflatex.
