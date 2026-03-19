# Symmetric Graph Drawing

Quotient-driven layouts for symmetric graphs, using automorphism group orbits to determine vertex placement on concentric rings with rotational symmetry.

## Key files

| File | Description |
|------|-------------|
| `gallery.pdf` | All 13 graphs on one page |
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

All graphs are **cubic** (3-regular). The collection is the complete [Foster census](https://en.wikipedia.org/wiki/Foster_census) of cubic arc-transitive (symmetric) graphs up to 20 vertices, plus selected entries up to 30. Cubic graphs are the simplest non-trivial regular graphs, and their automorphism groups are large relative to their size, making them a natural testbed for symmetry-aware drawing.

**The list is exhaustive for n ≤ 20**: every cubic arc-transitive graph on at most 20 vertices appears in the gallery.

Given a graph Γ, we compute G = Aut(Γ) and select a subgroup H ≤ G whose orbits give a clean quotient. The layout places each H-orbit on a concentric ring.

| Foster | Graph | n | G = Aut(Γ) | \|G\| | H | \|H\| | Orbits |
|:------:|-------|:-:|------------|:-----:|---|:-----:|--------|
| F004 | K₄ | 4 | S₄ | 24 | C₂ | 2 | [2, 2] |
| F006 | K₃,₃ | 6 | S₃ ≀ C₂ | 72 | C₂ | 2 | [3, 3] |
| F008 | Cube | 8 | S₄ × C₂ | 48 | C₄ | 4 | [4, 4] |
| F010 | Petersen | 10 | S₅ | 120 | D₁₀ | 10 | [5, 5] |
| F014 | Heawood | 14 | PGL(2,7) | 336 | C₂ × C₂ | 4 | [2, 2, 2, 4, 4] |
| F016 | Mobius-Kantor GP(8,3) | 16 | GL(2,3) ⋊ C₂ | 96 | C₈ | 8 | [8, 8] |
| F018 | Pappus | 18 | (C₃ × C₃) ⋊ S₃ | 216 | S₃ | 6 | [6, 6, 6] |
| F020A | Dodecahedron | 20 | A₅ × C₂ | 120 | C₅ | 5 | [5, 5, 5, 5] |
| F020B | Desargues GP(10,3) | 20 | S₅ × C₂ | 240 | C₁₀ | 10 | [10, 10] |
| F024 | Nauru GP(12,5) | 24 | S₄ × D₁₂ | 144 | C₁₂ | 12 | [12, 12] |
| F026A | F26A | 26 | C₁₃ ⋊ C₆ | 78 | — | — | — |
| F028 | Coxeter | 28 | PGL(2,7) | 336 | D₁₄ | 14 | [7, 7, 7, 7] |
| F030 | Tutte-Coxeter | 30 | S₆ × C₂ | 1440 | — | — | — |

Note: no cubic arc-transitive graph exists on 12 vertices.

## Usage

```bash
# Single graph
python3 draw_graph.py dodecahedron

# Full gallery
python3 draw_all.py
```

Requires [GAP](https://www.gap-system.org/) with the GRAPE package, Python 3, and pdflatex.
