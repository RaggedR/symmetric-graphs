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

| Graph | Vertices | Crossing number | Our layout |
|-------|:--------:|:---------------:|:----------:|
| Petersen | 10 | 2 | 5 |
| Heawood | 14 | ? | 14 |
| Pappus | 18 | 5 | 7 |
| Dodecahedron | 20 | 0 | 0 |
| Desargues | 20 | ? | 20 |
| Mobius-Kantor | 16 | ? | 16 |
| Cube | 8 | 0 | 0 |

## Usage

```bash
# Single graph
python3 draw_graph.py dodecahedron

# Full gallery
python3 draw_all.py
```

Requires [GAP](https://www.gap-system.org/) with the GRAPE package, Python 3, and pdflatex.
