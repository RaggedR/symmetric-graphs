# CLAUDE.md

## What This Is

Symmetry-aware graph drawing for cubic arc-transitive graphs using automorphism group quotients. Computes Aut(G) via GAP, selects optimal orbit quotient, places vertices on concentric rings with rotational symmetry, renders via TikZ/pdflatex.

## Stack

Python 3, GAP (with GRAPE package), LaTeX/TikZ, pdflatex.

## Commands

```bash
python3 draw_graph.py dodecahedron   # Single graph
python3 draw_all.py                  # Full gallery (all 13 Foster census graphs up to 30 vertices)
```

## Key Files

| File | Purpose |
|------|---------|
| `draw_graph.py` | Layout engine: GAP interface, crossing minimization, TikZ output |
| `draw_all.py` | Generates the gallery of all graphs |
| `symmetric-layout.g` | GAP script: Aut(G), orbit decompositions, quotient graphs |
| `gallery.pdf` | Rendered gallery |
| `SYMMETRIC_EMBEDDING.md` | Theory and implementation notes |

## Pipeline

1. GAP computes Aut(G), enumerates subgroups, finds orbit quotients
2. Python selects best quotient, detects cyclic vs cylindrical layout
3. Crossing minimization via path reordering + angular offset search
4. TikZ renders final PDF
