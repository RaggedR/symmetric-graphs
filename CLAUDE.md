# CLAUDE.md

## What This Is

Symmetry-aware graph drawing for cubic arc-transitive graphs using automorphism group quotients. Computes Aut(G) via GAP, selects optimal orbit quotient, places vertices on concentric rings with rotational symmetry, renders via TikZ/pdflatex.

## Stack

Python 3, GAP (with GRAPE package), LaTeX/TikZ, pdflatex.

## Commands

```bash
python3 src/draw_graph.py dodecahedron   # Single graph → individual_graphs/
python3 src/draw_all.py                  # Gallery → gallery.pdf (11 graphs; K4 and K3,3 omitted as trivially small)
python3 src/draw_quotients.py            # Orbit quotient diagrams → quotients/
```

GAP scripts that `Read("symmetric-layout.g")` (e.g. `clayworth_build.g`, `g2_cubic_search.g`) must be run from `src/`. The standalone census scripts in `src/census/` run from anywhere.

## Layout

| Path | Purpose |
|------|---------|
| `src/draw_graph.py` | Layout engine: GAP interface, crossing minimization, TikZ output (imported as a module by `draw_all.py`, `draw_quotients.py`, `graph_codes.py`) |
| `src/symmetric-layout.g` | GAP script: Aut(G), orbit decompositions, quotient graphs, GraphLibrary |
| `src/census/` | Standalone GAP research scripts (PSL(2,q), M11, Macbeath, …) |
| `docs/SYMMETRIC_EMBEDDING.md` | Theory and implementation notes — read this first |
| `docs/RESEARCH_PROGRAM.md` | The PSL(2,q)/SU(2) research line |
| `tex/` | Hand-written papers; compile with pdflatex from inside `tex/` |
| `individual_graphs/` | Mixed generated (`auto-*`) and hand-edited drawings — see its README before regenerating |
| `gallery.pdf` | Rendered gallery (root, showcase artifact) |

## Pipeline

1. GAP computes Aut(G), enumerates subgroups, finds orbit quotients
2. Python selects best quotient, detects cyclic vs cylindrical layout
3. Crossing minimization via path reordering + angular offset search
4. TikZ renders final PDF

Add new graphs in `src/symmetric-layout.g` GraphLibrary; always verify |Aut|.
