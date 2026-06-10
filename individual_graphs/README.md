# Individual Graph Drawings

TikZ drawings of single graphs, one `.tex`/`.pdf` pair each.

- `auto-*.tex` — **generated** by `python3 src/draw_graph.py <name>` (regenerating overwrites them)
- Everything else (e.g. `cage-pair.tex`, `F056A-minimal.tex`, `klein-z2-quotient.tex`) — **hand-edited**; do not regenerate over these

- `*.jpg` — **generated** from the PDFs by `python3 src/make_catalogue.py` (200 dpi, via pdftoppm)

Referenced by `tex/cubic-bi-circulant.tex` via `../individual_graphs/` paths, and by the root `catalogue.pdf` (one page per drawing with descriptions).
