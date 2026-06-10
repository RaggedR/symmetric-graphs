# Named Graphs

Symmetry-aware drawings of the graphs formalised in this project.
Layouts computed by `draw_graph.py` using automorphism group quotients.

| File | Graph | V | Sabidussi | g | CSS |
|------|-------|---|-----------|---|-----|
| `heawood-F014A.jpg` | Heawood (F014A) | 14 | Sab(G₄₂, C₃) | 1 | [[21, 2]] |
| `mobius-kantor-F016A.jpg` | Möbius-Kantor (F016A) | 16 | Sab(GL(2,3), C₃) | 2 | [[24, 4]] |
| `klein-F056A.jpg` | Klein (F056A) | 56 | Sab(PSL(2,7), C₃) | 3 | [[84, 6]] |
| `zhou-F182A.jpg` | Zhou (F182A) | 182 | Sab(PSL(2,13), S₃) | — | — |
| `tutte-12cage.jpg` | Tutte 12-cage (F126A) | 126 | from G₂(2) | — | context |
| `langer-63v.jpg` | Langer | 63 (6-reg) | Sab(G₂(2), H₁₉₂) — NOT a quotient of Clayworth (see `clayworth/PR_NOTE.md`) | — | — |

The Tutte 12-cage does not have a {p,3} surface tiling (12F = 2E = 378 gives F = 31.5,
not an integer). It is included as context: it is the incidence graph of GH(2,2) and
the parent of the Clayworth construction.

**Not drawn** (too large): Clayworth (4032 vertices, Sab(G₂(2), C₃), genus 505, [[6048, 1010]])

## TODO: regenerate all images in homogeneous format

The current images use inconsistent styles (different title formats, node sizes, legends).
They should be regenerated using `draw_graph.py` with a consistent template:

- Title: "Name (Foster ID)" e.g. "Klein (F056A)"
- Subtitle: Aut group, |Aut|, vertex count, girth
- Subtitle: vertex stabiliser, edge stabiliser
- Consistent node size and edge style
- No orbit legend (or all have orbit legend)

Source script: `/Users/robin/git/graphviz/draw_graph.py`

## References

- Breuckmann & Terhal, *Constructions and Noise Threshold of Hyperbolic Surface Codes*, arXiv:1506.04029
- Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
