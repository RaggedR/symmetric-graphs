# GAP Census Scripts

Standalone GAP scripts for constructing and classifying cubic arc-transitive graphs from group data. None depend on `symmetric-layout.g` — each runs directly with `gap -q < script.g`.

- `psl2_*.g` — censuses of cubic graphs on PSL(2,q) coset spaces (C₃, S₃, S₄ point stabilisers); see [docs/RESEARCH_PROGRAM.md](../../docs/RESEARCH_PROGRAM.md) for the motivation
- `m11_cubic_search.g`, `m11_component.g` — searches over M₁₁
- `macbeath_graph.g` — the Macbeath graph (PSL(2,8), Hurwitz)
- `bolza_graphs.g` — graphs related to the Bolza surface (GL(2,3))
- `klein_quartic_codes.g` — Klein quartic code constructions
- `tutte12cage.g` — Tutte 12-cage = incidence graph of GH(2,2); **generated** by `src/build_hexagon.py`
