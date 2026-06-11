# LESSONS

- When a request is ambiguous and one reading would displace approved images
  (`IMAGE_OVERRIDES` — the clayworth/lean canonical set, including the
  hot-linked `lean/named_graphs/`), ASK before rebuilding the catalogue.
  Never let new renders silently outrank `IMAGE_OVERRIDES`.
- Page names hide identities: the `auto-gh22` catalogue page shows
  `langer-63v.jpg`, and `auto-psl2_13` shows `zhou3-F182A.jpg`. "Langer" and
  "zhou" live behind override mappings — check `IMAGE_OVERRIDES` in
  `src/make_catalogue.py` before concluding what an image change will affect.
- `draw_graph.py --style dmt` output is NOT the canonical dmt style and never
  goes into the catalogue. Canonical: `clayworth/generate_named_graphs.py`;
  langer/zhou come from `clayworth/draw_langer.py` / `draw_zhou_quotient.py`.
- This shell is zsh: unquoted `$var` does not word-split. Build CLI flag
  lists with arrays or write flags literally.
