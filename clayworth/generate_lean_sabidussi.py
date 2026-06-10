#!/usr/bin/env python3
"""Generate Lean Sabidussi proof files from GAP data files."""

import re, os, sys

def parse_gap_data(filename):
    """Parse a GAP data file into a dict of arrays."""
    data = {}
    with open(filename) as f:
        content = f.read()
    # Parse "prefixFieldName := [ ... ];"
    for m in re.finditer(r'(\w+)\s*:=\s*(\[.*?\]|(\d+))\s*;', content, re.DOTALL):
        key = m.group(1)
        val_str = m.group(2).strip()
        if val_str.startswith('['):
            # Parse nested list
            val_str = val_str.replace('\n', ' ').replace('  ', ' ')
            try:
                val = eval(val_str)  # Safe for numeric lists
            except:
                val = val_str
        else:
            val = int(val_str)
        data[key] = val
    return data

def chunk_array(arr, per_line=30):
    lines = []
    for i in range(0, len(arr), per_line):
        chunk = arr[i:i+per_line]
        lines.append("  " + ", ".join(str(x) for x in chunk))
    return ",\n".join(lines)

def gen_lean(name, lean_graph_name, import_line, n, gen1_fwd, gen1_inv,
             gen2_fwd, gen2_inv, words, edges, group_desc, base=0,
             need_graph_def=False):
    """Generate complete Lean file content."""
    cap = name[0].upper() + name[1:]
    # Edge count for degree
    degree = None
    if edges:
        from collections import Counter
        deg_count = Counter()
        for u, v in edges:
            deg_count[u] += 1
            deg_count[v] += 1
        degree = deg_count[0] if 0 in deg_count else deg_count[min(deg_count)]

    lines = [f"""/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# {cap}: Sabidussi coset graph

{group_desc}
-/

import CSSCode.SabidussiWitness
import {import_line}
"""]

    # Graph definition if needed
    if need_graph_def:
        lines.append(f"/-! ## Graph definition -/\n")
        lines.append(f"private def {name}Edges : Array (Fin {n} × Fin {n}) := #[")
        edge_strs = [f"  ({u}, {v})" for u, v in edges]
        lines.append(",\n".join(edge_strs))
        lines.append(f"]\n")
        lines.append(f"""private def {name}AdjBool (u v : Fin {n}) : Bool :=
  {name}Edges.any fun (a, b) => (u == a && v == b) || (u == b && v == a)

def {lean_graph_name} : SimpleGraph (Fin {n}) where
  Adj u v := {name}AdjBool u v
  symm u v := by simp only [{name}AdjBool]; revert u v; native_decide
  loopless := ⟨fun u => by simp only [{name}AdjBool]; revert u; native_decide⟩

instance {name}DecAdj : DecidableRel {lean_graph_name}.Adj :=
  fun u v => inferInstanceAs (Decidable ({name}AdjBool u v))
""")

    # Generator arrays
    lines.append(f"/-! ## Generators -/\n")
    for i, (fwd, inv) in enumerate([(gen1_fwd, gen1_inv), (gen2_fwd, gen2_inv)], 1):
        lines.append(f"private def {name}Gen{i}Fwd : Array (Fin {n}) := #[\n{chunk_array(fwd)}\n]")
        lines.append(f"private def {name}Gen{i}Inv : Array (Fin {n}) := #[\n{chunk_array(inv)}\n]\n")

    # Size proofs + Equiv.Perm defs
    for i in [1, 2]:
        for suf in ["Fwd", "Inv"]:
            lines.append(f"private theorem {name}Gen{i}{suf}_size : {name}Gen{i}{suf}.size = {n} := by native_decide")
    lines.append("")

    for i in [1, 2]:
        lines.append(f"""def {name}Gen{i} : Equiv.Perm (Fin {n}) where
  toFun i := {name}Gen{i}Fwd[i.val]'(by have := {name}Gen{i}Fwd_size; omega)
  invFun i := {name}Gen{i}Inv[i.val]'(by have := {name}Gen{i}Inv_size; omega)
  left_inv := by native_decide
  right_inv := by native_decide
""")

    # Gens function
    lines.append(f"""def {name}Gens : Fin 2 → Equiv.Perm (Fin {n})
  | 0 => {name}Gen1
  | 1 => {name}Gen2
""")

    # Adjacency preservation
    for i in [1, 2]:
        lines.append(f"""theorem {name}Gen{i}_adj : ∀ u v : Fin {n},
    {lean_graph_name}.Adj u v → {lean_graph_name}.Adj ({name}Gen{i} u) ({name}Gen{i} v) := by
  native_decide

theorem {name}Gen{i}_inv_adj : ∀ u v : Fin {n},
    {lean_graph_name}.Adj u v → {lean_graph_name}.Adj ({name}Gen{i}.symm u) ({name}Gen{i}.symm v) := by
  native_decide
""")

    # Witness words
    lines.append(f"/-! ## BFS witness words -/\n")
    lines.append(f"private def {name}WitnessWords : Fin {n} → List (Fin 4) := fun v =>")
    lines.append(f"  match v with")
    for v in range(n):
        w = words[v]
        if w is None or (isinstance(w, str) and w == 'fail'):
            lines.append(f"  | {v} => []  -- unreachable")
        elif len(w) == 0:
            lines.append(f"  | {v} => []")
        else:
            wstr = ", ".join(str(x) for x in w)
            lines.append(f"  | {v} => [{wstr}]")
    # For n > 10, Lean may not recognize the match as exhaustive
    if n > 10:
        lines.append(f"  | _ => []  -- unreachable, needed for exhaustiveness")
    lines.append("")

    lines.append(f"""theorem {name}Witness_correct :
    ∀ v : Fin {n}, applyWordGen {name}Gens ({name}WitnessWords v) {base} = v := by
  native_decide
""")

    # SabidussiWitness
    lines.append(f"""/-! ## Sabidussi witness -/

def {name}SabidussiWitness : SabidussiWitness {lean_graph_name} where
  numGens := 2
  gens := {name}Gens
  gens_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact {name}Gen1_adj
    | ⟨1, _⟩ => exact {name}Gen2_adj
  gens_inv_preserve_adj := by
    intro i; match i with
    | ⟨0, _⟩ => exact {name}Gen1_inv_adj
    | ⟨1, _⟩ => exact {name}Gen2_inv_adj
  base := {base}
  witnessWords := {name}WitnessWords
  witness_correct := {name}Witness_correct

noncomputable def {name}SabidussiIso := {name}SabidussiWitness.toSabidussiIso
""")

    return "\n".join(lines)

# ============================================================
# Generate all files
# ============================================================

GRAPHS = [
    {
        "name": "dodec",
        "file": "/tmp/dodec_data.g",
        "prefix": "dodec",
        "lean_graph": "dodecGraph",
        "import": "CSSCode.DodecSabidussiGraph",
        "need_graph": True,
        "desc": "The **dodecahedron** (20 vertices, cubic) is Sab(A₅, C₃).\n|A₅| = 60, |C₃| = 3, 60/3 = 20 vertices.",
    },
    {
        "name": "petersen",
        "file": "/tmp/petersen_data.g",
        "prefix": "pet",
        "lean_graph": "petGraph",
        "import": "CSSCode.PetersenSabidussiGraph",
        "need_graph": True,
        "desc": "The **Petersen graph** (10 vertices, cubic) is Sab(A₅, S₃).\n|A₅| = 60, |S₃| = 6, 60/6 = 10 vertices.",
    },
    {
        "name": "heawood",
        "file": "/tmp/heawood_data.g",
        "prefix": "heaw",
        "lean_graph": "heawGraph",
        "import": "CSSCode.HeawoodSabidussiGraph",
        "need_graph": True,
        "desc": "The **Heawood graph** (14 vertices, cubic) is Sab(PGL(2,7), H₂₄).\n|PGL(2,7)| = 336, |H₂₄| = 24, 336/24 = 14 vertices.",
    },
    {
        "name": "mk",
        "file": "/tmp/mk_data.g",
        "prefix": "mk",
        "lean_graph": "mkSabGraph",
        "import": "CSSCode.MKSabidussiGraph",
        "need_graph": True,
        "desc": "The **Möbius-Kantor graph** (16 vertices, cubic) is Sab(GL(2,3), C₃).\n|GL(2,3)| = 48, |C₃| = 3, 48/3 = 16 vertices.",
    },
    {
        "name": "klein",
        "file": "/tmp/klein_data.g",
        "prefix": "klein",
        "lean_graph": "kleinSabGraph",
        "import": "CSSCode.KleinSabidussiGraph",
        "need_graph": True,
        "desc": "The **Klein quartic graph** (56 vertices, cubic) is Sab(PSL(2,7), C₃).\n|PSL(2,7)| = 168, |C₃| = 3, 168/3 = 56 vertices.",
    },
    {
        "name": "zhou3",
        "file": "/tmp/zhou3_data.g",
        "prefix": "zhou3",
        "lean_graph": "zhou3Graph",
        "import": "CSSCode.Zhou3SabidussiGraph",
        "need_graph": True,
        "desc": "The **Zhou-3 graph** (F182A, 182 vertices, cubic) is Sab(PSL(2,13), S₃).\n|PSL(2,13)| = 1092, |S₃| = 6, 1092/6 = 182 vertices.",
    },
    {
        "name": "zhou6",
        "file": "/tmp/zhou6_data.g",
        "prefix": "zhou6",
        "lean_graph": "zhou6Graph",
        "import": "CSSCode.Zhou6SabidussiGraph",
        "need_graph": True,
        "desc": "The **Zhou-6 graph** (91 vertices, 6-regular) is Sab(PSL(2,13), D₁₂).\n|PSL(2,13)| = 1092, |D₁₂| = 12, 1092/12 = 91 vertices.",
    },
]

outdir = os.path.join(os.path.dirname(__file__), "CSSCode")

for g in GRAPHS:
    print(f"Processing {g['name']}...")
    data = parse_gap_data(g["file"])
    prefix = g["prefix"]
    n = data[f"{prefix}N"]
    gen1_fwd = data[f"{prefix}Gen1Fwd"]
    gen1_inv = data[f"{prefix}Gen1Inv"]
    gen2_fwd = data[f"{prefix}Gen2Fwd"]
    gen2_inv = data[f"{prefix}Gen2Inv"]
    words = data[f"{prefix}Words"]
    edges = data.get(f"{prefix}Edges", [])

    lean_code = gen_lean(
        name=g["name"],
        lean_graph_name=g["lean_graph"],
        import_line=g["import"],
        n=n,
        gen1_fwd=gen1_fwd,
        gen1_inv=gen1_inv,
        gen2_fwd=gen2_fwd,
        gen2_inv=gen2_inv,
        words=words,
        edges=edges,
        group_desc=g["desc"],
        need_graph_def=g.get("need_graph", False),
    )

    # For graphs that define their own graph, the import should be self
    # Actually, we put everything in one file: graph def + sabidussi proof
    # So the import just needs SabidussiWitness
    lean_code = lean_code.replace(f"import {g['import']}\n", "")

    outfile = os.path.join(outdir, f"{g['name'].title().replace('_','')}Sabidussi.lean")
    with open(outfile, 'w') as f:
        f.write(lean_code)
    print(f"  -> {outfile} ({n} vertices, {len(edges)} edges)")

print("Done!")
