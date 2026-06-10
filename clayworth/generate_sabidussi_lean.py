#!/usr/bin/env python3
"""Generate Lean Sabidussi proof files from GAP output data.

Reads generator arrays and BFS witness words (as Python data),
outputs a complete Lean file proving the graph is a Sabidussi coset graph.
"""

import sys
import os

def chunk_array(arr, chunk_size=40):
    """Split array into lines of chunk_size elements for readability."""
    lines = []
    for i in range(0, len(arr), chunk_size):
        chunk = arr[i:i+chunk_size]
        lines.append(", ".join(str(x) for x in chunk))
    return lines

def generate_lean_file(
    graph_name: str,       # e.g., "dodecahedron"
    lean_graph: str,        # e.g., "dodecahedronGraph"
    lean_graph_import: str, # e.g., "CSSCode.PetersenQuotient"
    n_vertices: int,
    gen1_fwd: list[int],
    gen1_inv: list[int],
    gen2_fwd: list[int],
    gen2_inv: list[int],
    witness_words: list[list[int]],  # witness_words[v] = word from base to v
    group_name: str = "Aut",
    group_description: str = "",
    base: int = 0,
):
    """Generate a complete Lean file."""
    cap = graph_name[0].upper() + graph_name[1:]
    n = n_vertices

    lines = []
    lines.append(f"""/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# {cap}: Sabidussi coset graph

The **{graph_name} graph** ({n} vertices) is a Sabidussi coset graph.
{group_description}
-/

import CSSCode.SabidussiWitness
import {lean_graph_import}

/-! ## Generators -/
""")

    # Generator arrays
    for i, (fwd, inv, label) in enumerate([
        (gen1_fwd, gen1_inv, "gen1"),
        (gen2_fwd, gen2_inv, "gen2"),
    ], 1):
        lines.append(f"private def {graph_name}Gen{i}Fwd : Array (Fin {n}) := #[")
        for line in chunk_array(fwd):
            lines.append(f"  {line},")
        # Remove trailing comma from last line
        lines[-1] = lines[-1].rstrip(',')
        lines.append(f"]")
        lines.append(f"private def {graph_name}Gen{i}Inv : Array (Fin {n}) := #[")
        for line in chunk_array(inv):
            lines.append(f"  {line},")
        lines[-1] = lines[-1].rstrip(',')
        lines.append(f"]")
        lines.append("")

    # Size proofs
    for i in [1, 2]:
        for suf in ["Fwd", "Inv"]:
            lines.append(f"private theorem {graph_name}Gen{i}{suf}_size : {graph_name}Gen{i}{suf}.size = {n} := by native_decide")
    lines.append("")

    # Equiv.Perm definitions
    for i in [1, 2]:
        lines.append(f"def {graph_name}Gen{i} : Equiv.Perm (Fin {n}) where")
        lines.append(f"  toFun i := {graph_name}Gen{i}Fwd[i.val]'(by have := {graph_name}Gen{i}Fwd_size; omega)")
        lines.append(f"  invFun i := {graph_name}Gen{i}Inv[i.val]'(by have := {graph_name}Gen{i}Inv_size; omega)")
        lines.append(f"  left_inv := by native_decide")
        lines.append(f"  right_inv := by native_decide")
        lines.append("")

    # Gens function
    lines.append(f"def {graph_name}Gens : Fin 2 → Equiv.Perm (Fin {n})")
    lines.append(f"  | 0 => {graph_name}Gen1")
    lines.append(f"  | 1 => {graph_name}Gen2")
    lines.append("")

    # Adjacency preservation
    for i in [1, 2]:
        lines.append(f"theorem {graph_name}Gen{i}_adj : ∀ u v : Fin {n},")
        lines.append(f"    {lean_graph}.Adj u v → {lean_graph}.Adj ({graph_name}Gen{i} u) ({graph_name}Gen{i} v) := by")
        lines.append(f"  native_decide")
        lines.append("")
        lines.append(f"theorem {graph_name}Gen{i}_inv_adj : ∀ u v : Fin {n},")
        lines.append(f"    {lean_graph}.Adj u v → {lean_graph}.Adj ({graph_name}Gen{i}.symm u) ({graph_name}Gen{i}.symm v) := by")
        lines.append(f"  native_decide")
        lines.append("")

    # Witness words
    lines.append(f"/-! ## BFS witness words -/")
    lines.append("")
    lines.append(f"private def {graph_name}WitnessWords : Fin {n} → List (Fin 4) := fun v =>")
    lines.append(f"  match v with")
    for v in range(n):
        w = witness_words[v]
        if len(w) == 0:
            lines.append(f"  | {v} => []")
        else:
            wstr = ", ".join(str(x) for x in w)
            lines.append(f"  | {v} => [{wstr}]")
    lines.append("")

    lines.append(f"theorem {graph_name}Witness_correct :")
    lines.append(f"    ∀ v : Fin {n}, applyWordGen {graph_name}Gens ({graph_name}WitnessWords v) {base} = v := by")
    lines.append(f"  native_decide")
    lines.append("")

    # SabidussiWitness
    lines.append(f"/-! ## Sabidussi witness -/")
    lines.append("")
    lines.append(f"def {graph_name}SabidussiWitness : SabidussiWitness {lean_graph} where")
    lines.append(f"  numGens := 2")
    lines.append(f"  gens := {graph_name}Gens")
    lines.append(f"  gens_preserve_adj := by")
    lines.append(f"    intro i; match i with")
    lines.append(f"    | ⟨0, _⟩ => exact {graph_name}Gen1_adj")
    lines.append(f"    | ⟨1, _⟩ => exact {graph_name}Gen2_adj")
    lines.append(f"  gens_inv_preserve_adj := by")
    lines.append(f"    intro i; match i with")
    lines.append(f"    | ⟨0, _⟩ => exact {graph_name}Gen1_inv_adj")
    lines.append(f"    | ⟨1, _⟩ => exact {graph_name}Gen2_inv_adj")
    lines.append(f"  base := {base}")
    lines.append(f"  witnessWords := {graph_name}WitnessWords")
    lines.append(f"  witness_correct := {graph_name}Witness_correct")
    lines.append("")

    # Main theorem
    lines.append(f"noncomputable def {graph_name}SabidussiIso := {graph_name}SabidussiWitness.toSabidussiIso")

    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    print("This script is imported by other scripts to generate Lean files.")
    print("Use generate_lean_file() with the appropriate parameters.")
