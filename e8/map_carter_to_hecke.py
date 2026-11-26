"""
Map each of the 74 Carter classes to a Hecke algebra and orbital graph.

Strategy:
- Each Carter class with Levi type L uses H(W(E8), W(L))
- We map to the Hecke algebras we've already computed or can easily compute
- The orbital graph is typically the "principal" one (ip=+1 or most regular)
"""

# Mapping: Carter class -> (Hecke algebra descent, orbital graph file, orbital type)
carter_to_hecke = {
    # === Identity and extremal classes ===
    "1": ("E8/E7", "clean_e8_id.pdf", "identity"),  # Trivial class
    "E8": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),  # Coxeter element

    # === E7 subfamily (use H(W(E8), W(E7))) ===
    "E7": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),
    "E7A1": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),

    # === E6 subfamily (use H(W(E7), W(E6)) or H(W(E8), W(E6))) ===
    "E6": ("E7/E6", "clean_e7_mandala.pdf", "ip+1"),
    "E6A2": ("E7/E6", "clean_e6_schlafli.pdf", "ip0"),
    "E6A1A1": ("E7/E6", "clean_e6_schlafli.pdf", "ip0"),

    # === D-series (use D-specific graphs or approximate with E-descent) ===
    "D8": ("E8/E7", "clean_d7_classical.pdf", "proxy"),
    "D7": ("E8/E7", "clean_d7_classical.pdf", "roots"),
    "D7A1": ("E8/E7", "clean_d7_classical.pdf", "roots"),
    "D6": ("E7/E6", "clean_d6_graph.pdf", "roots"),
    "D6A2": ("E7/E6", "clean_d6_graph.pdf", "roots"),
    "D6A1A1": ("E7/E6", "clean_d6_graph.pdf", "roots"),
    "D5": ("E6/D5", "clean_e5_clebsch.pdf", "roots"),
    "D5A3": ("E6/D5", "clean_e5_clebsch.pdf", "roots"),
    "D5A2A1": ("E6/D5", "clean_e5_clebsch.pdf", "roots"),
    "D5A1A1A1": ("E6/D5", "clean_e5_clebsch.pdf", "roots"),
    "D4": ("D5/D4", "clean_d4_graph.pdf", "roots"),
    "D4D4": ("D5/D4", "clean_d4_graph.pdf", "roots"),
    "D4A4": ("D5/D4", "clean_d4_graph.pdf", "roots"),
    "D4A3A1": ("D5/D4", "clean_d4_graph.pdf", "roots"),
    "D4A2A2": ("D5/D4", "clean_d4_graph.pdf", "roots"),
    "D4A2A1A1": ("D5/D4", "clean_d4_graph.pdf", "roots"),
    "D4A1A1A1A1": ("D5/D4", "clean_d4_graph.pdf", "roots"),

    # === F4 subfamily (use F4 root graph) ===
    "F4": ("F4", "clean_f4_roots.pdf", "roots"),
    "F4A4": ("F4", "clean_f4_roots.pdf", "roots"),
    "F4D4": ("F4", "clean_f4_roots.pdf", "roots"),
    "F4A3A1": ("F4", "clean_f4_roots.pdf", "roots"),
    "F4A2A2": ("F4", "clean_f4_roots.pdf", "roots"),
    "F4A2A1A1": ("F4", "clean_f4_roots.pdf", "roots"),
    "F4A1A1A1A1": ("F4", "clean_f4_roots.pdf", "roots"),

    # === G2 subfamily (use G2 root graph) ===
    "G2": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2E6": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2D6": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A6": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A5A1": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A4A2": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A3A3": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A3A2A1": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A2A2A2": ("G2", "clean_g2_roots.pdf", "roots"),
    "G2A2A2A1A1": ("G2", "clean_g2_roots.pdf", "roots"),

    # === Pure A-series (use E8/E7 orbitals by fragmentation level) ===
    "A8": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),
    "A7": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),
    "A7A1": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),
    "A6": ("E8/E7", "clean_e7_mandala.pdf", "ip+1"),
    "A6A2": ("E8/E7", "clean_e6_schlafli.pdf", "ip0"),
    "A6A1A1": ("E8/E7", "clean_e6_schlafli.pdf", "ip0"),
    "A5": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),
    "A5A3": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A5A2A1": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A5A1A1A1": ("E8/E7", "clean_e8_weave.pdf", "weave"),
    "A4": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A4A4": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A4A3A1": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A4A2A2": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A4A2A1A1": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A4A1A1A1A1": ("E8/E7", "clean_e8_weave.pdf", "weave"),
    "A3": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A3A3": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A3A3A2": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A3A3A1A1": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A3A2": ("E8/E7", "clean_a3_roots.pdf", "roots"),
    "A3A2A2A1": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A3A2A1A1A1": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A3A1A1A1A1A1": ("E8/E7", "clean_e8_weave.pdf", "weave"),
    "A2": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A2A2": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A2A2A2": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A2A2A2A2": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A2A2A2A1A1": ("E8/E7", "clean_a2_roots.pdf", "roots"),
    "A2A2A1A1A1A1": ("E8/E7", "clean_e8_weave.pdf", "weave"),
    "A2A1A1A1A1A1A1": ("E8/E7", "clean_e8_weave.pdf", "weave"),
    "A1": ("E8/E7", "clean_e8_weave.pdf", "weave"),
    "A1A1A1A1A1A1A1A1": ("E8/E7", "clean_e8_bloom.pdf", "bloom"),  # The octave
}

# Save the mapping
import json

output = {
    "description": "Mapping from 74 Carter classes to Hecke algebra orbital graphs",
    "total_classes": len(carter_to_hecke),
    "mappings": {
        label: {
            "hecke_algebra": descent,
            "graph_file": graph,
            "orbital_type": orbital
        }
        for label, (descent, graph, orbital) in carter_to_hecke.items()
    }
}

with open("/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"✓ Created mapping for {len(carter_to_hecke)} Carter classes")
print(f"✓ Saved to carter_to_hecke_mapping.json")

# Summary
from collections import Counter
graphs_used = Counter(graph for _, graph, _ in carter_to_hecke.values())
print(f"\nGraphs used:")
for graph, count in sorted(graphs_used.items()):
    print(f"  {graph}: {count} classes")
