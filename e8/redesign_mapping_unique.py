"""
Redesign the Carter class mapping to ensure 1-1 correspondence:
74 Carter classes → 74 unique orbital graphs

Strategy: For classes using the same Hecke algebra, assign different orbitals.
"""

# Available orbitals for each Hecke algebra
orbitals_available = {
    "E8/E7": [
        ("clean_e8_id.pdf", "ip+2: identity"),
        ("clean_e8_mandala.pdf", "ip+1: mandala (56-reg)"),
        ("clean_e8_bloom.pdf", "ip0: bloom (126-reg)"),
        ("clean_e8_weave.pdf", "ip-1: weave (56-reg)"),
        ("clean_e8_antipodal.pdf", "ip-2: antipodal (K2)"),
    ],
    "E7/E6": [
        ("clean_e7_mandala.pdf", "Gosset graph (27-reg)"),
        # We could compute more E7/E6 orbitals if needed
    ],
    "E6/D5": [
        ("clean_e6_schlafli.pdf", "Schläfli graph SRG(27,10,1,5)"),
        # complement exists too
    ],
    "D-series": [
        ("clean_d7_classical.pdf", "D7 roots"),
        ("clean_d6_graph.pdf", "D6 roots"),
        ("clean_d4_graph.pdf", "D4 roots"),
    ],
    "A-series": [
        ("clean_a3_roots.pdf", "A3 (tetrahedron)"),
        ("clean_a2_roots.pdf", "A2 (triangle)"),
    ],
    "F4": [
        ("clean_f4_roots.pdf", "F4 roots (48 roots)"),
    ],
    "G2": [
        ("clean_g2_roots.pdf", "G2 roots (12 roots, hexagon)"),
    ],
    "E5": [
        ("clean_e5_clebsch.pdf", "D5/E5 Clebsch graph"),
    ],
}

# Carter classes organized by their natural Hecke algebra
# We'll assign unique orbitals to each
classes_by_hecke = {
    "E8/E7": [
        # Use all 5 E8/E7 orbitals
        ("1", 0),           # identity
        ("E8", 1),          # mandala (Coxeter element)
        ("E7", 1),          # mandala
        ("E7A1", 3),        # weave
        ("A8", 1),          # mandala
        ("A7", 2),          # bloom
        ("A7A1", 3),        # weave
        ("A5", 1),          # mandala
        ("A5A1A1A1", 3),    # weave
        ("A4A1A1A1A1", 2),  # bloom
        ("A3A1A1A1A1A1", 3), # weave
        ("A2A2A1A1A1A1", 3), # weave
        ("A2A1A1A1A1A1A1", 3), # weave
        ("A1", 3),          # weave
        ("A1A1A1A1A1A1A1A1", 2), # bloom (the octave)
    ],

    "E7/E6": [
        ("E6", 0),
        ("A6", 0),
    ],

    "E6/D5": [
        ("E6A2", 0),
        ("E6A1A1", 0),
        ("A6A2", 0),
        ("A6A1A1", 0),
    ],

    "E5": [
        ("D5", 0),
        ("D5A3", 0),
        ("D5A2A1", 0),
        ("D5A1A1A1", 0),
    ],

    "D-series": [
        # D8, D7 use D7 graph
        ("D8", 0),
        ("D7", 0),
        ("D7A1", 0),
        # D6 family
        ("D6", 1),
        ("D6A2", 1),
        ("D6A1A1", 1),
        # D4 family
        ("D4", 2),
        ("D4D4", 2),
        ("D4A4", 2),
        ("D4A3A1", 2),
        ("D4A2A2", 2),
        ("D4A2A1A1", 2),
        ("D4A1A1A1A1", 2),
    ],

    "A3": [
        ("A5A3", 0),
        ("A5A2A1", 0),
        ("A4", 0),
        ("A4A4", 0),
        ("A4A3A1", 0),
        ("A3", 0),
        ("A3A3", 0),
        ("A3A3A2", 0),
        ("A3A3A1A1", 0),
        ("A3A2", 0),
    ],

    "A2": [
        ("A4A2A2", 0),
        ("A4A2A1A1", 0),
        ("A3A2A2A1", 0),
        ("A3A2A1A1A1", 0),
        ("A2", 0),
        ("A2A2", 0),
        ("A2A2A2", 0),
        ("A2A2A2A2", 0),
        ("A2A2A2A1A1", 0),
    ],

    "F4": [
        ("F4", 0),
        ("F4A4", 0),
        ("F4D4", 0),
        ("F4A3A1", 0),
        ("F4A2A2", 0),
        ("F4A2A1A1", 0),
        ("F4A1A1A1A1", 0),
    ],

    "G2": [
        ("G2", 0),
        ("G2E6", 0),
        ("G2D6", 0),
        ("G2A6", 0),
        ("G2A5A1", 0),
        ("G2A4A2", 0),
        ("G2A3A3", 0),
        ("G2A3A2A1", 0),
        ("G2A2A2A2", 0),
        ("G2A2A2A1A1", 0),
    ],
}

# Build the 1-1 mapping
carter_to_hecke_unique = {}
used_graphs = set()

for hecke_type, classes in classes_by_hecke.items():
    if hecke_type == "E8/E7":
        orbitals = orbitals_available["E8/E7"]
    elif hecke_type == "E7/E6":
        orbitals = orbitals_available["E7/E6"]
    elif hecke_type == "E6/D5":
        orbitals = orbitals_available["E6/D5"]
    elif hecke_type == "E5":
        orbitals = orbitals_available["E5"]
    elif hecke_type == "D-series":
        orbitals = orbitals_available["D-series"]
    elif hecke_type == "A3":
        orbitals = orbitals_available["A-series"][:1]  # just A3
    elif hecke_type == "A2":
        orbitals = orbitals_available["A-series"][1:]  # just A2
    elif hecke_type == "F4":
        orbitals = orbitals_available["F4"]
    elif hecke_type == "G2":
        orbitals = orbitals_available["G2"]
    else:
        orbitals = []

    for label, orbital_idx in classes:
        graph_file = orbitals[orbital_idx % len(orbitals)][0]
        orbital_type = orbitals[orbital_idx % len(orbitals)][1]

        # Check for uniqueness
        if graph_file in used_graphs and hecke_type == "E8/E7":
            # For E8/E7, we have 5 orbitals, so try to spread them out
            # Find an unused orbital
            for i, (candidate_graph, candidate_type) in enumerate(orbitals):
                if candidate_graph not in used_graphs:
                    graph_file = candidate_graph
                    orbital_type = candidate_type
                    break

        carter_to_hecke_unique[label] = {
            "hecke_algebra": hecke_type,
            "graph_file": graph_file,
            "orbital_type": orbital_type,
        }
        used_graphs.add(graph_file)

# Check for 1-1 mapping
unique_graphs = set(v["graph_file"] for v in carter_to_hecke_unique.values())
print(f"Classes: {len(carter_to_hecke_unique)}")
print(f"Unique graphs: {len(unique_graphs)}")
print(f"1-1 mapping: {len(carter_to_hecke_unique) == len(unique_graphs)}")

if len(carter_to_hecke_unique) != len(unique_graphs):
    from collections import Counter
    graph_counts = Counter(v["graph_file"] for v in carter_to_hecke_unique.values())
    print("\nDuplicates:")
    for graph, count in graph_counts.items():
        if count > 1:
            classes = [k for k, v in carter_to_hecke_unique.items() if v["graph_file"] == graph]
            print(f"  {graph}: {count} classes - {classes}")

# Save
import json
output = {
    "description": "1-1 mapping from 74 Carter classes to unique orbital graphs",
    "total_classes": len(carter_to_hecke_unique),
    "unique_graphs": len(unique_graphs),
    "is_one_to_one": len(carter_to_hecke_unique) == len(unique_graphs),
    "mappings": carter_to_hecke_unique
}

with open("/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\n✓ Saved to carter_to_hecke_mapping.json")
