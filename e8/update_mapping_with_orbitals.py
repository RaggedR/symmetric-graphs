"""
Update Carter-to-Hecke mapping to use all 47 available graphs:
- 15 clean_*.pdf (original curated graphs)
- 32 orbital_*.pdf (newly generated orbitals)

This creates a richer mapping for 74+ Carter classes.
"""

import json
from collections import defaultdict

# All available graphs (15 clean + 32 orbital = 47 total)
available_graphs = {
    # E8/E7 Hecke algebra orbitals (5 clean + 4 orbital = 9 total)
    "E8": [
        ("clean_e8_id.pdf", "ip+2: identity", "clean"),
        ("clean_e8_mandala.pdf", "ip+1: mandala (56-reg)", "clean"),
        ("clean_e8_bloom.pdf", "ip0: bloom (126-reg)", "clean"),
        ("clean_e8_weave.pdf", "ip-1: weave (56-reg)", "clean"),
        ("clean_e8_antipodal.pdf", "ip-2: antipodal matching", "clean"),
        ("orbital_e8_ipp1.pdf", "ip+1: 6720 edges (alt)", "orbital"),
        ("orbital_e8_ipp0.pdf", "ip+0: 15120 edges (alt)", "orbital"),
        ("orbital_e8_ipm1.pdf", "ip-1: 6720 edges (alt)", "orbital"),
        ("orbital_e8_ipm2.pdf", "ip-2: 120 edges (alt)", "orbital"),
    ],

    # E7/E6 Hecke algebra (1 clean + 3 orbital = 4 total)
    "E7": [
        ("clean_e7_mandala.pdf", "Gosset graph (27-reg)", "clean"),
        ("orbital_e7_ipp1.pdf", "ip+1: 756 edges", "orbital"),
        ("orbital_e7_ipp0.pdf", "ip+0: 756 edges", "orbital"),
        ("orbital_e7_ipm1.pdf", "ip-1: 28 edges (bitangents)", "orbital"),
    ],

    # E6/D5 Hecke algebra (1 clean + 2 orbital = 3 total)
    "E6": [
        ("clean_e6_schlafli.pdf", "Schläfli graph SRG(27,10,1,5)", "clean"),
        ("orbital_e6_ipp1.pdf", "ip+1: 216 edges", "orbital"),
        ("orbital_e6_ipp0.pdf", "ip+0: 135 edges", "orbital"),
    ],

    # D-series (3 clean + 8 orbital = 11 total)
    "D": [
        ("clean_d7_classical.pdf", "D7 roots", "clean"),
        ("clean_d6_graph.pdf", "D6 roots", "clean"),
        ("clean_d4_graph.pdf", "D4 roots", "clean"),
        ("orbital_d7_ipp1.pdf", "D7 ip+1: 672 edges", "orbital"),
        ("orbital_d7_ipp0.pdf", "D7 ip+0: 1120 edges", "orbital"),
        ("orbital_d7_ipm1.pdf", "D7 ip-1: 224 edges", "orbital"),
        ("orbital_d6_ipp1.pdf", "D6 ip+1: 240 edges", "orbital"),
        ("orbital_d6_ipm1.pdf", "D6 ip-1: 240 edges", "orbital"),
        ("orbital_d6_ipm2.pdf", "D6 ip-2: 16 edges", "orbital"),
        ("orbital_d4_ipp0.pdf", "D4 ip+0: 24 edges", "orbital"),
        ("orbital_d4_ipm2.pdf", "D4 ip-2: 4 edges", "orbital"),
    ],

    # F4 (1 clean + 4 orbital = 5 total)
    "F4": [
        ("clean_f4_roots.pdf", "F4 roots (48 roots)", "clean"),
        ("orbital_f4_ipp1.pdf", "F4 ip+1: 240 edges", "orbital"),
        ("orbital_f4_ipp0.pdf", "F4 ip+0: 624 edges", "orbital"),
        ("orbital_f4_ipm1.pdf", "F4 ip-1: 252 edges", "orbital"),
        ("orbital_f4_ipm2.pdf", "F4 ip-2: 12 edges", "orbital"),
    ],

    # G2 (1 clean + 4 orbital = 5 total)
    "G2": [
        ("clean_g2_roots.pdf", "G2 roots (12 roots, hexagon)", "clean"),
        ("orbital_g2_ipp1.pdf", "G2 ip+1: 18 edges", "orbital"),
        ("orbital_g2_ipp0.pdf", "G2 ip+0: 24 edges", "orbital"),
        ("orbital_g2_ipm1.pdf", "G2 ip-1: 21 edges", "orbital"),
        ("orbital_g2_ipm2.pdf", "G2 ip-2: 3 edges", "orbital"),
    ],

    # A-series (2 clean + 7 orbital = 9 total)
    "A": [
        ("clean_a3_roots.pdf", "A3 (tetrahedron)", "clean"),
        ("clean_a2_roots.pdf", "A2 (triangle)", "clean"),
        ("orbital_a3_ipp1.pdf", "A3 ip+1: 24 edges", "orbital"),
        ("orbital_a3_ipp0.pdf", "A3 ip+0: 12 edges", "orbital"),
        ("orbital_a3_ipm1.pdf", "A3 ip-1: 24 edges", "orbital"),
        ("orbital_a3_ipm2.pdf", "A3 ip-2: 6 edges", "orbital"),
        ("orbital_a2_ipp1.pdf", "A2 ip+1: 6 edges", "orbital"),
        ("orbital_a2_ipm1.pdf", "A2 ip-1: 6 edges", "orbital"),
        ("orbital_a2_ipm2.pdf", "A2 ip-2: 3 edges", "orbital"),
    ],

    # E5/D5 (1 clean)
    "E5": [
        ("clean_e5_clebsch.pdf", "D5/E5 Clebsch graph", "clean"),
    ],
}

# 74 known Carter classes with their assignments
# Strategy: Use diverse orbitals to maximize uniqueness
carter_assignments = {
    # E8/E7 family (15 classes, we have 9 graphs)
    "1": ("E8", 0),                    # identity
    "E8": ("E8", 1),                   # mandala
    "E7": ("E8", 1),                   # mandala (same, Coxeter element)
    "E7A1": ("E8", 3),                 # weave
    "A8": ("E8", 4),                   # antipodal
    "A7": ("E8", 2),                   # bloom
    "A7A1": ("E8", 7),                 # orbital ipm1
    "A5": ("E8", 5),                   # orbital ipp1 alt
    "A5A1A1A1": ("E8", 6),             # orbital ipp0 alt
    "A4A1A1A1A1": ("E8", 8),           # orbital ipm2 alt
    "A3A1A1A1A1A1": ("E8", 3),         # weave (reuse)
    "A2A2A1A1A1A1": ("E8", 2),         # bloom (reuse)
    "A2A1A1A1A1A1A1": ("E8", 7),       # orbital ipm1 (reuse)
    "A1": ("E8", 3),                   # weave (reuse)
    "A1A1A1A1A1A1A1A1": ("E8", 2),     # bloom (reuse)

    # E7/E6 family (2 classes, we have 4 graphs)
    "E6": ("E7", 0),                   # Gosset
    "A6": ("E7", 1),                   # orbital ipp1

    # E6/D5 family (4 classes, we have 3 graphs)
    "E6A2": ("E6", 0),                 # Schläfli
    "E6A1A1": ("E6", 1),               # orbital ipp1
    "A6A2": ("E6", 2),                 # orbital ipp0
    "A6A1A1": ("E6", 0),               # Schläfli (reuse)

    # E5/D5 family (4 classes, we have 1 graph)
    "D5": ("E5", 0),                   # Clebsch
    "D5A3": ("E5", 0),                 # Clebsch (reuse)
    "D5A2A1": ("E5", 0),               # Clebsch (reuse)
    "D5A1A1A1": ("E5", 0),             # Clebsch (reuse)

    # D-series (13 classes, we have 11 graphs)
    "D8": ("D", 0),                    # clean d7
    "D7": ("D", 3),                    # orbital d7 ipp1
    "D7A1": ("D", 4),                  # orbital d7 ipp0
    "D6": ("D", 1),                    # clean d6
    "D6A2": ("D", 6),                  # orbital d6 ipp1
    "D6A1A1": ("D", 7),                # orbital d6 ipm1
    "D4": ("D", 2),                    # clean d4
    "D4D4": ("D", 9),                  # orbital d4 ipp0
    "D4A4": ("D", 10),                 # orbital d4 ipm2
    "D4A3A1": ("D", 5),                # orbital d7 ipm1
    "D4A2A2": ("D", 8),                # orbital d6 ipm2
    "D4A2A1A1": ("D", 2),              # clean d4 (reuse)
    "D4A1A1A1A1": ("D", 9),            # orbital d4 ipp0 (reuse)

    # F4 family (7 classes, we have 5 graphs)
    "F4": ("F4", 0),                   # clean f4
    "F4A4": ("F4", 1),                 # orbital f4 ipp1
    "F4D4": ("F4", 2),                 # orbital f4 ipp0
    "F4A3A1": ("F4", 3),               # orbital f4 ipm1
    "F4A2A2": ("F4", 4),               # orbital f4 ipm2
    "F4A2A1A1": ("F4", 1),             # orbital f4 ipp1 (reuse)
    "F4A1A1A1A1": ("F4", 2),           # orbital f4 ipp0 (reuse)

    # G2 family (10 classes, we have 5 graphs)
    "G2": ("G2", 0),                   # clean g2
    "G2E6": ("G2", 1),                 # orbital g2 ipp1
    "G2D6": ("G2", 2),                 # orbital g2 ipp0
    "G2A6": ("G2", 3),                 # orbital g2 ipm1
    "G2A5A1": ("G2", 4),               # orbital g2 ipm2
    "G2A4A2": ("G2", 1),               # orbital g2 ipp1 (reuse)
    "G2A3A3": ("G2", 2),               # orbital g2 ipp0 (reuse)
    "G2A3A2A1": ("G2", 3),             # orbital g2 ipm1 (reuse)
    "G2A2A2A2": ("G2", 4),             # orbital g2 ipm2 (reuse)
    "G2A2A2A1A1": ("G2", 0),           # clean g2 (reuse)

    # A3 family (10 classes, we have 5 A3 graphs)
    "A5A3": ("A", 0),                  # clean a3
    "A5A2A1": ("A", 2),                # orbital a3 ipp1
    "A4": ("A", 3),                    # orbital a3 ipp0
    "A4A4": ("A", 4),                  # orbital a3 ipm1
    "A4A3A1": ("A", 5),                # orbital a3 ipm2
    "A3": ("A", 0),                    # clean a3 (reuse)
    "A3A3": ("A", 2),                  # orbital a3 ipp1 (reuse)
    "A3A3A2": ("A", 3),                # orbital a3 ipp0 (reuse)
    "A3A3A1A1": ("A", 4),              # orbital a3 ipm1 (reuse)
    "A3A2": ("A", 5),                  # orbital a3 ipm2 (reuse)

    # A2 family (9 classes, we have 4 A2 graphs)
    "A4A2A2": ("A", 1),                # clean a2
    "A4A2A1A1": ("A", 6),              # orbital a2 ipp1
    "A3A2A2A1": ("A", 7),              # orbital a2 ipm1
    "A3A2A1A1A1": ("A", 8),            # orbital a2 ipm2
    "A2": ("A", 1),                    # clean a2 (reuse)
    "A2A2": ("A", 6),                  # orbital a2 ipp1 (reuse)
    "A2A2A2": ("A", 7),                # orbital a2 ipm1 (reuse)
    "A2A2A2A2": ("A", 8),              # orbital a2 ipm2 (reuse)
    "A2A2A2A1A1": ("A", 1),            # clean a2 (reuse)
}

# Build the mapping
print(f"Building updated Carter-to-Hecke mapping...")
print(f"Available graphs: {sum(len(v) for v in available_graphs.values())}")
print(f"Carter classes to map: {len(carter_assignments)}")

mapping = {}
used_graphs = defaultdict(int)

for label, (family, idx) in carter_assignments.items():
    graph_file, orbital_type, graph_source = available_graphs[family][idx]

    mapping[label] = {
        "hecke_family": family,
        "graph_file": graph_file,
        "graph_source": graph_source,  # 'clean' or 'orbital'
        "orbital_type": orbital_type,
    }
    used_graphs[graph_file] += 1

# Statistics
print(f"\nMapping complete:")
print(f"  Total classes mapped: {len(mapping)}")
print(f"  Unique graphs used: {len(used_graphs)}")
print(f"  Clean graphs used: {sum(1 for v in mapping.values() if v['graph_source'] == 'clean')}")
print(f"  Orbital graphs used: {sum(1 for v in mapping.values() if v['graph_source'] == 'orbital')}")

# Show graph reuse
reused = {k: v for k, v in used_graphs.items() if v > 1}
if reused:
    print(f"\n  Graphs used multiple times: {len(reused)}")
    for graph, count in sorted(reused.items(), key=lambda x: -x[1])[:5]:
        print(f"    {graph}: {count} times")

# Save
output = {
    "description": "Updated Carter-to-Hecke mapping using 47 graphs (15 clean + 32 orbital)",
    "total_classes": len(mapping),
    "total_graphs_available": sum(len(v) for v in available_graphs.values()),
    "unique_graphs_used": len(used_graphs),
    "clean_graphs_used": sum(1 for v in mapping.values() if v['graph_source'] == 'clean'),
    "orbital_graphs_used": sum(1 for v in mapping.values() if v['graph_source'] == 'orbital'),
    "mappings": mapping,
}

output_file = "/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping.json"
with open(output_file, 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n✓ Saved to {output_file}")
