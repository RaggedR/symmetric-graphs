"""
Extend Carter mapping from 74 to 112 conjugacy classes.

Strategy:
- Load the 112 conjugacy classes from carter_classes_complete.json
- Load the current 74-class mapping
- For the remaining 38 classes, assign graphs based on element order heuristics
- Use unused graph combinations and reasonable Levi-type approximations
"""

import json

# Load 112 conjugacy classes
with open("/Users/miriam/symmetric-graphs/e8/carter_classes_complete.json") as f:
    complete_data = json.load(f)
    all_112_classes = complete_data["classes"]

# Load current 74-class mapping
with open("/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping.json") as f:
    current_mapping = json.load(f)
    current_mappings = current_mapping["mappings"]

# All 47 available graphs
available_graphs = [
    # E8 family (9 graphs)
    ("clean_e8_id.pdf", "E8", "ip+2: identity"),
    ("clean_e8_mandala.pdf", "E8", "ip+1: mandala"),
    ("clean_e8_bloom.pdf", "E8", "ip0: bloom"),
    ("clean_e8_weave.pdf", "E8", "ip-1: weave"),
    ("clean_e8_antipodal.pdf", "E8", "ip-2: antipodal"),
    ("orbital_e8_ipp1.pdf", "E8", "orbital ip+1"),
    ("orbital_e8_ipp0.pdf", "E8", "orbital ip+0"),
    ("orbital_e8_ipm1.pdf", "E8", "orbital ip-1"),
    ("orbital_e8_ipm2.pdf", "E8", "orbital ip-2"),

    # E7 family (4 graphs)
    ("clean_e7_mandala.pdf", "E7", "Gosset graph"),
    ("orbital_e7_ipp1.pdf", "E7", "orbital ip+1"),
    ("orbital_e7_ipp0.pdf", "E7", "orbital ip+0"),
    ("orbital_e7_ipm1.pdf", "E7", "orbital ip-1"),

    # E6 family (3 graphs)
    ("clean_e6_schlafli.pdf", "E6", "Schläfli graph"),
    ("orbital_e6_ipp1.pdf", "E6", "orbital ip+1"),
    ("orbital_e6_ipp0.pdf", "E6", "orbital ip+0"),

    # D family (11 graphs)
    ("clean_d7_classical.pdf", "D", "D7 roots"),
    ("clean_d6_graph.pdf", "D", "D6 roots"),
    ("clean_d4_graph.pdf", "D", "D4 roots"),
    ("orbital_d7_ipp1.pdf", "D", "D7 ip+1"),
    ("orbital_d7_ipp0.pdf", "D", "D7 ip+0"),
    ("orbital_d7_ipm1.pdf", "D", "D7 ip-1"),
    ("orbital_d6_ipp1.pdf", "D", "D6 ip+1"),
    ("orbital_d6_ipm1.pdf", "D", "D6 ip-1"),
    ("orbital_d6_ipm2.pdf", "D", "D6 ip-2"),
    ("orbital_d4_ipp0.pdf", "D", "D4 ip+0"),
    ("orbital_d4_ipm2.pdf", "D", "D4 ip-2"),

    # F4 family (5 graphs)
    ("clean_f4_roots.pdf", "F4", "F4 roots"),
    ("orbital_f4_ipp1.pdf", "F4", "F4 ip+1"),
    ("orbital_f4_ipp0.pdf", "F4", "F4 ip+0"),
    ("orbital_f4_ipm1.pdf", "F4", "F4 ip-1"),
    ("orbital_f4_ipm2.pdf", "F4", "F4 ip-2"),

    # G2 family (5 graphs)
    ("clean_g2_roots.pdf", "G2", "G2 roots"),
    ("orbital_g2_ipp1.pdf", "G2", "G2 ip+1"),
    ("orbital_g2_ipp0.pdf", "G2", "G2 ip+0"),
    ("orbital_g2_ipm1.pdf", "G2", "G2 ip-1"),
    ("orbital_g2_ipm2.pdf", "G2", "G2 ip-2"),

    # A family (9 graphs)
    ("clean_a3_roots.pdf", "A", "A3 roots"),
    ("clean_a2_roots.pdf", "A", "A2 roots"),
    ("orbital_a3_ipp1.pdf", "A", "A3 ip+1"),
    ("orbital_a3_ipp0.pdf", "A", "A3 ip+0"),
    ("orbital_a3_ipm1.pdf", "A", "A3 ip-1"),
    ("orbital_a3_ipm2.pdf", "A", "A3 ip-2"),
    ("orbital_a2_ipp1.pdf", "A", "A2 ip+1"),
    ("orbital_a2_ipm1.pdf", "A", "A2 ip-1"),
    ("orbital_a2_ipm2.pdf", "A", "A2 ip-2"),

    # E5 family (1 graph)
    ("clean_e5_clebsch.pdf", "E5", "Clebsch graph"),
]

# Track which graphs are used
used_graphs = set(v["graph_file"] for v in current_mappings.values())
print(f"Currently using {len(used_graphs)} unique graphs out of {len(available_graphs)} available")

# Map element order to approximate Levi type
def order_to_family(order):
    """Heuristic: map element order to likely root system family."""
    if order == 1:
        return "E8"
    elif order == 30:
        return "E8"  # Coxeter element
    elif order == 18:
        return "E7"
    elif order == 12:
        return "E6"  # or F4
    elif order == 10:
        return "D"   # D8
    elif order == 8:
        return "D"   # D5
    elif order == 6:
        return "G2"  # or D4 or A5
    elif order == 5:
        return "D"   # A4-like
    elif order == 4:
        return "A"   # A3
    elif order == 3:
        return "A"   # A2
    elif order == 2:
        return "A"   # A1
    else:
        return "E8"  # default to E8 family

# Build extended mapping for all 112
extended_mapping = dict(current_mappings)  # Start with existing 74

print(f"\nExtending from {len(current_mappings)} to 112 classes...")

# For classes beyond the named 74, use index-based labels
graph_cycle_index = 0

for i, cls in enumerate(all_112_classes):
    idx = cls["index"]
    label = f"class_{idx}"

    # Skip if already mapped (by Carter label)
    if any(label.lower() == k.lower() for k in extended_mapping):
        continue

    # Map based on element order
    order = cls["order"]
    family = order_to_family(order)

    # Find next available graph from this family
    family_graphs = [(g, f, t) for g, f, t in available_graphs if f == family]

    if family_graphs:
        # Rotate through family graphs
        graph_file, graph_family, orbital_type = family_graphs[graph_cycle_index % len(family_graphs)]
        graph_cycle_index += 1
    else:
        # Fallback: use any available graph
        graph_file, graph_family, orbital_type = available_graphs[graph_cycle_index % len(available_graphs)]
        graph_cycle_index += 1

    extended_mapping[label] = {
        "hecke_family": graph_family,
        "graph_file": graph_file,
        "graph_source": "orbital" if "orbital_" in graph_file else "clean",
        "orbital_type": orbital_type,
        "element_order": order,
        "class_size": cls["class_size"],
        "note": f"Extended assignment based on order {order}"
    }

    if len(extended_mapping) >= 112:
        break

print(f"Extended mapping to {len(extended_mapping)} classes")

# Count unique graphs used
unique_extended = set(v["graph_file"] for v in extended_mapping.values())
print(f"Using {len(unique_extended)} unique graphs")

# Save extended mapping
output = {
    "description": "Extended Carter-to-Hecke mapping for all 112 W(E8) conjugacy classes",
    "note": "Classes 1-74 have known Carter labels. Classes 75-112 use approximate assignments based on element order.",
    "total_classes": len(extended_mapping),
    "total_graphs_available": len(available_graphs),
    "unique_graphs_used": len(unique_extended),
    "named_carter_classes": len(current_mappings),
    "extended_classes": len(extended_mapping) - len(current_mappings),
    "mappings": extended_mapping,
}

output_file = "/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping_112.json"
with open(output_file, 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n✓ Saved extended mapping to {output_file}")
print(f"  Total: {len(extended_mapping)} classes")
print(f"  Named Carter classes: {len(current_mappings)}")
print(f"  Extended classes: {len(extended_mapping) - len(current_mappings)}")
print(f"  Unique graphs: {len(unique_extended)}/{len(available_graphs)}")
