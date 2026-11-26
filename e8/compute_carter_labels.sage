#!/usr/bin/env sage

"""
Compute Carter labels (Levi types) for all 112 conjugacy classes of W(E8).

The Carter label of a conjugacy class is determined by the Levi type of its
centralizer, which is a parabolic subgroup.
"""

print("Loading W(E8)...")
W = WeylGroup("E8", prefix="s")
print(f"|W(E8)| = {W.cardinality()}\n")

print("Getting conjugacy classes...")
cc_list = W.conjugacy_classes_representatives()
print(f"Found {len(cc_list)} representatives\n")

print("Computing Carter labels (this will take a while)...")
print("Progress:")

carter_data = []

for i, w in enumerate(cc_list):
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/112...")

    # Get basic properties
    order = w.order()

    # Try to get Carter label
    # In Sage, we can use the standard Carter labeling
    try:
        # The Carter label can be extracted from the support
        # (the set of simple reflections appearing in minimal factorizations)
        supp = set()
        red = w.reduced_word()
        for j in red:
            supp.add(j)

        # This is a rough approximation - actual Carter labels require
        # analyzing the cent ralizer structure
        carter_label = f"class_{i+1}_order{order}"

    except Exception as e:
        carter_label = f"class_{i+1}"

    carter_data.append({
        "index": int(i + 1),
        "carter_label": str(carter_label),
        "element_order": int(order),
        "reduced_word": [int(x) for x in w.reduced_word()],
        "length": int(w.length()),
    })

print(f"\n✓ Processed all {len(carter_data)} classes\n")

# Save
import json

output = {
    "group": "W(E8)",
    "num_classes": len(carter_data),
    "classes": carter_data,
    "note": "Carter labels are placeholder - full computation requires centralizer analysis"
}

output_file = "/Users/miriam/symmetric-graphs/e8/carter_labels_preliminary.json"
with open(output_file, 'w') as f:
    json.dump(output, f, indent=2)

print(f"✓ Saved to {output_file}")
print("\nNOTE: To get true Carter labels, we need to:")
print("1. Compute centralizer C_W(w) for each representative w")
print("2. Identify C_W(w) as a parabolic W(L)")
print("3. Extract the Levi type L")
print("\nThis requires deeper Coxeter group machinery...")
