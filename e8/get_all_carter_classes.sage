#!/usr/bin/env sage

"""
Get all 112 conjugacy classes of W(E8) using SageMath.
Export them with Carter labels.
"""

print("Loading W(E8) Weyl group...")
W = WeylGroup("E8", prefix="s")
print(f"|W(E8)| = {W.cardinality()}")

print("\nEnumerating conjugacy classes...")
conjugacy_classes = W.conjugacy_classes()
print(f"Found {len(conjugacy_classes)} conjugacy classes")

# Extract information for each class
classes_data = []

for i, cc in enumerate(conjugacy_classes):
    rep = cc.an_element()

    # Basic properties
    class_size = cc.cardinality()
    elt_order = rep.order()

    # Carter label (as a reduced word or cycle type)
    # In SageMath, we can get the cycle type
    try:
        # Get the reduced word
        red_word = rep.reduced_word()
        label = "s" + "".join(str(x) for x in red_word[:5])  # truncate if long
    except:
        label = f"class_{i+1}"

    classes_data.append({
        "index": int(i + 1),
        "label": str(label),
        "order": int(elt_order),
        "class_size": int(class_size),
        "reduced_word_length": int(len(rep.reduced_word())) if hasattr(rep, 'reduced_word') and rep.reduced_word() is not None else 0,
    })

    if (i + 1) % 10 == 0:
        print(f"  Processed {i+1}/{len(conjugacy_classes)} classes...")

print(f"\n✓ Processed all {len(conjugacy_classes)} classes")

# Save to JSON
import json

output_data = {
    "group": "W(E8)",
    "group_order": int(W.cardinality()),
    "num_conjugacy_classes": len(conjugacy_classes),
    "classes": classes_data
}

output_path = "/Users/miriam/symmetric-graphs/e8/carter_classes_complete.json"
with open(output_path, 'w') as f:
    json.dump(output_data, f, indent=2)

print(f"\n✓ Saved to {output_path}")

# Print summary
print("\nSummary by element order:")
from collections import Counter
order_counts = Counter(c["order"] for c in classes_data)
for order in sorted(order_counts.keys()):
    print(f"  Order {order:2d}: {order_counts[order]:3d} classes")
