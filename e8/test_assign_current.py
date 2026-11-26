"""
Test the assignment algorithm with current 60 orbitals.
This will show how much better the distribution is.
"""

import json
import os
from collections import defaultdict, Counter

# Import the assignment functions
import sys
sys.path.insert(0, '/Users/miriam/symmetric-graphs/e8')
from assign_one_to_one import (
    parse_levi_type, levi_rank, get_family, order_to_family,
    assign_carter_to_orbitals_one_to_one
)


def collect_all_available_orbitals():
    """Collect all currently available orbital PDFs."""

    orbitals = []

    # Scan directory for all orbital PDFs
    base_dir = '/Users/miriam/symmetric-graphs/e8'

    for filename in os.listdir(base_dir):
        if not filename.endswith('.pdf'):
            continue

        if filename.startswith('clean_') or filename.startswith('orbital_'):
            # Parse filename to extract info
            # Examples: clean_e8_mandala.pdf, orbital_e8_ipp1.pdf, orbital_a4_ipm2.pdf

            name = filename.replace('.pdf', '')

            # Determine family and type
            if 'e8' in name:
                family = 'E8'
                levi_type = 'E8'
            elif 'e7' in name:
                family = 'E7'
                levi_type = 'E7'
            elif 'e6' in name:
                family = 'E6'
                levi_type = 'E6'
            elif 'e5' in name or 'd5' in name:
                family = 'E5'
                levi_type = 'D5'
            elif 'd7' in name:
                family = 'D'
                levi_type = 'D7'
            elif 'd6' in name:
                family = 'D'
                levi_type = 'D6'
            elif 'd4' in name:
                family = 'D'
                levi_type = 'D4'
            elif 'f4' in name:
                family = 'F4'
                levi_type = 'F4'
            elif 'g2' in name:
                family = 'G2'
                levi_type = 'G2'
            elif 'a5' in name:
                family = 'A'
                levi_type = 'A5'
            elif 'a4' in name:
                family = 'A'
                levi_type = 'A4'
            elif 'a3' in name:
                family = 'A'
                levi_type = 'A3'
            elif 'a2' in name:
                family = 'A'
                levi_type = 'A2'
            elif 'a1' in name:
                family = 'A'
                levi_type = 'A1'
            else:
                family = 'unknown'
                levi_type = 'unknown'

            # Determine ip value
            if 'ipp2' in name:
                ip_value = 2
            elif 'ipp1' in name:
                ip_value = 1
            elif 'ipp0' in name:
                ip_value = 0
            elif 'ipm1' in name:
                ip_value = -1
            elif 'ipm2' in name:
                ip_value = -2
            else:
                ip_value = 1  # default for clean graphs

            orbitals.append({
                'graph_file': filename,
                'levi_type': levi_type,
                'family': family,
                'ip_value': ip_value
            })

    return orbitals


def main():
    """Test assignment algorithm with current orbitals."""

    print("="*70)
    print("Testing Assignment Algorithm with Current 60 Orbitals")
    print("="*70)
    print()

    # Collect available orbitals
    print("Collecting available orbitals...")
    orbitals = collect_all_available_orbitals()

    print(f"Found {len(orbitals)} orbital PDFs")
    print()

    # Show distribution by family
    family_counts = Counter(o['family'] for o in orbitals)
    print("Orbitals by family:")
    for family, count in sorted(family_counts.items(), key=lambda x: -x[1]):
        print(f"  {family:10s}: {count:3d} orbitals")
    print()

    # Load 112 Carter classes
    with open('/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping_112.json') as f:
        data = json.load(f)

    carter_classes = []
    for label, info in data['mappings'].items():
        carter_classes.append({
            'label': label,
            'levi_type': label if not label.startswith('class_') else info.get('hecke_family', 'unknown'),
            'element_order': info.get('element_order', 2)
        })

    print(f"Loaded {len(carter_classes)} Carter classes")
    print()

    # Run assignment algorithm
    print("="*70)
    print("Running 1-to-1 Assignment Algorithm")
    print("="*70)
    print()

    assignments = assign_carter_to_orbitals_one_to_one(carter_classes, orbitals)

    # Analyze results
    print()
    print("="*70)
    print("ANALYSIS")
    print("="*70)
    print()

    # Count usage
    usage_counts = Counter(assignments.values())

    print(f"Orbital usage statistics:")
    print(f"  Total assignments: {len(assignments)}")
    print(f"  Unique orbitals used: {len(usage_counts)}")
    print(f"  Orbitals available but unused: {len(orbitals) - len(usage_counts)}")
    print()

    # Show most reused
    print("Most reused orbitals:")
    for graph, count in usage_counts.most_common(10):
        classes = [k for k, v in assignments.items() if v == graph]
        print(f"  {graph}: {count} classes")
        if count <= 3:
            print(f"    {', '.join(classes)}")
    print()

    # Compare to old mapping
    old_usage = Counter(info['graph_file'] for info in data['mappings'].values())

    print("Comparison with old mapping:")
    print(f"  Old: {len(set(old_usage.keys()))} unique graphs, max reuse = {max(old_usage.values())}")
    print(f"  New: {len(usage_counts)} unique graphs, max reuse = {max(usage_counts.values())}")
    print(f"  Improvement: {len(usage_counts) - len(set(old_usage.keys()))} more unique graphs used")
    print()

    # Save new mapping
    output = {
        'description': 'Improved Carter-to-Hecke mapping using 60 orbitals and 1-to-1 algorithm',
        'algorithm_version': '1.0',
        'total_classes': len(assignments),
        'orbitals_available': len(orbitals),
        'unique_orbitals_used': len(usage_counts),
        'max_reuse': max(usage_counts.values()),
        'is_one_to_one': len(assignments) == len(usage_counts),
        'mappings': assignments
    }

    output_file = '/Users/miriam/symmetric-graphs/e8/carter_to_orbital_improved.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"✓ Saved improved mapping to {output_file}")
    print()

    return assignments, orbitals


if __name__ == '__main__':
    main()
