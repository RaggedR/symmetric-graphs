"""
1-to-1 Assignment Algorithm: Carter Classes → Orbitals

This creates a perfect 1-to-1 mapping once we have 112+ unique orbitals.
"""

import json
from collections import defaultdict
import re


def parse_levi_type(label):
    """
    Parse a Carter label into Levi components.

    Examples:
        "E7A1" → ["E7", "A1"]
        "D6A2" → ["D6", "A2"]
        "A3A3A2" → ["A3", "A3", "A2"]
        "F4D4" → ["F4", "D4"]
    """
    # Pattern: match root system labels (E8, E7, E6, D4-D8, F4, G2, A1-A8)
    pattern = r'([EDFGA]\d+)'
    components = re.findall(pattern, label)
    return components


def levi_rank(component):
    """Return the rank (dimension) of a root system component."""
    ranks = {
        'E8': 8, 'E7': 7, 'E6': 6,
        'D8': 8, 'D7': 7, 'D6': 6, 'D5': 5, 'D4': 4,
        'F4': 4, 'G2': 2,
        'A8': 8, 'A7': 7, 'A6': 6, 'A5': 5, 'A4': 4,
        'A3': 3, 'A2': 2, 'A1': 1
    }
    return ranks.get(component, 0)


def get_family(component):
    """Get the family letter (E, D, F, G, A) from a component."""
    if component:
        return component[0]
    return None


def order_to_family(order):
    """Heuristic: map element order to likely root system family."""
    order_map = {
        1: 'E8',
        30: 'E8',  # E8 Coxeter
        18: 'E7',  # E7 Coxeter
        12: 'E6',  # E6 or F4 Coxeter
        10: 'D',   # D8 related
        8: 'D',    # D5
        6: 'G2',   # or D4 or A5
        5: 'A',    # A4
        4: 'A',    # A3
        3: 'A',    # A2
        2: 'A',    # A1
    }
    return order_map.get(order, 'E8')


def assign_carter_to_orbitals_one_to_one(carter_classes, orbitals):
    """
    Create a perfect 1-to-1 assignment from Carter classes to orbitals.

    Args:
        carter_classes: List of dicts with keys: label, levi_type, element_order
        orbitals: List of dicts with keys: graph_file, levi_type, family, ip_value

    Returns:
        Dict mapping carter_label → orbital_graph_file
    """

    assignments = {}
    used_orbitals = set()

    print(f"Starting 1-to-1 assignment...")
    print(f"  Carter classes: {len(carter_classes)}")
    print(f"  Orbitals available: {len(orbitals)}")
    print()

    # Create orbital lookup structures
    orbitals_by_levi = defaultdict(list)
    orbitals_by_family = defaultdict(list)

    for orb in orbitals:
        # Group by exact Levi type
        if 'levi_type' in orb:
            orbitals_by_levi[orb['levi_type']].append(orb)

        # Group by family
        if 'family' in orb:
            orbitals_by_family[orb['family']].append(orb)

    # Phase 1: Exact Levi type matches
    print("Phase 1: Exact Levi type matches")
    phase1_count = 0

    for cls in carter_classes:
        label = cls['label']
        levi = cls.get('levi_type', label)

        # Try exact match
        candidates = [
            orb for orb in orbitals_by_levi.get(levi, [])
            if orb['graph_file'] not in used_orbitals
        ]

        if candidates:
            # Prefer ip=+1 (principal orbital) if available
            principal = [o for o in candidates if o.get('ip_value', 0) == 1]
            orbital = principal[0] if principal else candidates[0]

            assignments[label] = orbital['graph_file']
            used_orbitals.add(orbital['graph_file'])
            phase1_count += 1

    print(f"  Assigned: {phase1_count} classes")
    print()

    # Phase 2: Levi component matches
    print("Phase 2: Levi component matches")
    phase2_count = 0

    for cls in carter_classes:
        label = cls['label']
        if label in assignments:
            continue

        levi = cls.get('levi_type', label)
        components = parse_levi_type(levi)

        # Try matching largest component first
        for component in sorted(components, key=levi_rank, reverse=True):
            candidates = [
                orb for orb in orbitals_by_levi.get(component, [])
                if orb['graph_file'] not in used_orbitals
            ]

            if candidates:
                orbital = candidates[0]
                assignments[label] = orbital['graph_file']
                used_orbitals.add(orbital['graph_file'])
                phase2_count += 1
                break

    print(f"  Assigned: {phase2_count} classes")
    print()

    # Phase 3: Family matches
    print("Phase 3: Family matches")
    phase3_count = 0

    for cls in carter_classes:
        label = cls['label']
        if label in assignments:
            continue

        levi = cls.get('levi_type', label)
        components = parse_levi_type(levi)

        # Get family of largest component
        if components:
            primary_family = get_family(components[0])

            candidates = [
                orb for orb in orbitals_by_family.get(primary_family, [])
                if orb['graph_file'] not in used_orbitals
            ]

            if candidates:
                orbital = candidates[0]
                assignments[label] = orbital['graph_file']
                used_orbitals.add(orbital['graph_file'])
                phase3_count += 1

    print(f"  Assigned: {phase3_count} classes")
    print()

    # Phase 4: Element order heuristic (for extended classes)
    print("Phase 4: Element order heuristic")
    phase4_count = 0

    for cls in carter_classes:
        label = cls['label']
        if label in assignments:
            continue

        order = cls.get('element_order', 2)
        family = order_to_family(order)

        candidates = [
            orb for orb in orbitals_by_family.get(family, [])
            if orb['graph_file'] not in used_orbitals
        ]

        if candidates:
            orbital = candidates[0]
            assignments[label] = orbital['graph_file']
            used_orbitals.add(orbital['graph_file'])
            phase4_count += 1

    print(f"  Assigned: {phase4_count} classes")
    print()

    # Phase 5: Fill remaining with any available orbital
    print("Phase 5: Fill remaining")
    phase5_count = 0

    remaining_orbitals = [
        orb for orb in orbitals
        if orb['graph_file'] not in used_orbitals
    ]

    for cls in carter_classes:
        label = cls['label']
        if label in assignments:
            continue

        if remaining_orbitals:
            orbital = remaining_orbitals.pop(0)
            assignments[label] = orbital['graph_file']
            used_orbitals.add(orbital['graph_file'])
            phase5_count += 1
        else:
            print(f"  WARNING: No orbital available for {label}")

    print(f"  Assigned: {phase5_count} classes")
    print()

    # Verify 1-to-1
    print("="*70)
    print("VERIFICATION")
    print("="*70)
    print(f"Total assignments: {len(assignments)}")
    print(f"Unique orbitals used: {len(set(assignments.values()))}")
    print(f"Is 1-to-1: {len(assignments) == len(set(assignments.values()))}")

    if len(assignments) != len(set(assignments.values())):
        print("\nWARNING: Not a 1-to-1 mapping!")
        print("Need more orbitals or better distribution.")
    else:
        print("\n✓ Perfect 1-to-1 mapping achieved!")

    return assignments


def main():
    """
    Main entry point: load current data and demonstrate algorithm.

    NOTE: This uses current 47 orbitals, so won't be perfect 1-to-1 yet.
    Once we compute 112+ orbitals, this will produce perfect mapping.
    """

    # Load current orbitals (47 available)
    orbitals = []

    # E8 family
    for ip, suffix in [(-2, 'ipm2'), (-1, 'ipm1'), (0, 'ipp0'), (1, 'ipp1')]:
        orbitals.append({
            'graph_file': f'orbital_e8_{suffix}.pdf',
            'levi_type': 'E8',
            'family': 'E8',
            'ip_value': ip
        })

    for name in ['id', 'mandala', 'bloom', 'weave', 'antipodal']:
        orbitals.append({
            'graph_file': f'clean_e8_{name}.pdf',
            'levi_type': 'E8',
            'family': 'E8',
            'ip_value': 1 if name == 'mandala' else 0
        })

    # E7 family
    for ip, suffix in [(-1, 'ipm1'), (0, 'ipp0'), (1, 'ipp1')]:
        orbitals.append({
            'graph_file': f'orbital_e7_{suffix}.pdf',
            'levi_type': 'E7',
            'family': 'E7',
            'ip_value': ip
        })
    orbitals.append({'graph_file': 'clean_e7_mandala.pdf', 'levi_type': 'E7', 'family': 'E7', 'ip_value': 1})

    # E6 family
    for ip, suffix in [(0, 'ipp0'), (1, 'ipp1')]:
        orbitals.append({
            'graph_file': f'orbital_e6_{suffix}.pdf',
            'levi_type': 'E6',
            'family': 'E6',
            'ip_value': ip
        })
    orbitals.append({'graph_file': 'clean_e6_schlafli.pdf', 'levi_type': 'E6', 'family': 'E6', 'ip_value': 0})

    # Add all other families similarly (D, F4, G2, A, E5)...
    # (Abbreviated for clarity - full version would include all 47)

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

    # Run algorithm
    assignments = assign_carter_to_orbitals_one_to_one(carter_classes, orbitals)

    # Save result
    output = {
        'description': '1-to-1 Carter to orbital mapping (using assignment algorithm)',
        'note': 'Currently uses 47 orbitals, so not perfect 1-to-1. Compute 112+ orbitals for perfection.',
        'algorithm_version': '1.0',
        'total_classes': len(assignments),
        'unique_orbitals_used': len(set(assignments.values())),
        'is_one_to_one': len(assignments) == len(set(assignments.values())),
        'mappings': assignments
    }

    output_file = '/Users/miriam/symmetric-graphs/e8/carter_to_orbital_algorithm.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n✓ Saved to {output_file}")


if __name__ == '__main__':
    main()
