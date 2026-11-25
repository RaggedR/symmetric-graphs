"""
Enumerate possible Carter class labels for W(E8).

Carter classes are labeled by Levi subgroups, which are products of
irreducible root systems of total rank ≤ 8.

Not all combinations occur - we need to verify against Carter's tables.
"""

def generate_partitions_with_types():
    """
    Generate all possible ways to partition rank 8 into irreducible components.
    Format: list of (type, rank) pairs
    """
    types = ['A', 'D', 'E']  # Main types (B, C are for different contexts)

    results = []

    # Single exceptional types
    results.append([('E', 8)])  # E8
    results.append([('E', 7), ('A', 1)])  # E7A1
    results.append([('E', 6), ('A', 2)])  # E6A2
    results.append([('E', 6), ('A', 1), ('A', 1)])  # E6A1²

    # D-series dominance
    for d_rank in range(4, 9):  # D4 through D8
        remaining = 8 - d_rank
        results.append([('D', d_rank)])  # Just Dn

        # Dn + A-series
        if remaining > 0:
            # Single A
            if remaining >= 1:
                results.append([('D', d_rank), ('A', remaining)])

            # Two A's
            for a1 in range(1, remaining):
                a2 = remaining - a1
                if a1 <= a2:  # Avoid duplicates
                    results.append([('D', d_rank), ('A', a2), ('A', a1)])

            # Three A's
            if remaining >= 3:
                for a1 in range(1, remaining - 1):
                    for a2 in range(a1, remaining - a1):
                        a3 = remaining - a1 - a2
                        if a2 <= a3:
                            results.append([('D', d_rank), ('A', a3), ('A', a2), ('A', a1)])

    # Pure A-series
    def partition_rank(n, max_part=None):
        """Generate integer partitions of n"""
        if max_part is None:
            max_part = n
        if n == 0:
            yield []
            return
        for i in range(min(n, max_part), 0, -1):
            for p in partition_rank(n - i, i):
                yield [i] + p

    for partition in partition_rank(8):
        if len(partition) > 0:
            components = [('A', r) for r in partition]
            results.append(components)

    # E6 + small combinations
    results.append([('E', 6)])
    results.append([('E', 6), ('A', 1), ('A', 1)])  # E6A1²

    # E7 combinations
    results.append([('E', 7)])
    results.append([('E', 7), ('A', 1)])

    # F4 (rank 4)
    results.append([('F', 4), ('A', 4)])
    results.append([('F', 4), ('A', 3), ('A', 1)])
    results.append([('F', 4), ('A', 2), ('A', 2)])
    results.append([('F', 4), ('A', 2), ('A', 1), ('A', 1)])
    results.append([('F', 4), ('D', 4)])

    # G2 (rank 2)
    results.append([('G', 2), ('A', 6)])
    results.append([('G', 2), ('D', 6)])
    results.append([('G', 2), ('E', 6)])
    results.append([('G', 2), ('A', 3), ('A', 3)])

    # De-duplicate
    unique = []
    seen = set()
    for r in results:
        # Sort components by (type, -rank) for canonical form
        canonical = tuple(sorted(r, key=lambda x: (x[0], -x[1])))
        if canonical not in seen:
            seen.add(canonical)
            unique.append(canonical)

    return sorted(unique, key=lambda x: (-sum(r for t, r in x), x))


def format_label(components):
    """Format a list of components as a Carter label"""
    from collections import Counter

    # Group by (type, rank)
    counts = Counter(components)

    # Sort: E first, then D, then A, by decreasing rank
    priority = {'E': 0, 'F': 1, 'G': 2, 'D': 3, 'A': 4}
    sorted_items = sorted(counts.items(), key=lambda x: (priority[x[0][0]], -x[0][1]))

    parts = []
    for (typ, rank), count in sorted_items:
        part = f"{typ}{rank}"
        if count > 1:
            part += f"^{count}" if count <= 3 else f"^{count}"
        parts.append(part)

    return "".join(parts)


if __name__ == "__main__":
    possibilities = generate_partitions_with_types()

    print(f"Generated {len(possibilities)} possible Carter class labels")
    print("(Note: Not all of these actually occur in W(E8)!)\n")

    # Group by dominant type
    exceptional = []
    d_dominant = []
    a_dominant = []
    mixed = []

    for components in possibilities:
        label = format_label(components)

        if any(t == 'E' for t, r in components):
            exceptional.append(label)
        elif any(t == 'D' for t, r in components):
            d_dominant.append(label)
        elif all(t == 'A' for t, r in components):
            a_dominant.append(label)
        else:
            mixed.append(label)

    print("## Exceptional Classes:")
    for label in exceptional[:20]:
        print(f"  {label}")
    print(f"  ... ({len(exceptional)} total)\n")

    print("## D-Dominant Classes:")
    for label in d_dominant[:20]:
        print(f"  {label}")
    print(f"  ... ({len(d_dominant)} total)\n")

    print("## Pure A-Series:")
    for label in a_dominant[:20]:
        print(f"  {label}")
    print(f"  ... ({len(a_dominant)} total)\n")

    print(f"\nTotal candidates: {len(possibilities)}")
    print("Actual W(E8) has exactly 112 conjugacy classes")
