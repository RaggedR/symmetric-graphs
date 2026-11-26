"""
Dynkin diagram TikZ generator for Carter classes.
Handles composite types like A2A2A2, D4A3A1, etc.
"""

import re


def parse_levi_type(levi_str):
    """
    Parse a Levi type string into components.

    Examples:
        "E8" -> [("E", 8)]
        "A2A2A2" -> [("A", 2), ("A", 2), ("A", 2)]
        "D4A3A1" -> [("D", 4), ("A", 3), ("A", 1)]
        "E7A1" -> [("E", 7), ("A", 1)]
        "A1A1A1A1A1A1A1A1" -> [("A", 1)] * 8
    """
    if levi_str in ["1", "class_1"] or levi_str.startswith("class_"):
        return []  # Trivial or unknown class

    components = []

    # Pattern: Letter followed by optional digit(s)
    # Handles E8, A2, D4, etc.
    pattern = r'([A-Z])(\d+)'

    matches = re.findall(pattern, levi_str)
    for letter, rank_str in matches:
        rank = int(rank_str)
        components.append((letter, rank))

    return components


def single_dynkin_diagram(root_type, rank, scale=0.4):
    """
    Generate TikZ code for a single Dynkin diagram.
    Returns (tikz_lines, width) where width is approximate diagram width.
    """
    tikz_lines = []

    if root_type == "A":
        # A_n: o---o---o---...---o (simple chain of n nodes)
        nodes = [(i, 0) for i in range(rank)]
        edges = [(i, i+1) for i in range(rank-1)]
        width = rank - 1 if rank > 1 else 0.5

    elif root_type == "D":
        # D_n: branching at the end
        #      o
        #     /
        # o---o---...---o
        #     \
        #      o
        if rank >= 4:
            nodes = [(i, 0) for i in range(rank-2)] + [(rank-3, 1), (rank-3, -1)]
            edges = [(i, i+1) for i in range(rank-3)] + [(rank-3, rank-2), (rank-3, rank-1)]
            width = rank - 2
        else:
            # Degenerate case
            nodes = [(i, 0) for i in range(rank)]
            edges = [(i, i+1) for i in range(rank-1)]
            width = rank - 1

    elif root_type == "E":
        # E_6, E_7, E_8: main chain with branch
        if rank == 6:
            # E6: o---o---o---o---o
            #         |
            #         o
            nodes = [(i, 0) for i in range(5)] + [(2, -1)]
            edges = [(i, i+1) for i in range(4)] + [(2, 5)]
            width = 4
        elif rank == 7:
            # E7: o---o---o---o---o---o
            #             |
            #             o
            nodes = [(i, 0) for i in range(6)] + [(2, -1)]
            edges = [(i, i+1) for i in range(5)] + [(2, 6)]
            width = 5
        elif rank == 8:
            # E8: o---o---o---o---o---o---o
            #                 |
            #                 o
            nodes = [(i, 0) for i in range(7)] + [(3, -1)]
            edges = [(i, i+1) for i in range(6)] + [(3, 7)]
            width = 6
        else:
            # Unknown E type
            nodes = [(0, 0)]
            edges = []
            width = 0.5

    elif root_type == "F":
        # F4: o---o==>o---o (double bond between 2nd and 3rd)
        if rank == 4:
            nodes = [(i, 0) for i in range(4)]
            edges = [(0, 1), (2, 3)]  # Single bonds on the sides
            # Draw double bond as two parallel lines (more visible)
            offset = 0.06
            tikz_lines.append(fr"\draw[white!60!bg,line width=0.5pt] (1,{-offset})--(2,{-offset});")  # Bottom
            tikz_lines.append(fr"\draw[white!60!bg,line width=0.5pt] (1,{offset})--(2,{offset});")  # Top
            width = 3
        else:
            nodes = [(i, 0) for i in range(rank)]
            edges = [(i, i+1) for i in range(rank-1)]
            width = rank - 1

    elif root_type == "G":
        # G2: o≡≡≡o (triple bond - draw as 3 parallel lines)
        if rank == 2:
            nodes = [(0, 0), (1.2, 0)]
            # Draw triple bond as three parallel lines (more visible than TikZ double)
            offset = 0.08  # Spacing between lines
            tikz_lines.append(fr"\draw[white!60!bg,line width=0.5pt] (0,{-offset})--(1.2,{-offset});")  # Bottom
            tikz_lines.append(fr"\draw[white!60!bg,line width=0.5pt] (0,0)--(1.2,0);")  # Middle
            tikz_lines.append(fr"\draw[white!60!bg,line width=0.5pt] (0,{offset})--(1.2,{offset});")  # Top
            edges = []  # Already drawn
            width = 1.2
        else:
            nodes = [(i, 0) for i in range(rank)]
            edges = [(i, i+1) for i in range(rank-1)]
            width = rank - 1
    else:
        # Unknown type - single node
        nodes = [(0, 0)]
        edges = []
        width = 0.5

    # Draw edges (unless already drawn for special bonds)
    for i, j in edges:
        x1, y1 = nodes[i]
        x2, y2 = nodes[j]
        tikz_lines.append(fr"\draw[white!60!bg,line width=0.6pt] ({x1},{y1})--({x2},{y2});")

    # Draw nodes
    for x, y in nodes:
        tikz_lines.append(fr"\fill[white!90!bg] ({x},{y}) circle (0.15);")

    return tikz_lines, width


def calculate_diagram_height(levi_type):
    """
    Calculate the vertical extent of the Dynkin diagram.
    Returns the height in coordinate units (before scaling).
    """
    components = parse_levi_type(levi_type)

    if not components:
        return 0.5  # Small label

    num_components = len(components)

    # Determine spacing
    if num_components >= 6:
        vertical_spacing = 1.3
    elif num_components >= 4:
        vertical_spacing = 1.6
    else:
        vertical_spacing = 2.0

    # Height is number of components times spacing, plus some for the diagram itself
    # Each diagram is about 1.5 units tall (including branches)
    total_height = (num_components - 1) * vertical_spacing + 1.5

    # Apply scale of 0.4
    return total_height * 0.4


def dynkin_diagram_tikz(levi_type, scale=0.4, shift_x=-4.0, shift_y=4.0, rotate=0):
    """
    Generate TikZ code for Dynkin diagram(s) based on Levi type.
    Handles composite types by showing multiple disjoint diagrams.

    Args:
        levi_type: String like "E8", "A2A2A2", "D4A3A1", etc.
        scale: Overall scale factor (kept constant at 0.4 for consistency)
        shift_x, shift_y: Position of top-left corner

    Returns:
        TikZ code as string
    """
    components = parse_levi_type(levi_type)

    if not components:
        # Trivial or unknown - just show label
        tikz_lines = []
        tikz_lines.append(r"\begin{scope}[shift={(" + f"{shift_x},{shift_y}" + r")}]")
        safe_label = levi_type.replace("_", r"\_")
        tikz_lines.append(fr"\node[white!70!bg,font=\small] at (0,0) {{{safe_label}}};")
        tikz_lines.append(r"\end{scope}")
        return "\n".join(tikz_lines)

    # Adaptive spacing based on number of components
    # Keep scale constant at 0.4, only adjust spacing
    num_components = len(components)
    if num_components >= 6:
        # Many components (e.g., A1^8) - use tighter spacing
        vertical_spacing = 1.3
    elif num_components >= 4:
        # Several components - moderate spacing
        vertical_spacing = 1.6
    else:
        # Few components - normal spacing
        vertical_spacing = 2.0

    tikz_lines = []
    # Add rotation if specified (for distinguishing dual/related classes)
    if rotate != 0:
        tikz_lines.append(r"\begin{scope}[shift={(" + f"{shift_x},{shift_y}" + r")},scale=" + f"{scale}" + r",rotate=" + f"{rotate}]")
    else:
        tikz_lines.append(r"\begin{scope}[shift={(" + f"{shift_x},{shift_y}" + r")},scale=" + f"{scale}]")

    # Layout components vertically with spacing
    current_y = 0

    for root_type, rank in components:
        # Generate single diagram
        tikz_lines.append(fr"\begin{{scope}}[shift={{(0,{current_y})}}]")
        diagram_lines, width = single_dynkin_diagram(root_type, rank, scale)
        tikz_lines.extend(diagram_lines)
        tikz_lines.append(r"\end{scope}")

        # Move down for next component
        current_y -= vertical_spacing

    tikz_lines.append(r"\end{scope}")

    return "\n".join(tikz_lines)


if __name__ == "__main__":
    # Test parsing
    test_cases = [
        "E8",
        "A2",
        "A2A2A2",
        "D4A3A1",
        "E7A1",
        "A1A1A1A1A1A1A1A1",
        "F4",
        "G2",
        "D6A2",
        "E6A1A1",
    ]

    print("Testing Levi type parsing:")
    for levi in test_cases:
        components = parse_levi_type(levi)
        print(f"  {levi:20s} -> {components}")
