"""
Generate all 74 Carter class images: orbital graph + Dynkin diagram(s).
"""
import json
import subprocess
from pathlib import Path

# Load the mapping
with open("/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping.json") as f:
    mapping_data = json.load(f)

mappings = mapping_data["mappings"]

# Dynkin diagram TikZ code for each type
def dynkin_diagram_tikz(levi_type, scale=0.4, shift_x=-4.0, shift_y=4.0):
    """Generate TikZ code for Dynkin diagram(s) based on Levi type."""

    # Parse the Levi type (e.g., "E7×A1", "D6×A2", "A1⁸")
    # For now, handle the main cases

    tikz_lines = []
    tikz_lines.append(r"\begin{scope}[shift={(" + f"{shift_x},{shift_y}" + r")},scale=" + f"{scale}]")

    if levi_type == "E8":
        # E8: o---o---o---o---o---o---o
        #                 |
        #                 o
        nodes = [(i, 0) for i in range(7)] + [(2, -1)]
        edges = [(i, i+1) for i in range(6)] + [(2, 7)]

    elif levi_type == "E7":
        # E7: o---o---o---o---o---o
        #             |
        #             o
        nodes = [(i, 0) for i in range(6)] + [(2, -1)]
        edges = [(i, i+1) for i in range(5)] + [(2, 6)]

    elif levi_type == "E6":
        # E6: o---o---o---o---o
        #         |
        #         o
        nodes = [(i, 0) for i in range(5)] + [(2, -1)]
        edges = [(i, i+1) for i in range(4)] + [(2, 5)]

    elif levi_type.startswith("E7") and "A1" in levi_type:
        # E7 + A1 (stacked)
        nodes = [(i, 0) for i in range(6)] + [(2, -1)] + [(7, 0)]  # E7 + isolated A1
        edges = [(i, i+1) for i in range(5)] + [(2, 6)]

    elif levi_type.startswith("E6") and "A2" in levi_type:
        # E6 + A2
        nodes = [(i, 0) for i in range(5)] + [(2, -1)] + [(i+6, 0) for i in range(2)]
        edges = [(i, i+1) for i in range(4)] + [(2, 5)] + [(6, 7)]

    elif levi_type.startswith("D"):
        # D_n: branching at the end
        rank = int(levi_type[1]) if len(levi_type) > 1 and levi_type[1].isdigit() else 4
        nodes = [(i, 0) for i in range(rank-1)] + [(rank-3, 1)]  # Branch
        edges = [(i, i+1) for i in range(rank-2)] + [(rank-3, rank-1)]

    elif levi_type.startswith("F4"):
        # F4: o---o==>o---o (double bond)
        nodes = [(i, 0) for i in range(4)]
        edges = [(0, 1), (1, 2), (2, 3)]  # Simplified

    elif levi_type.startswith("G2"):
        # G2: o≡≡≡o (triple bond)
        nodes = [(0, 0), (1, 0)]
        edges = [(0, 1)]

    elif levi_type.startswith("A"):
        # A_n: simple chain
        # Handle A1^8, A2, etc.
        if "A1" in levi_type and "A1A1" in levi_type:
            # Multiple A1's - show as disconnected nodes
            count = levi_type.count("A1")
            nodes = [(i*0.7, 0) for i in range(min(count, 8))]
            edges = []
        else:
            rank = int(levi_type[1]) if len(levi_type) > 1 and levi_type[1].isdigit() else 1
            nodes = [(i, 0) for i in range(rank)]
            edges = [(i, i+1) for i in range(rank-1)]

    else:
        # Default: single node
        nodes = [(0, 0)]
        edges = []

    # Draw edges
    for i, j in edges:
        x1, y1 = nodes[i]
        x2, y2 = nodes[j]
        tikz_lines.append(fr"\draw[white!60!bg,line width=0.6pt] ({x1},{y1})--({x2},{y2});")

    # Draw nodes
    for x, y in nodes:
        tikz_lines.append(fr"\fill[white!90!bg] ({x},{y}) circle (0.15);")

    tikz_lines.append(r"\end{scope}")

    return "\n".join(tikz_lines)


def create_carter_image(label, levi_type, graph_file):
    """Create a Carter class image: graph + Dynkin diagram."""

    output_file = f"carter_class_{label.replace('×', '_').replace('^', '')}.pdf"
    tex_file = f"carter_class_{label.replace('×', '_').replace('^', '')}.tex"

    # Build TikZ document
    tikz = [
        r"\documentclass[tikz,border=6pt]{standalone}",
        r"\usepackage{graphicx}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}",
        r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",

        # Include the graph as background, clipped to remove caption
        r"\begin{scope}",
        r"\clip (-5.5,-5.5) rectangle (5.5,5.5);",  # Clip to remove bottom caption
        r"\node at (0,0) {\includegraphics[width=10cm]{" + graph_file + r"}};",
        r"\end{scope}",

        # Add Dynkin diagram in top left
        dynkin_diagram_tikz(levi_type),

        r"\end{tikzpicture}",
        r"\end{document}",
    ]

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / tex_file
    with open(tex_path, 'w') as f:
        f.write("\n".join(tikz))

    return tex_path, output_file


# Generate all images
print(f"Generating {len(mappings)} Carter class images...\n")

generated = []
failed = []

for label, data in sorted(mappings.items()):  # ALL 74 classes
    levi = label  # For now, use label as Levi type
    graph = data["graph_file"]

    try:
        print(f"Creating {label}...", end=" ")
        tex_path, pdf_name = create_carter_image(label, levi, graph)

        # Compile immediately
        result = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', str(tex_path)],
            cwd='/Users/miriam/symmetric-graphs/e8',
            capture_output=True,
            timeout=30
        )

        if result.returncode == 0:
            print("✓")
            generated.append((label, tex_path, pdf_name))
        else:
            print("✗ (compile failed)")
            failed.append(label)
    except Exception as e:
        print(f"✗ ({e})")
        failed.append(label)

print(f"\n✓ Generated {len(generated)}/{len(mappings)} images")
if failed:
    print(f"✗ Failed: {len(failed)} images")
    print(f"  {', '.join(failed[:10])}" + ("..." if len(failed) > 10 else ""))

# Save list
with open("/Users/miriam/symmetric-graphs/e8/carter_images_generated.txt", "w") as f:
    for label, tex, pdf in generated:
        f.write(f"{label}\t{pdf}\n")

print(f"\n✓ Saved list to carter_images_generated.txt")
