"""
Generate all 112 Carter class images using the extended mapping.
Each image shows: orbital graph + Dynkin diagram (for named classes) or simple label (for extended classes).
"""
import json
import subprocess
from pathlib import Path

# Load the 112-class mapping
with open("/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping_112.json") as f:
    mapping_data = json.load(f)

mappings = mapping_data["mappings"]

print(f"Generating {len(mappings)} Carter class images...")
print(f"  Named classes: {mapping_data['named_carter_classes']}")
print(f"  Extended classes: {mapping_data['extended_classes']}")
print()

def dynkin_diagram_tikz(label, levi_type=None, scale=0.4, shift_x=-4.0, shift_y=4.0):
    """Generate TikZ code for Dynkin diagram(s) based on Carter label."""

    # For extended classes (class_N), just show a text label
    if label.startswith("class_"):
        tikz_lines = []
        tikz_lines.append(r"\begin{scope}[shift={(" + f"{shift_x},{shift_y}" + r")}]")
        # Escape underscores for LaTeX
        safe_label = label.replace("_", r"\_")
        tikz_lines.append(fr"\node[white!70!bg,font=\small] at (0,0) {{{safe_label}}};")
        tikz_lines.append(r"\end{scope}")
        return "\n".join(tikz_lines)

    # For named Carter classes, draw Dynkin diagram
    tikz_lines = []
    tikz_lines.append(r"\begin{scope}[shift={(" + f"{shift_x},{shift_y}" + r")},scale=" + f"{scale}]")

    # Parse the label to determine diagram structure
    if label == "1":
        # Identity - just a single node
        tikz_lines.append(fr"\node[white!70!bg,font=\small] at (0,0) {{id}};")

    elif label == "E8":
        # E8 Dynkin diagram
        nodes = [(i, 0) for i in range(7)] + [(2, -1)]
        edges = [(i, i+1) for i in range(6)] + [(2, 7)]

    elif label.startswith("E7"):
        # E7 Dynkin diagram
        nodes = [(i, 0) for i in range(6)] + [(2, -1)]
        edges = [(i, i+1) for i in range(5)] + [(2, 6)]

    elif label.startswith("E6"):
        # E6 Dynkin diagram
        nodes = [(i, 0) for i in range(5)] + [(2, -1)]
        edges = [(i, i+1) for i in range(4)] + [(2, 5)]

    elif label.startswith("D"):
        # D_n Dynkin diagram - branching at end
        # Extract rank from label
        rank = 8 if "D8" in label else 7 if "D7" in label else 6 if "D6" in label else 5 if "D5" in label else 4
        if rank >= 4:
            nodes = [(i, 0) for i in range(rank-1)] + [(rank-3, 1)]
            edges = [(i, i+1) for i in range(rank-2)] + [(rank-3, rank-1)]
        else:
            nodes = [(i, 0) for i in range(rank)]
            edges = [(i, i+1) for i in range(rank-1)]

    elif label.startswith("F4"):
        # F4 Dynkin diagram: o---o==>o---o
        nodes = [(i, 0) for i in range(4)]
        edges = [(0, 1), (1, 2), (2, 3)]

    elif label.startswith("G2"):
        # G2 Dynkin diagram: o≡≡≡o
        nodes = [(0, 0), (1, 0)]
        edges = [(0, 1)]

    elif label.startswith("A"):
        # A_n Dynkin diagram - simple chain
        # Handle different A-type labels
        if "A1A1A1A1A1A1A1A1" in label:
            # A1^8 - show 8 disconnected nodes
            nodes = [(i*0.7, 0) for i in range(8)]
            edges = []
        elif "A1" in label:
            # Count A1's
            count = label.count("A1")
            nodes = [(i*0.7, 0) for i in range(min(count, 8))]
            edges = []
        else:
            # Extract rank from first A_n
            import re
            match = re.search(r'A(\d+)', label)
            rank = int(match.group(1)) if match else 2
            nodes = [(i, 0) for i in range(rank)]
            edges = [(i, i+1) for i in range(rank-1)]

    else:
        # Default - single node with label
        tikz_lines.append(fr"\node[white!70!bg,font=\small] at (0,0) {{{label}}};")
        tikz_lines.append(r"\end{scope}")
        return "\n".join(tikz_lines)

    # Draw edges if we have a diagram
    if 'edges' in locals() and 'nodes' in locals():
        for i, j in edges:
            x1, y1 = nodes[i]
            x2, y2 = nodes[j]
            tikz_lines.append(fr"\draw[white!60!bg,line width=0.6pt] ({x1},{y1})--({x2},{y2});")

        # Draw nodes
        for x, y in nodes:
            tikz_lines.append(fr"\fill[white!90!bg] ({x},{y}) circle (0.15);")

    tikz_lines.append(r"\end{scope}")
    return "\n".join(tikz_lines)


def create_carter_image(label, graph_file):
    """Create a Carter class image: graph + Dynkin diagram or label."""

    output_file = f"carter_class_{label.replace('×', '_').replace('^', '').replace('⁸', '8')}.pdf"
    tex_file = f"carter_class_{label.replace('×', '_').replace('^', '').replace('⁸', '8')}.tex"

    # Build TikZ document
    tikz = [
        r"\documentclass[tikz,border=6pt]{standalone}",
        r"\usepackage{graphicx}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}",
        r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",

        # Include the graph as background
        r"\begin{scope}",
        r"\clip (-5.5,-5.5) rectangle (5.5,5.5);",
        r"\node at (0,0) {\includegraphics[width=10cm]{" + graph_file + r"}};",
        r"\end{scope}",

        # Add Dynkin diagram or label in top left
        dynkin_diagram_tikz(label),

        r"\end{tikzpicture}",
        r"\end{document}",
    ]

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / tex_file
    with open(tex_path, 'w') as f:
        f.write("\n".join(tikz))

    return tex_path, output_file


# Generate all images
generated = []
failed = []
skipped = []

for i, (label, data) in enumerate(sorted(mappings.items()), 1):
    graph = data["graph_file"]

    # Check if graph file exists
    graph_path = Path("/Users/miriam/symmetric-graphs/e8") / graph
    if not graph_path.exists():
        print(f"[{i:3d}/112] {label:20s} ✗ (graph {graph} not found)")
        skipped.append(label)
        continue

    try:
        print(f"[{i:3d}/112] {label:20s} ", end="")
        tex_path, pdf_name = create_carter_image(label, graph)

        # Compile with pdflatex
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

    # Progress update every 20 images
    if i % 20 == 0:
        print(f"  Progress: {i}/112 ({i*100//112}%)\n")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"✓ Generated: {len(generated)}/{len(mappings)} images")
if skipped:
    print(f"⊘ Skipped: {len(skipped)} (missing graphs)")
    print(f"  {', '.join(skipped[:5])}" + ("..." if len(skipped) > 5 else ""))
if failed:
    print(f"✗ Failed: {len(failed)} images")
    print(f"  {', '.join(failed[:10])}" + ("..." if len(failed) > 10 else ""))

# Save manifest
with open("/Users/miriam/symmetric-graphs/e8/carter_images_112_manifest.txt", "w") as f:
    f.write(f"# Generated {len(generated)} Carter class images\n")
    f.write(f"# Named classes: {mapping_data['named_carter_classes']}\n")
    f.write(f"# Extended classes: {mapping_data['extended_classes']}\n\n")
    for label, tex, pdf in generated:
        f.write(f"{label}\t{pdf}\n")

print(f"\n✓ Saved manifest to carter_images_112_manifest.txt")
