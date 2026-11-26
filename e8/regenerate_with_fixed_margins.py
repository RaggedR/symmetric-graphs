"""
Regenerate Carter class images with FIXED MARGINS and adaptive graph sizing.
Ensures Dynkin diagrams have consistent margins from image edges.
"""

import json
import subprocess
from pathlib import Path
from dynkin_diagrams import dynkin_diagram_tikz, calculate_diagram_height


def create_carter_image_fixed_margins(label, graph_file, hecke_family=None):
    """Create a Carter class image with fixed margins and adaptive graph sizing."""

    output_file = f"carter_class_{label.replace('×', '_').replace('^', '').replace('⁸', '8')}.pdf"
    tex_file = f"carter_class_{label.replace('×', '_').replace('^', '').replace('⁸', '8')}.tex"

    # For class_* labels, use hecke_family for the diagram
    diagram_label = hecke_family if (label.startswith('class_') and hecke_family) else label

    # Fixed parameters
    CANVAS_SIZE = 11.0  # Total canvas is 11x11 units (-5.5 to 5.5)
    DIAGRAM_MARGIN_TOP = 0.6   # Margin from top edge to diagram (increased for visibility)
    DIAGRAM_MARGIN_LEFT = 0.6  # Margin from left edge to diagram (increased for visibility)

    # Calculate Dynkin diagram height
    diagram_height = calculate_diagram_height(diagram_label)

    # Position diagram with FIXED margins from edge
    # Top edge is at y=5.5, so diagram starts at 5.5 - DIAGRAM_MARGIN_TOP
    diagram_y = 5.5 - DIAGRAM_MARGIN_TOP
    # Left edge is at x=-5.5, so diagram starts at -5.5 + DIAGRAM_MARGIN_LEFT
    diagram_x = -5.5 + DIAGRAM_MARGIN_LEFT

    # Calculate graph size to avoid overlap
    # Diagram bottom edge in canvas coordinates
    diagram_bottom = diagram_y - diagram_height

    # Default graph is 10cm (fills most of 11-unit canvas)
    graph_width_cm = 10.0

    # If diagram is very tall and extends below y=1.0, shrink graph
    if diagram_bottom < 1.0:
        # Graph is centered at (0,0)
        # We want graph top to stay below diagram bottom with 0.5 unit buffer
        max_allowed_top = diagram_bottom - 0.5
        # Graph of width W cm extends to approximately ±W/2 units
        # So we want W/2 < max_allowed_top
        graph_width_cm = max_allowed_top * 2.0
        # Clamp to reasonable range
        graph_width_cm = max(7.0, min(graph_width_cm, 10.0))

    # Determine rotation for class_* to distinguish duals
    # Rotate specific classes to create visual distinction
    rotation = 0
    if label in ['class_6', 'class_10', 'class_22']:
        rotation = 90  # Rotate these to distinguish from their duals

    # Build TikZ document
    tikz = [
        r"\documentclass[tikz,border=6pt]{standalone}",
        r"\usepackage{graphicx}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}",
        r"\definecolor{bg}{rgb}{0.04,0.05,0.09}",

        # BLACK BACKGROUND for full image
        r"\fill[black] (-5.5,-5.5) rectangle (5.5,5.5);",

        # Include the graph with adaptive sizing
        r"\begin{scope}",
        r"\clip (-5.5,-5.5) rectangle (5.5,5.5);",
        r"\node at (0,0) {\includegraphics[width=" + f"{graph_width_cm:.1f}cm" + r"]{" + graph_file + r"}};",
        r"\end{scope}",

        # Add Dynkin diagram(s) with FIXED margins and rotation for duals
        dynkin_diagram_tikz(diagram_label, shift_x=diagram_x, shift_y=diagram_y, rotate=rotation),

        r"\end{tikzpicture}",
        r"\end{document}",
    ]

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / tex_file
    with open(tex_path, 'w') as f:
        f.write("\n".join(tikz))

    return tex_path, output_file, graph_width_cm


def main():
    """Regenerate all 60 Carter images with fixed margins."""

    # Load improved mapping
    with open('/Users/miriam/symmetric-graphs/e8/carter_to_orbital_improved.json') as f:
        mapping_data = json.load(f)

    mappings = mapping_data['mappings']

    # Load full 112 mapping to get hecke_family for class_* labels
    with open('/Users/miriam/symmetric-graphs/e8/carter_to_hecke_mapping_112.json') as f:
        full_data = json.load(f)
    full_mappings = full_data['mappings']

    print(f"Regenerating {len(mappings)} Carter class images with FIXED MARGINS...")
    print(f"Background: Black (full image)")
    print(f"Diagram margins: 0.6 units from top and left edges (larger for visibility)")
    print(f"Graph sizing: Adaptive to avoid overlap with diagrams")
    print()

    # Generate images
    generated = []
    failed = []
    skipped = []

    for i, (label, graph_file) in enumerate(sorted(mappings.items()), 1):
        # Check if graph exists
        graph_path = Path("/Users/miriam/symmetric-graphs/e8") / graph_file
        if not graph_path.exists():
            print(f"[{i:3d}/{len(mappings)}] {label:25s} ⊘ (graph not found)")
            skipped.append(label)
            continue

        try:
            # Get hecke_family for class_* labels
            hecke_family = None
            if label.startswith('class_') and label in full_mappings:
                hecke_family = full_mappings[label].get('hecke_family')

            tex_path, pdf_name, graph_size = create_carter_image_fixed_margins(label, graph_file, hecke_family)

            # Compile
            result = subprocess.run(
                ['pdflatex', '-interaction=nonstopmode', str(tex_path)],
                cwd='/Users/miriam/symmetric-graphs/e8',
                capture_output=True,
                timeout=30
            )

            if result.returncode == 0:
                size_note = f"(graph {graph_size:.1f}cm)" if graph_size < 10.0 else ""
                print(f"[{i:3d}/{len(mappings)}] {label:25s} ✓ {size_note}")
                generated.append((label, graph_file, pdf_name))
            else:
                print(f"[{i:3d}/{len(mappings)}] {label:25s} ✗ (compile failed)")
                failed.append(label)
        except Exception as e:
            print(f"[{i:3d}/{len(mappings)}] {label:25s} ✗ ({e})")
            failed.append(label)

    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"✓ Generated: {len(generated)}/{len(mappings)} images")
    print(f"  Black background: Full image for diagram visibility")
    print(f"  FIXED margins: 0.6 units from edges (larger for clarity)")
    print(f"  Adaptive graph sizing: Smaller when needed to avoid overlap")
    print()

    if skipped:
        print(f"⊘ Skipped: {len(skipped)} (missing graphs)")

    if failed:
        print(f"✗ Failed: {len(failed)}")

    # Clean up aux files
    import os
    for f in os.listdir('/Users/miriam/symmetric-graphs/e8'):
        if f.startswith('carter_class_') and (f.endswith('.aux') or f.endswith('.log')):
            os.remove(f'/Users/miriam/symmetric-graphs/e8/{f}')

    print()
    print("✓ All images now have consistent diagram positioning!")


if __name__ == "__main__":
    main()
