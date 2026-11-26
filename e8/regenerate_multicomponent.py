"""
Regenerate Carter class images with many components that have positioning issues.
"""

import json
import subprocess
import re
from pathlib import Path
from dynkin_diagrams import dynkin_diagram_tikz


def create_carter_image_with_dynkin(label, graph_file):
    """Create a Carter class image with proper Dynkin diagrams."""

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

        # Add Dynkin diagram(s) with improved positioning
        dynkin_diagram_tikz(label),

        r"\end{tikzpicture}",
        r"\end{document}",
    ]

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / tex_file
    with open(tex_path, 'w') as f:
        f.write("\n".join(tikz))

    return tex_path, output_file


def main():
    """Regenerate multi-component Carter images with better positioning."""

    # Load improved mapping
    with open('/Users/miriam/symmetric-graphs/e8/carter_to_orbital_improved.json') as f:
        mapping_data = json.load(f)

    mappings = mapping_data['mappings']

    # Find classes with multiple components (3+) or specifically mentioned ones
    to_regenerate = []
    for label, graph_file in mappings.items():
        components = re.findall(r'[A-Z]\d+', label)
        if len(components) >= 3:  # 3+ components likely need better positioning
            to_regenerate.append((label, graph_file, len(components)))

    to_regenerate.sort(key=lambda x: -x[2])  # Sort by component count

    print(f"Regenerating {len(to_regenerate)} Carter class images with improved positioning...")
    print(f"These classes have 3+ components that may be too close to margins")
    print()

    # Generate images
    generated = []
    failed = []
    skipped = []

    for i, (label, graph_file, num_comp) in enumerate(to_regenerate, 1):
        # Check if graph exists
        graph_path = Path("/Users/miriam/symmetric-graphs/e8") / graph_file
        if not graph_path.exists():
            print(f"[{i:3d}/{len(to_regenerate)}] {label:25s} ⊘ (graph not found)")
            skipped.append(label)
            continue

        try:
            print(f"[{i:3d}/{len(to_regenerate)}] {label:25s} ({num_comp} comp) ", end="")
            tex_path, pdf_name = create_carter_image_with_dynkin(label, graph_file)

            # Compile
            result = subprocess.run(
                ['pdflatex', '-interaction=nonstopmode', str(tex_path)],
                cwd='/Users/miriam/symmetric-graphs/e8',
                capture_output=True,
                timeout=30
            )

            if result.returncode == 0:
                print(f"✓")
                generated.append((label, graph_file, pdf_name))
            else:
                print("✗ (compile failed)")
                failed.append(label)
        except Exception as e:
            print(f"✗ ({e})")
            failed.append(label)

    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"✓ Regenerated: {len(generated)}/{len(to_regenerate)} images")
    print(f"  Improved positioning for multi-component diagrams")
    print()

    if skipped:
        print(f"⊘ Skipped: {len(skipped)} (missing graphs)")

    if failed:
        print(f"✗ Failed: {len(failed)}")

    print()
    print("Scaling adjustments:")
    print("  6+ components: 60% scale, 1.3 spacing")
    print("  4-5 components: 75% scale, 1.5 spacing")
    print("  2-3 components: 100% scale, 2.0 spacing")
    print()

    # Clean up aux files
    import os
    for f in os.listdir('/Users/miriam/symmetric-graphs/e8'):
        if f.startswith('carter_class_') and (f.endswith('.aux') or f.endswith('.log')):
            os.remove(f'/Users/miriam/symmetric-graphs/e8/{f}')


if __name__ == "__main__":
    main()
