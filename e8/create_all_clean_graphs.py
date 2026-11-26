"""
Create clean (no caption) versions of all orbital/root graphs used in Carter class images.
Uses TikZ clipping to remove bottom caption area from source PDFs.
"""
import subprocess
from pathlib import Path

# All source graphs that need clean versions
graphs_to_clean = [
    "g2_roots.pdf",
    "f4_roots.pdf",
    "a2_roots.pdf",
    "a3_roots.pdf",
    "d4_graph.pdf",
    "d6_graph.pdf",
    "d7_classical.pdf",
    "e5_clebsch.pdf",
    "e6_schlafli.pdf",
    "e7_mandala.pdf",
]

def create_clean_version(source_pdf):
    """Create a clean version by clipping bottom caption area."""

    base_name = source_pdf.replace(".pdf", "")
    clean_name = f"clean_{base_name}"

    # Create TikZ wrapper that includes the PDF but clips the caption
    tikz = [
        r"\documentclass[tikz,border=0pt]{standalone}",
        r"\usepackage{graphicx}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{tikzpicture}",

        # Clip to remove bottom caption area (adjust clip height as needed)
        r"\clip (-8.5,-7.5) rectangle (8.5,8.5);",

        # Include the source PDF
        r"\node at (0,0) {\includegraphics[width=16cm]{" + source_pdf + r"}};",

        r"\end{tikzpicture}",
        r"\end{document}",
    ]

    # Write tex file
    tex_path = Path("/Users/miriam/symmetric-graphs/e8") / f"{clean_name}.tex"
    with open(tex_path, 'w') as f:
        f.write("\n".join(tikz))

    return tex_path, f"{clean_name}.pdf"


print(f"Creating clean versions of {len(graphs_to_clean)} graphs...\n")

created = []
failed = []

for source in graphs_to_clean:
    print(f"Processing {source}...", end=" ")

    try:
        tex_path, pdf_name = create_clean_version(source)

        # Compile
        result = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', str(tex_path)],
            cwd='/Users/miriam/symmetric-graphs/e8',
            capture_output=True,
            timeout=30
        )

        if result.returncode == 0:
            print("✓")
            created.append(pdf_name)
        else:
            print("✗ (compile failed)")
            failed.append(source)

    except Exception as e:
        print(f"✗ ({e})")
        failed.append(source)

print(f"\n✓ Created {len(created)}/{len(graphs_to_clean)} clean graphs")

if failed:
    print(f"✗ Failed: {', '.join(failed)}")

# Save list
with open("/Users/miriam/symmetric-graphs/e8/clean_graphs_list.txt", "w") as f:
    for name in created:
        f.write(f"{name}\n")

print(f"\n✓ Saved list to clean_graphs_list.txt")
