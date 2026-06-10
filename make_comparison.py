#!/usr/bin/env python3
"""
Generate a 2x2 comparison image:
  Top row:    dodecahedron-download.png  |  pappus-download.png
  Bottom row: auto-dodecahedron (clean)  |  auto-pappus (clean)

Strips titles, vertex labels, and legends from the auto-*.tex files,
recompiles to PDF, converts to PNG, and combines all four.
"""

import os
import re
import subprocess
import tempfile
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__))
INDIV = os.path.join(BASE, "individual_graphs")


def strip_tex(tex_path):
    """Remove title/subtitle nodes, vertex number labels, and legend from a .tex file."""
    with open(tex_path) as f:
        lines = f.readlines()

    out = []
    in_legend = False
    for line in lines:
        # Skip title and subtitle nodes
        if r'\node[title]' in line or r'\node[subtitle]' in line:
            continue
        # Skip legend block
        if r'\node[draw, rounded corners, fill=gray!5' in line:
            in_legend = True
            continue
        if in_legend:
            if line.strip() == '};':
                in_legend = False
            continue
        # Remove vertex number labels: {1} -> {} but keep the node
        line = re.sub(r'(\(v\d+\) at \([^)]+\)) \{\d+\}', r'\1 {}', line)
        # Also make vertices smaller since no labels
        line = line.replace('minimum size=6.5mm', 'minimum size=3.5mm')
        out.append(line)

    return ''.join(out)


def compile_tex(tex_content, name, workdir):
    """Compile TeX to PDF, return PDF path."""
    tex_path = os.path.join(workdir, f"{name}.tex")
    pdf_path = os.path.join(workdir, f"{name}.pdf")
    with open(tex_path, 'w') as f:
        f.write(tex_content)
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", f"{name}.tex"],
        capture_output=True, text=True, cwd=workdir
    )
    if result.returncode != 0:
        print(f"LaTeX failed for {name}:")
        print(result.stdout[-500:])
        raise RuntimeError(f"pdflatex failed for {name}")
    return pdf_path


def pdf_to_png(pdf_path, dpi=300):
    """Convert PDF to PNG using pdftoppm, return PNG path."""
    out_prefix = pdf_path.replace('.pdf', '')
    subprocess.run(
        ["pdftoppm", "-png", "-r", str(dpi), "-singlefile", pdf_path, out_prefix],
        check=True
    )
    return out_prefix + ".png"


def flatten_to_white(img):
    """Convert RGBA image to RGB with white background."""
    if img.mode == 'RGBA':
        bg = Image.new('RGB', img.size, 'white')
        bg.paste(img, mask=img.split()[3])
        return bg
    return img.convert('RGB')


def main():
    from PIL import ImageDraw, ImageFont

    # Downloaded reference images
    dodec_dl = os.path.join(BASE, "dodecahedron-download.png")
    pappus_dl = os.path.join(BASE, "pappus-download.png")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Generate clean versions of auto-dodecahedron and auto-pappus
        clean_images = {}
        for name in ["dodecahedron", "pappus"]:
            tex_path = os.path.join(INDIV, f"auto-{name}.tex")
            clean_tex = strip_tex(tex_path)
            pdf_path = compile_tex(clean_tex, f"clean-{name}", tmpdir)
            png_path = pdf_to_png(pdf_path)
            clean_images[name] = flatten_to_white(Image.open(png_path))
            print(f"  Generated clean {name}: {clean_images[name].size}")

        # Load downloaded images (flatten transparency to white)
        dodec_img = flatten_to_white(Image.open(dodec_dl))
        pappus_img = flatten_to_white(Image.open(pappus_dl))
        print(f"  Downloaded dodecahedron: {dodec_img.size}")
        print(f"  Downloaded pappus: {pappus_img.size}")

        # 2x2 grid: columns are graphs, rows are drawing styles
        # Col 1: Dodecahedron, Col 2: Pappus
        # Row 1: Standard, Row 2: Symmetry-aware
        padding = 50
        gap = 30
        label_height = 50
        col_label_width = 0  # will compute

        # Target: uniform cell size based on a common height
        target_cell_h = 600

        # Scale all 4 images to target_cell_h
        def scale_to_height(img, h):
            scale = h / img.height
            return img.resize((int(img.width * scale), h), Image.LANCZOS)

        cells = [
            [scale_to_height(dodec_img, target_cell_h),
             scale_to_height(pappus_img, target_cell_h)],
            [scale_to_height(clean_images["dodecahedron"], target_cell_h),
             scale_to_height(clean_images["pappus"], target_cell_h)],
        ]

        # Make columns equal width (max width in each column)
        col_widths = [max(cells[r][c].width for r in range(2)) for c in range(2)]
        total_w = padding * 2 + sum(col_widths) + gap
        total_h = padding * 2 + target_cell_h * 2 + gap + label_height * 2

        canvas = Image.new('RGB', (total_w, total_h), 'white')
        draw = ImageDraw.Draw(canvas)

        try:
            font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 32)
            font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 28)
        except (OSError, IOError):
            font = ImageFont.load_default()
            font_small = font

        # Column headers
        col_labels = ["Dodecahedron", "Pappus graph"]
        for ci in range(2):
            cx = padding + sum(col_widths[:ci]) + ci * gap + col_widths[ci] // 2
            bbox = draw.textbbox((0, 0), col_labels[ci], font=font)
            tw = bbox[2] - bbox[0]
            draw.text((cx - tw // 2, padding // 2), col_labels[ci], fill='black', font=font)

        # Row labels and images
        row_labels = ["Standard", "Symmetry-aware"]
        for ri in range(2):
            y_top = padding + label_height + ri * (target_cell_h + gap)
            # Row label (left-aligned above row? or let's skip row labels
            # and just put them as a subtitle)

            for ci in range(2):
                x_left = padding + sum(col_widths[:ci]) + ci * gap
                img = cells[ri][ci]
                # Center image in its column
                x_offset = (col_widths[ci] - img.width) // 2
                canvas.paste(img, (x_left + x_offset, y_top))

        # Row labels on the left side - actually let's put them centered below each row
        for ri in range(2):
            y_label = padding + label_height + ri * (target_cell_h + gap) + target_cell_h + 2
            bbox = draw.textbbox((0, 0), row_labels[ri], font=font_small)
            tw = bbox[2] - bbox[0]
            draw.text(((total_w - tw) // 2, y_label), row_labels[ri],
                      fill='gray', font=font_small)

        out_path = os.path.join(BASE, "comparison.png")
        canvas.save(out_path, dpi=(300, 300))
        print(f"\nSaved: {out_path}")
        print(f"Size: {canvas.size}")


if __name__ == "__main__":
    main()
