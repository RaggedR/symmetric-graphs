#!/usr/bin/env python3
"""Build catalogue.pdf: every drawing in individual_graphs/, one per page.

For each PDF in individual_graphs/ this script:
  1. renders a JPEG next to it (pdftoppm, same convention as auto-gh22.jpg)
  2. emits a catalogue page: title, the JPEG, a brief description, provenance

Output: catalogue.tex + catalogue.pdf in the repo root (beside gallery.pdf).
Descriptions are curated below, sourced from docs/CATALOGUE.md; any PDF on
disk that has no entry still gets a page, flagged "(no description yet)".

Usage: python3 src/make_catalogue.py
"""

import glob
import os
import subprocess
import sys

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SRC_DIR)
INDIV = os.path.join(ROOT, "individual_graphs")
JPEG_DPI = 200

# (section title, [(file stem, page title, description), ...])
# Foster numbering follows the standard census (F020A = Desargues).
SECTIONS = [
    ("Foster Census Classics", [
        ("auto-k4", r"$K_4$ (F004A)",
         r"The complete graph on 4 vertices — the tetrahedron skeleton and the "
         r"smallest cubic arc-transitive graph. $\mathrm{Aut} = S_4$ of order 24."),
        ("auto-k33", r"$K_{3,3}$ (F006A)",
         r"The complete bipartite graph on $3+3$ vertices (the `utility graph'), "
         r"Levi graph of the degenerate digon geometry. "
         r"$\mathrm{Aut} = (S_3 \times S_3) \rtimes C_2$ of order 72."),
        ("auto-cube", r"Cube $Q_3$ = GP(4,1) (F008A)",
         r"The 3-dimensional hypercube on 8 vertices, drawn as two concentric "
         r"4-rings. $\mathrm{Aut} = C_2 \times S_4$ of order 48."),
        ("auto-petersen", r"Petersen = GP(5,2) (F010A)",
         r"The Petersen graph on 10 vertices: Kneser graph $K(5,2)$, "
         r"$\mathrm{Aut} = S_5 \cong \mathrm{PGL}(2,5)$ of order 120. Drawn on two "
         r"$D_{10}$-orbit rings of 5 — the classic pentagon/pentagram picture."),
        ("auto-heawood", r"Heawood (F014A)",
         r"The Heawood graph on 14 vertices: Levi (point--line incidence) graph of "
         r"the Fano plane $\mathrm{PG}(2,2)$ and the $(3,6)$-cage. "
         r"$\mathrm{Aut} = \mathrm{PGL}(2,7)$ of order 336."),
        ("auto-moebiuskantor", r"M\"obius--Kantor = GP(8,3) (F016A)",
         r"The M\"obius--Kantor graph on 16 vertices, generalized Petersen "
         r"$\mathrm{GP}(8,3)$, drawn on two $C_8$-orbit rings. "
         r"$\mathrm{Aut} = \mathrm{GL}(2,3) \rtimes C_2$ of order 96."),
        ("auto-pappus", r"Pappus (F018A)",
         r"The Pappus graph on 18 vertices: Levi graph of the Pappus configuration "
         r"$(9_3, 9_3)$, tied to the affine plane $\mathrm{AG}(2,3)$. "
         r"$\mathrm{Aut}$ has order 216; drawn as three $S_3$-orbit rings of 6."),
        ("auto-desargues", r"Desargues = GP(10,3) (F020A)",
         r"The Desargues graph on 20 vertices: Levi graph of the Desargues "
         r"configuration $(10_3, 10_3)$. $\mathrm{Aut} = C_2 \times S_5$ of order "
         r"240; two $C_{10}$-orbit rings of 10."),
        ("auto-dodecahedron", r"Dodecahedron (F020B)",
         r"The skeleton of the regular dodecahedron, 20 vertices. "
         r"$\mathrm{Aut} = C_2 \times A_5$ of order 120; four $C_5$-orbit rings "
         r"of 5 give the Schlegel-like layout."),
        ("auto-nauru", r"Nauru = GP(12,5) (F024A)",
         r"The Nauru graph on 24 vertices, generalized Petersen "
         r"$\mathrm{GP}(12,5)$. $\mathrm{Aut} = S_4 \times S_3$ of order 144; two "
         r"$C_{12}$-orbit rings of 12."),
        ("auto-f26a", r"$F_{26A}$",
         r"The F26A graph on 26 vertices, $\mathrm{Aut} = C_{13} \rtimes C_6$ of "
         r"order 78 — the smallest Foster-census entry whose full automorphism "
         r"group is solvable of this metacyclic shape. Two rings of 13."),
        ("auto-coxeter", r"Coxeter (F028A)",
         r"The Coxeter graph on 28 vertices, a famous non-hamiltonian cubic graph. "
         r"$\mathrm{Aut} = \mathrm{PGL}(2,7)$ of order 336 (the same group as "
         r"Heawood); four $D_{14}$-orbit rings of 7."),
        ("auto-tuttecoxeter", r"Tutte--Coxeter (F030A)",
         r"The Tutte--Coxeter graph on 30 vertices: Levi graph of the generalized "
         r"quadrangle $\mathrm{GQ}(2,2)$ and the $(3,8)$-cage. "
         r"$\mathrm{Aut} \cong \mathrm{Aut}(S_6)$ of order 1440 — the home of "
         r"$S_6$'s outer automorphism. Three rings of 10."),
    ]),
    ("PSL(2,q) Family", [
        ("auto-psl2_5_c3", r"$\mathrm{PSL}(2,5)/C_3$ — double cover of Petersen",
         r"Coset graph of $\mathrm{PSL}(2,5) \cong A_5$ on the 20 cosets of a "
         r"$C_3$ point stabilizer: a double cover of the Petersen graph. First "
         r"rung of the $C_3 \subset S_3 \subset S_4$ tower of covers."),
        ("auto-psl2_7_c3", r"$\mathrm{PSL}(2,7)/C_3$ — 56-vertex coset graph",
         r"Cubic coset graph of $\mathrm{PSL}(2,7)$ (order 168) acting on the 56 "
         r"cosets of a $C_3$ stabilizer, from the PSL(2,q) census programme."),
        ("auto-psl2_13", r"$\mathrm{PSL}(2,13)/S_3$ coset graph",
         r"Cubic coset graph of $\mathrm{PSL}(2,13)$ with $S_3$ point stabilizer "
         r"(the Foster-182 line of the census). Part of the evidence that the "
         r"$S_3$-stabilizer family is Ramanujan."),
        ("auto-biggssmith", r"Biggs--Smith (F102A)",
         r"The Biggs--Smith graph on 102 vertices, "
         r"$\mathrm{Aut} = \mathrm{PSL}(2,17)$ with $S_4$ stabilizer — "
         r"distance-transitive, and in this programme's slogan ``a picture of "
         r"$SU(2)$ over $\mathbb{F}_{17}$'' (see docs/RESEARCH\_PROGRAM.md)."),
        ("auto-biggssmith2", r"Biggs--Smith — alternative layout",
         r"A second rendering of the Biggs--Smith graph with a different orbit "
         r"quotient choice, for comparison with the previous page."),
        ("auto-cr19", r"Cubic Ramanujan graph from $\mathrm{PSL}(2,19)$",
         r"A cubic Ramanujan graph constructed from $\mathrm{PSL}(2,19)$ in the "
         r"census programme; spectral bounds checked in src/spectral\_analysis.py."),
        ("auto-cr31", r"Cubic Ramanujan graph from $\mathrm{PSL}(2,31)$",
         r"A cubic Ramanujan graph constructed from $\mathrm{PSL}(2,31)$ — the "
         r"largest member of the rendered Ramanujan series."),
    ]),
    ("Incidence Geometries", [
        ("auto-pg23", r"Levi graph of $\mathrm{PG}(2,3)$",
         r"Point--line incidence graph of the projective plane of order 3: "
         r"$13+13 = 26$ vertices, valency 4, acted on by $\mathrm{PSL}(3,3)$. The "
         r"$q=3$ sibling of the Heawood graph."),
        ("auto-gh22", r"Generalized hexagon $\mathrm{GH}(2,2)$ on 63 points",
         r"Collinearity drawing of the split Cayley hexagon $\mathrm{GH}(2,2)$: "
         r"the 63 points of the $G_2(2)$ hexagon with its collinearity relation. "
         r"Built by src/build\_hexagon.py."),
        ("auto-tuttecage12", r"Tutte 12-cage (F126A)",
         r"The Tutte 12-cage on 126 vertices: incidence graph of "
         r"$\mathrm{GH}(2,2)$, $\mathrm{Aut} = G_2(2)$ of order 12096. Notably "
         r"\emph{not} vertex-transitive — points and lines form orbits of "
         r"$63+63$, since the hexagon is self-dual only for $q = 3^h$."),
        ("tutte12cage-clean", r"Tutte 12-cage — clean rendering",
         r"Hand-tuned presentation drawing of the Tutte 12-cage (labels, legend "
         r"and clutter stripped) for use in papers and slides."),
        ("auto-meinhold72", r"Meinhold-72 orbital graph",
         r"Orbital graph of $G_2(2)$ acting on the 72 cosets of "
         r"$\mathrm{PSL}(3,2)$. Its rank-8 Hecke algebra satisfies $B^2 = 13I$ — "
         r"the ``why 13'' story ($13 = 7 + 2\cdot 3$ via the 36+36 chiral Fano "
         r"structures); see the why13 paper in the hecke repo."),
        ("F056A", r"$F_{056A}$",
         r"The Foster-census graph F056A on 56 vertices, cubic arc-transitive."),
        ("F056A-minimal", r"$F_{056A}$ — minimal rendering",
         r"Hand-edited minimal drawing of F056A, decluttered for presentation."),
    ]),
    ("Voltage Graphs \\& Variants", [
        ("auto-petersen_voltage", r"Petersen — voltage graph construction",
         r"The Petersen graph presented as a voltage graph: the base quotient "
         r"with $\mathbb{Z}_5$ voltage labels whose derived graph reconstructs "
         r"the full graph."),
        ("auto-heawood_voltage", r"Heawood — voltage graph construction",
         r"Voltage-graph presentation of the Heawood graph: a small quotient with "
         r"cyclic voltage labels encoding the 14-vertex cover."),
        ("auto-mk_voltage", r"M\"obius--Kantor — voltage graph construction",
         r"Voltage-graph presentation of the M\"obius--Kantor graph "
         r"$\mathrm{GP}(8,3)$."),
        ("tuttecoxeter-clean", r"Tutte--Coxeter — clean rendering",
         r"Hand-tuned presentation drawing of the Tutte--Coxeter graph "
         r"(F030A, the $(3,8)$-cage)."),
        ("cage-pair", r"Cage pair",
         r"Hand-edited side-by-side of two cages from the generalized-polygon "
         r"series at $s=2$: the $(3,8)$-cage (Tutte--Coxeter, $\mathrm{GQ}(2,2)$) "
         r"and the $(3,12)$-cage (Tutte 12-cage, $\mathrm{GH}(2,2)$)."),
    ]),
    ("Klein Quartic \\& Quotients", [
        ("auto-klein-c7", r"Klein quartic graph — $C_7$ layout",
         r"Graph derived from the Klein quartic's $\mathrm{PSL}(2,7)$ cell "
         r"structure (genus 3, Hurwitz-maximal with 168 automorphisms), laid out "
         r"on $C_7$ orbit rings."),
        ("auto-klein-minimal", r"Klein quartic graph — minimal rendering",
         r"Decluttered drawing of the Klein quartic graph."),
        ("klein-z2-quotient", r"Klein graph — $\mathbb{Z}_2$ quotient",
         r"Hand-edited drawing of a $\mathbb{Z}_2$ quotient of the Klein quartic "
         r"graph, from the quotient-surfaces side-programme."),
        ("zhou-z2-quotient", r"Zhou $\mathbb{Z}_2$-quotient",
         r"The Zhou $\mathbb{Z}_2$-quotient construction from the clayworth "
         r"CSS-code programme (quotients of graph tilings yielding quantum "
         r"error-correcting codes)."),
        ("zhou-quotient-spectral", r"Zhou quotient — spectral layout",
         r"The same Zhou quotient drawn with a spectral (eigenvector) layout for "
         r"comparison with the orbit-ring method."),
    ]),
]

PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\setlength{\parindent}{0pt}
\begin{document}

\begin{titlepage}
\centering
{\Huge\bfseries Catalogue of Symmetric Graph Drawings\par}
\vspace{1cm}
{\large One page per drawing in \texttt{individual\_graphs/}\par}
\vspace{2cm}
\begin{minipage}{0.8\textwidth}
Each graph was laid out by the orbit-quotient pipeline
(\texttt{src/draw\_graph.py} + \texttt{src/symmetric-layout.g}): GAP computes
$\mathrm{Aut}(\Gamma)$ and an orbit decomposition under a chosen subgroup, and
each orbit is placed on a concentric ring with rotational symmetry. Pages
marked \emph{hand-edited} were tuned manually after generation. Group facts
follow \texttt{docs/CATALOGUE.md}; the research context is in
\texttt{docs/RESEARCH\_PROGRAM.md}.
\end{minipage}
\vfill
{\small Generated by \texttt{src/make\_catalogue.py} --- regenerate after adding drawings.\par}
\end{titlepage}

\tableofcontents
\newpage
"""

PAGE = r"""\section{%(title)s}
\begin{center}
\includegraphics[width=0.92\textwidth,height=0.66\textheight,keepaspectratio]{%(image)s}
\end{center}
\vfill
%(description)s

\medskip
{\small\color{gray}\texttt{individual\_graphs/%(stem_tex)s.pdf}}
\newpage
"""


def pdf_to_jpeg(stem):
    """Render individual_graphs/<stem>.pdf to <stem>.jpg beside it."""
    pdf = os.path.join(INDIV, stem + ".pdf")
    out_prefix = os.path.join(INDIV, stem)
    subprocess.run(
        ["pdftoppm", "-jpeg", "-r", str(JPEG_DPI), "-singlefile", pdf, out_prefix],
        check=True,
    )
    return out_prefix + ".jpg"


def tex_escape_stem(stem):
    return stem.replace("_", r"\_")


def main():
    on_disk = {os.path.splitext(os.path.basename(p))[0]
               for p in glob.glob(os.path.join(INDIV, "*.pdf"))}
    catalogued = {stem for _, entries in SECTIONS for stem, _, _ in entries}

    missing_pdf = catalogued - on_disk
    if missing_pdf:
        sys.exit(f"ERROR: catalogue entries with no PDF on disk: {sorted(missing_pdf)}")

    parts = [PREAMBLE]
    sections = [(title, list(entries)) for title, entries in SECTIONS]

    # Uncatalogued drawings still get a page, in a trailing section.
    uncatalogued = sorted(on_disk - catalogued)
    if uncatalogued:
        sections.append((
            "Uncatalogued",
            [(stem, tex_escape_stem(stem), r"\emph{(no description yet --- add an "
              r"entry in src/make\_catalogue.py)}") for stem in uncatalogued],
        ))
        print(f"WARNING: {len(uncatalogued)} drawing(s) lack descriptions: {uncatalogued}")

    n_pages = 0
    for sec_title, entries in sections:
        parts.append("\\addtocontents{toc}{\\protect\\vspace{0.5em}"
                     "\\protect\\noindent\\textbf{%s}\\protect\\par}\n" % sec_title)
        for stem, title, desc in entries:
            print(f"  rendering {stem}.jpg")
            pdf_to_jpeg(stem)
            parts.append(PAGE % {
                "title": title,
                "image": f"individual_graphs/{stem}.jpg",
                "description": desc,
                "stem_tex": tex_escape_stem(stem),
            })
            n_pages += 1

    parts.append("\\end{document}\n")

    tex_path = os.path.join(ROOT, "catalogue.tex")
    with open(tex_path, "w") as f:
        f.write("".join(parts))
    print(f"Wrote {tex_path} ({n_pages} graph pages)")

    for _ in range(2):  # twice for the table of contents
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "catalogue.tex"],
            cwd=ROOT, capture_output=True, text=True,
        )
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-15:])
        sys.exit(f"pdflatex failed:\n{tail}")

    for ext in ("aux", "log", "out", "toc"):
        try:
            os.remove(os.path.join(ROOT, f"catalogue.{ext}"))
        except FileNotFoundError:
            pass
    print(f"Wrote {os.path.join(ROOT, 'catalogue.pdf')}")


if __name__ == "__main__":
    main()
