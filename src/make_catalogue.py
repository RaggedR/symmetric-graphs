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
        ("zhou-z2-quotient", r"Zhou $\mathbb{Z}_2$ quotient — 91 vertices",
         r"The 2:1 quotient $\mathrm{Sab}(\mathrm{PSL}(2,13), D_{12})$ of the "
         r"previous page's F182A: fibres are the size-2 blocks of "
         r"$\mathrm{PSL}(2,13)$, fusing the cover's 14 $Z_{13}$-rings into 7 "
         r"pairs. Vertices inherit the cover's colours (fill = outer fibre "
         r"ring, border = inner), so this page and its cover read as a pair."),
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
        ("meinhold72-rings", r"Meinhold-72 orbital graph",
         r"Orbital graph of $G_2(2)$ acting on the 72 cosets of "
         r"$\mathrm{PSL}(3,2)$, drawn by a $C_7 < G_2(2)$ ring decomposition: "
         r"10 rings of 7 (plus the 2 fixed points at the centre), whose "
         r"offset-optimised layout aligns into seven radial spokes. Its rank-8 "
         r"Hecke algebra satisfies $B^2 = 13I$ — the ``why 13'' story "
         r"($13 = 7 + 2\cdot 3$ via the 36+36 chiral Fano structures); see the "
         r"why13 paper in the hecke repo."),
        ("dmt", r"Twisted triality hexagon $\mathrm{GH}(2,8)$ — ${}^3D_4(2)$",
         r"Collinearity graph of the twisted triality hexagon "
         r"$\mathrm{GH}(2,8)$: 819 points, valency 18, distance-regular of "
         r"diameter 3, acted on by ${}^3D_4(2)$ of order 211{,}341{,}312 — the "
         r"triality twist of $D_4$, whose $q$ vs $q^3$ parameter asymmetry is "
         r"the scar of the threefold Dynkin symmetry. Laid out as 117 "
         r"concentric $C_7$-rings by surfaces/draw\_triality\_hexagon.py; "
         r"essay in surfaces/triality\_hexagon.md."),
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
    ]),
]

# Hecke-algebra stories, keyed by file stem.  Appended to the page
# description as a "Hecke." paragraph.  All facts computed in the hecke repo
# (~/git/hecke/): hecke_spherical.g/.out, hecke_irreps*.g/.out,
# hecke_catalogue_gaps.g/.out; written up in hecke-irreps.pdf (and why13.pdf
# for Meinhold-72).  T_i = orbital basis sorted by suborbit length.
HECKE = {
    # ── Foster census classics ──
    "auto-k4": r"Every Hecke story in this catalogue starts here. $S_4$ on "
        r"4 points is 2-transitive, so the algebra is the smallest "
        r"possible: rank 2, $\mathcal H\cong\C^2$, spherical functions "
        r"$(1,3)$ and $(1,-1)$. Rank 2 means the pair distinguishes only "
        r"``equal'' and ``adjacent'' — there is no geometry yet.",
    "auto-k33": r"One step up: rank 3, $\mathcal H\cong\C^3$. The new row "
        r"is the bipartite sign character $(1,2,-3)$ — the algebra sees "
        r"the two sides of the bipartition before any drawing does. The "
        r"remaining rows are $(1,2,3)$ and $(1,-1,0)$.",
    "auto-cube": r"Rank 4, $\mathcal H\cong\C^4$, spherical table fully "
        r"rational: $(1,\pm1,\pm3,3)$ and $(1,\pm1,\mp1,-1)$. The "
        r"antipodal map supplies the extra linear character: the cube is "
        r"the first graph here whose algebra detects ``distance 3 = "
        r"opposite vertex''.",
    "auto-petersen": r"The Petersen graph is the Gelfand pair "
        r"$(A_5,S_3)$ wearing a drawing: rank 3, spherical rows "
        r"$(1,-2,1)_{4}$ and $(1,1,-2)_{5}$ on suborbits $(1,3,6)$, "
        r"reassembling the spectrum $3,1^{(5)},(-2)^{(4)}$. The identical "
        r"table is realized by the coset presentation "
        r"$\mathrm{PSL}(2,5)/S_3$.",
    "auto-heawood": r"The first irrational number in the catalogue. "
        r"Rank 4: two rational rows from the bipartition, and a degree-6 "
        r"pair $(1,\mp\sqrt2,\pm\sqrt2,-1)$. That $\pm\sqrt2$ is "
        r"$\pm\sqrt q$ at $q=2$ — the arithmetic fingerprint of the Fano "
        r"plane $\mathrm{PG}(2,2)$ whose incidence graph this is.",
    "auto-moebiuskantor": r"Rank 6; the degree-4 rows carry $\pm\sqrt3$ on "
        r"the two length-3 suborbits and everything else is rational. "
        r"$\mathrm{GP}(8,3)$ trades Heawood's Fano $\sqrt2$ for $\sqrt3$.",
    "auto-pappus": r"Rank 5, with degree-6 rows "
        r"$(1,-1,\pm\sqrt3,0,\mp\sqrt3)$: the $\sqrt3$ of the affine plane "
        r"$\mathrm{AG}(2,3)$ — the affine sibling of Heawood's projective "
        r"$\sqrt2$, one prime up.",
    "auto-desargues": r"Rank 6, $\mathcal H\cong\C^6$, all rational. The "
        r"constituent degrees $1,1,4,4,5,5$ are Petersen's $1,4,5$ "
        r"doubled by the bipartite sign — exactly what a bipartite double "
        r"cover should do, and a reminder that doubling can be "
        r"arithmetically tame.",
    "auto-dodecahedron": r"Under the full $C_2\times A_5$ the pair is "
        r"Gelfand of rank 6, and the famous golden-ratio $\pm\sqrt5$ sits "
        r"in the degree-3 spherical rows. But restrict to "
        r"$\mathrm{PSL}(2,5)$ with $C_3$ stabilizer and the \emph{same "
        r"graph} refines to a rank-8 non-Gelfand pair (next section): the "
        r"first hint that Hecke algebras belong to pairs, not graphs. "
        r"\emph{Non-trivial representations: none for this pair (Gelfand, "
        r"all irreps $1$-dimensional); the refined pair's $2\times2$ "
        r"matrices are on the $\mathrm{PSL}(2,5)/C_3$ page.}",
    "auto-nauru": r"Rank 7, $\mathcal H\cong\C^7$, all-rational table "
        r"(degrees $1,1,3,3,4,6,6$) — rich in suborbits, arithmetically "
        r"flat. The largest commutative algebra among the small Foster "
        r"graphs.",
    "auto-f26a": r"The first rebel. F26A is the only Foster graph "
        r"$\le 30$ whose Hecke algebra is noncommutative: "
        r"$\mathcal H\cong\C^2\oplus M_2\oplus M_2$, rank 10, because "
        r"$\mathrm{Aut}=C_{13}{:}C_6$ of order $78=3\cdot26$ is "
        r"\emph{sharply} arc-transitive — no room for multiplicity-"
        r"freeness. The two $M_2$ blocks are Galois conjugate over "
        r"$\Q(\sqrt{13})$, and their traceless valency-3 generators "
        r"satisfy $T^2=\frac{5\mp\sqrt{13}}{2}I$ with $N_+N_-=3$, the "
        r"valency. The $\sqrt{13}$ comes from Gauss periods of "
        r"$\Q(\zeta_{13})$ — the same numeral as Meinhold-72, through a "
        r"different door. Non-trivial representations: two. With "
        r"$\eta_0,\eta_1=\frac{-1\pm\sqrt{13}}2$: "
        r"the first block is generated by $\rho(T_2)=\bigl(\begin{smallmatrix}"
        r"1&0\\0&-1\end{smallmatrix}\bigr)$ and $\rho(T_3)=\bigl("
        r"\begin{smallmatrix}\eta_0/2&1\\-2\eta_0-\frac54\eta_1&-\eta_0/2"
        r"\end{smallmatrix}\bigr)$, with relations $\rho(T_2)^2=I$, "
        r"$\rho(T_3)^2=\frac{5-\sqrt{13}}2I$, $\rho(T_7)^2=\frac{5+\sqrt{13}}2I$; "
        r"the second block is its Galois conjugate ($\eta_0\leftrightarrow\eta_1$).",
    "auto-coxeter": r"Gelfand of rank 5 via $\mathrm{PGL}(2,7)/D_{12}$ "
        r"(suborbits $1,3,6,6,12$). The valency-3 orbital takes spherical "
        r"values $3,\ -1\pm\sqrt2,\ -1,\ 2$, reassembling the Coxeter "
        r"spectrum $\{3,\ 2^{8},\ (-1)^{7},\ (-1\pm\sqrt2)^{6}\}$ — "
        r"$\sqrt2$ again, but shifted by $-1$: same field as Heawood, "
        r"different geometry.",
    "auto-tuttecoxeter": r"The $(3,8)$-cage is the incidence graph of the "
        r"generalized quadrangle $\mathrm{GQ}(2,2)$: rank 5 on suborbits "
        r"$(1,3,6,8,12)$, spectrum $\{\pm3,\ \pm2^{(9)},\ 0^{(10)}\}$. "
        r"Here $\pm2=\pm\sqrt{st}$ is rational — the quadrangle keeps its "
        r"surds hidden at $s=t=2$.",
    # ── PSL(2,q) family ──
    "auto-psl2_5_c3": r"The first non-Gelfand rung of the tower, and a "
        r"20-point miniature of Meinhold-72: $\mathcal H\cong\C^4\oplus "
        r"M_2$, rank 8, same shape one level down. In the $M_2$ block "
        r"$A^2=I$, $C^2+C+3I=0$ (determinant $=3=$ valency) and "
        r"$(\rho(T_6)-I)^2=4I$ — so the role of Meinhold's 13 is played "
        r"here by 4, and there is \emph{no} universal formula: the SRG "
        r"mechanism behind 13 is special. The dodecahedral $\pm\sqrt5$ "
        r"survives, but exiled to the 1-dimensional irreps. "
        r"The non-trivial irrep (dimension 2, on $\chi_4$, degree 4) is "
        r"generated by $A=\rho(T_2)=\bigl(\begin{smallmatrix}1&0\\0&-1"
        r"\end{smallmatrix}\bigr)$ and $C=\rho(T_5)=\bigl(\begin{smallmatrix}"
        r"1/2&1\\-15/4&-3/2\end{smallmatrix}\bigr)$, with relations $A^2=I$, "
        r"$C^2+C+3I=0$; then $\rho(T_6)=CA$, $\rho(T_8)=AC$, and "
        r"$\rho(T_3)=\rho(T_4)=\bigl(\begin{smallmatrix}-2&0\\0&0"
        r"\end{smallmatrix}\bigr)$.",
    "auto-psl2_7_c3": r"$\mathcal H\cong\C^3\oplus M_2\oplus M_2\oplus "
        r"M_3$, rank 20. Both $M_2$ blocks satisfy $B^2=2I$ where "
        r"$2=\mathrm{Nm}\bigl(\frac{-1+\sqrt{-7}}{2}\bigr)$ is the prime "
        r"that splits in $\Q(\sqrt{-7})$ — the field of the ambient "
        r"$\mathrm{PSL}(2,7)$. One Hecke generator has finite order 3, "
        r"and the $M_3$ block is integral with genuinely cubic "
        r"eigenvalues ($x^3+x-6$): the first block in the catalogue that "
        r"no quadratic surd can describe. Non-trivial irreps, explicitly: "
        r"first $M_2$ (with $\alpha=\frac{-1+\sqrt{-7}}2$, so "
        r"$\alpha\bar\alpha=2$): $D=\rho(T_3)=\bigl(\begin{smallmatrix}"
        r"\alpha&0\\0&\bar\alpha\end{smallmatrix}\bigr)$, "
        r"$S=\rho(T_5)=\bigl(\begin{smallmatrix}0&\bar\alpha\\\alpha&0"
        r"\end{smallmatrix}\bigr)$; relations $D^2+D+2I=0$, $S^2=2I$, "
        r"$SDS^{-1}=\bar D$. Second $M_2$: $B=\rho(T_5)=\bigl("
        r"\begin{smallmatrix}-1/2&1\\7/4&1/2\end{smallmatrix}\bigr)$ with "
        r"$B^2=2I$, and $R=\rho(T_{10})=\bigl(\begin{smallmatrix}1/2&1\\"
        r"-7/4&-3/2\end{smallmatrix}\bigr)$ of order 3 ($R^2+R+I=0$). "
        r"$M_3$ block: $\rho(T_2)=\mathrm{diag}(1,-1,-1)$ and "
        r"$\rho(T_5)=\bigl(\begin{smallmatrix}1&1&0\\-2&0&-2\\-4&0&-1"
        r"\end{smallmatrix}\bigr)$ (charpoly $x^3+x-6$) generate; "
        r"$\rho(T_{11})$ has eigenvalues $0,\frac{-1\pm\sqrt{17}}2$.",
    "auto-psl2_13": r"The biggest algebra fully worked out: "
        r"$\mathcal H\cong\C^4\oplus M_2^{\oplus3}\oplus M_3^{\oplus2}$, "
        r"rank 34. The spherical values live in $\Q(\sqrt{13})$ via the "
        r"\emph{F26A} Gauss periods (a different 13 than Meinhold's!); "
        r"the three $M_2$ blocks are Galois conjugate over "
        r"$\Q(\zeta_7)^+$; the integral Steinberg $M_3$ block has "
        r"$\rho(T_3)^3=-I$ and charpoly $x^3+x-26$ (echoing "
        r"$\mathrm{PSL}(2,7)$'s $x^3+x-6$). The punchline is the "
        r"degree-14 block: $\rho(T_3)^3=5I$ — the first \emph{cubic} "
        r"analogue of Meinhold's $B^2=13I$. Non-trivial irreps, "
        r"explicitly (the $M_2$ blocks have Gauss-period entries, "
        r"conjugated by $c_k\mapsto c_{2k}$; the $M_3$ blocks are clean): "
        r"Steinberg $M_3$ (deg 13): $\rho(T_2)=\mathrm{diag}(1,1,-1)$, "
        r"$\rho(T_3)=\bigl(\begin{smallmatrix}2&-1&-2\\2&-2&-3\\1&0&0"
        r"\end{smallmatrix}\bigr)$ with $\rho(T_3)^3=-I$, and "
        r"$\rho(T_{21})=\bigl(\begin{smallmatrix}-2&1&-3\\-5&4&-1\\3&-1&-2"
        r"\end{smallmatrix}\bigr)$ with charpoly $x^3+x-26$. Degree-14 "
        r"$M_3$: $\rho(T_3)=\bigl(\begin{smallmatrix}1/2&1/2&-11/4\\2&1&3\\"
        r"1&0&-3/2\end{smallmatrix}\bigr)$, charpoly $x^3-5$, so "
        r"$\rho(T_3)^3=5I$ on the nose.",
    "auto-biggssmith": r"Gelfand of rank 8 — diameter 7, the deepest "
        r"distance-transitive graph in the catalogue (suborbits "
        r"$1,3,6,8,12,24,24,24$). Its spherical table is a tour of small "
        r"number fields: a $\Q(\sqrt{17})$ pair (degree 9, values "
        r"$\frac{1\pm\sqrt{17}}{2}$ on the valency-3 orbital), a cubic "
        r"triple over $\Q(\zeta_9)^+$ (degree 16, $\beta^3-3\beta+1=0$), "
        r"and rational rows (degrees 17, 18) — together reassembling the "
        r"spectrum $\{3,\ \frac{1\pm\sqrt{17}}{2}^{(9)},\ "
        r"(t^3+3t^2-3)^{(16\times3)},\ 0^{(17)},\ 2^{(18)}\}$.",
    "auto-biggssmith2": r"Same graph, same pair, same algebra — see the "
        r"previous page. The Hecke algebra is blind to layout; only the "
        r"drawing changed. \emph{Non-trivial representations: none "
        r"(Gelfand).}",
    "auto-cr19": r"The F26A mechanism, one prime up. "
        r"$\mathrm{Aut}=C_{19}{:}C_6$ is again sharply arc-transitive "
        r"with $C_3$ vertex stabilizer, and the pair has rank 14, "
        r"$\mathcal H\cong\C^2\oplus M_2^{\oplus3}$: one $2{\times}2$ "
        r"block for each of the $(19-1)/6=3$ orbits of sextic Gauss "
        r"periods, permuted by the Galois group of the cubic period "
        r"field inside $\Q(\zeta_{19})$. The generator charpolys have "
        r"period coefficients — the same arithmetic that makes the graph "
        r"Ramanujan makes its Hecke algebra wild. Non-trivial "
        r"representations: three. In each block "
        r"$\rho(T_2)=\bigl(\begin{smallmatrix}1&0\\0&-1\end{smallmatrix}"
        r"\bigr)$ and the valency-3 generators take the shape "
        r"$\rho(T_3)=\bigl(\begin{smallmatrix}\eta_a/2&1\\ *&\eta_a/2"
        r"\end{smallmatrix}\bigr)$, where $\eta_a=\sum_{s\in H}\zeta^{as}$ "
        r"($H=\langle8\rangle\le(\mathbb Z/19)^\times$) is a sextic Gauss "
        r"period; the partner $\rho(T_4)$ is traceless with "
        r"$\rho(T_4)^2=N_aI$ for an algebraic integer $N_a$ of the cubic "
        r"period field, and $\mathrm{Gal}$ permutes the three blocks by "
        r"$a\mapsto 2a$.",
    "auto-cr31": r"The largest member of the Frobenius family: "
        r"$\mathrm{Aut}=C_{31}{:}C_6$, rank 22, "
        r"$\mathcal H\cong\C^2\oplus M_2^{\oplus5}$ — five Galois-"
        r"conjugate $2{\times}2$ blocks, one per orbit of sextic Gauss "
        r"periods, conjugate over the quintic period field of conductor "
        r"31. F26A $\to$ CR(19) $\to$ CR(31) is one story told over "
        r"deeper and deeper subfields of $\Q(\zeta_p)$: the block count "
        r"is $(p-1)/6$. The five non-trivial irreps have the same period "
        r"shape as CR(19)'s: $\rho(T_2)=\bigl(\begin{smallmatrix}1&0\\"
        r"0&-1\end{smallmatrix}\bigr)$, valency-3 generators with sextic-"
        r"period diagonal entries $\eta_a/2$ "
        r"($H=\langle26\rangle\le(\mathbb Z/31)^\times$), and traceless "
        r"partners squaring to scalars $N_aI$ in the quintic period "
        r"field; $a\mapsto2a$ cycles the five blocks.",
    # ── Incidence geometries ──
    "auto-pg23": r"The point action $\mathrm{PSL}(3,3)/P_1$ on 13 points "
        r"is 2-transitive, so like $K_4$ it tells the shortest story: "
        r"rank 2, spherical functions $(1,12)$ and $(1,-1)$. The real "
        r"algebra of a projective plane lives on its flags — see the "
        r"Heawood page for $q=2$.",
    "auto-gh22": r"Home ground of the why-13 question. The point action "
        r"$G_2(2)/P_1$ is Gelfand of rank 4 (suborbits $1,6,24,32$), "
        r"spherical rows $(1,-3,6,-4)_{14}$, $(1,3,0,-4)_{21}$, "
        r"$(1,-1,-4,4)_{27}$ — the distance-regular eigenvalues "
        r"$6,3,-1,-4$ of the hexagon. That $-4$ returns as the $-4$ in "
        r"Meinhold's $AB+BA=-4I$ two pages on.",
    "auto-tuttecage12": r"A graph that is two coset spaces in a trench "
        r"coat. $\mathrm{Aut}=G_2(2)$ has \emph{two} orbits of 63 "
        r"(points and lines of $\mathrm{GH}(2,2)$ — the cage is not "
        r"vertex-transitive), and the two orbits give \emph{literally "
        r"identical} spherical tables, both equal to the $G_2(2)/P_1$ "
        r"table — even though no automorphism exchanges them, since the "
        r"hexagon is self-dual only for $q=3^h$. \emph{Non-trivial "
        r"representations: none on either orbit — both pairs are Gelfand, "
        r"all irreps $1$-dimensional.}",
    "tutte12cage-clean": r"Same 126 vertices as the previous page — and "
        r"stripped of labels the graph still refuses to be vertex-"
        r"transitive. The twin spherical tables are on the F126A page. "
        r"\emph{Non-trivial representations: none (both orbit pairs "
        r"Gelfand).}",
    "meinhold72-rings": r"The headline act. Rank 8, $\mathcal H\cong\C^2"
        r"\oplus M_2\oplus\C^2$, and in the $M_2$ block three relations "
        r"hold: $A^2=I$, $B^2=13I$, $AB+BA=-4I$, with the directed "
        r"orbital pair $BA$, $AB$ of charpoly $x^2+4x+13$ — the Gaussian "
        r"integers $-2\pm3i$ of norm 13. Why 13? Because "
        r"$13=7+2\cdot3$: $NN^{\mathsf T}=7I+2A$ for the "
        r"$\mathrm{SRG}(36,21,12,12)$ on the 36 chiral Fano structures, "
        r"and $3=\sqrt{k-\mu}$. The $-4$ is the $U_3(3)$-graph "
        r"eigenvalue, inherited from the hexagon two pages back. "
        r"Explicitly, the non-trivial irrep (dimension 2, on the two "
        r"degree-6 constituents) is generated by $A=\rho(T_2)=\bigl("
        r"\begin{smallmatrix}1&0\\0&-1\end{smallmatrix}\bigr)$ and "
        r"$B=\rho(T_4)=\bigl(\begin{smallmatrix}-2&-1\\-9&2"
        r"\end{smallmatrix}\bigr)$; then $\rho(T_3)=BA=\bigl("
        r"\begin{smallmatrix}-2&1\\-9&-2\end{smallmatrix}\bigr)$, "
        r"$\rho(T_5)=AB$, $\rho(T_6)=\bigl(\begin{smallmatrix}-2&1\\9&2"
        r"\end{smallmatrix}\bigr)$, $\rho(T_7)=3A$, $\rho(T_8)=3I$. Full "
        r"story: why13.pdf; independent re-derivation: hecke-irreps.pdf.",
    "dmt": r"One level up the hexagonal tower. ${}^3D_4(2)$ acting on "
        r"the 819 points of $\mathrm{GH}(2,8)$ is Gelfand of rank 4, "
        r"$\mathcal H\cong\C^4$ with constituents $1\oplus26\oplus324"
        r"\oplus468$, spectrum $\{18,\ 5^{(324)},\ (-3)^{(468)},\ "
        r"(-9)^{(26)}\}$. Because $st=16$ is a perfect square the "
        r"hexagon surds $s-1\pm\sqrt{st}=5,-3$ are rational, and the "
        r"lonely $-9=-(t+1)$ pairs with the 26-dimensional irrep. Same "
        r"$\C^4$ shape as $\mathrm{GH}(2,2)$; whether a Meinhold-style "
        r"72-fold doubling exists at $q=8$ is open.",
    "F056A": r"This is the $\mathrm{PSL}(2,7)/C_3$ coset graph — see "
        r"that page for its Hecke algebra ($\C^3\oplus M_2\oplus "
        r"M_2\oplus M_3$, $B^2=2I$). \emph{Non-trivial representations: "
        r"three (two $2$-dimensional, one $3$-dimensional); explicit "
        r"matrices and relations on the $\mathrm{PSL}(2,7)/C_3$ page.}",
    "F056A-minimal": r"Same graph as the previous page; the algebra "
        r"($\C^3\oplus M_2\oplus M_2\oplus M_3$) is on the "
        r"$\mathrm{PSL}(2,7)/C_3$ page. \emph{Non-trivial "
        r"representations: three — matrices there.}",
    # ── Voltage graphs & variants ──
    "auto-petersen_voltage": r"Covering theory is the Hecke story in "
        r"reverse: the $\mathbb{Z}_5$ voltage quotient compresses "
        r"Petersen to a 2-vertex base, and the spectrum of the cover "
        r"decomposes along the 5 characters of the voltage group. The "
        r"pair and its rank-3 algebra are on the Petersen page. "
        r"\emph{Non-trivial representations: none (Gelfand).}",
    "auto-heawood_voltage": r"The $\pm\sqrt2$ of the Heawood spherical "
        r"table (see the Heawood page) is invisible in the voltage base "
        r"graph — it only appears when the cyclic cover is unrolled. "
        r"\emph{Non-trivial representations: none (Gelfand).}",
    "auto-mk_voltage": r"Voltage presentation of $\mathrm{GP}(8,3)$; the "
        r"rank-6 algebra with its $\pm\sqrt3$ rows is on the "
        r"M\"obius--Kantor page. \emph{Non-trivial representations: none "
        r"(Gelfand).}",
    "tuttecoxeter-clean": r"Same $(3,8)$-cage; rank-5 algebra and "
        r"spectrum $\{\pm3,\pm2^{(9)},0^{(10)}\}$ on the Tutte--Coxeter "
        r"page. \emph{Non-trivial representations: none (Gelfand).}",
    "cage-pair": r"Two cages, two algebras. The $(3,8)$-cage carries the "
        r"rank-5 $\C^5$ of $\mathrm{GQ}(2,2)$ on one vertex-transitive "
        r"orbit; the $(3,12)$-cage carries two \emph{identical} copies "
        r"of the rank-4 $\C^4$ of $\mathrm{GH}(2,2)$, one per 63-orbit, "
        r"because it is not vertex-transitive. The Hecke tables can tell "
        r"the cages apart before you count girth. \emph{Non-trivial "
        r"representations: none on either cage (all pairs Gelfand).}",
    # ── Klein quartic & quotients ──
    "auto-klein-c7": r"The cubic Klein graph — the 1-skeleton of the "
        r"$\{7,3\}$ tiling of the genus-3 Klein quartic, on the 56 "
        r"cosets of $C_3$ in $\mathrm{PSL}(2,7)$ but along a girth-7 "
        r"orbital (F056A, girth 6, uses a different one). Under the full "
        r"automorphism group $\mathrm{PGL}(2,7)$ (order 336, stabilizer "
        r"$S_3$) the pair has rank 13, "
        r"$\mathcal H\cong\C^5\oplus M_2\oplus M_2$: the 1-dimensional "
        r"rows carry $\pm\sqrt2$, one block contains an orbital of "
        r"finite order 3 ($x^2+x+1$), and the other has a traceless "
        r"valency-6 generator with $T^2=8I$ — Hurwitz-maximal symmetry, "
        r"yet wilder than every Foster graph except F26A. The two "
        r"non-trivial irreps, explicitly: the first is integral, "
        r"generated by $\rho(T_2)=\bigl(\begin{smallmatrix}1&0\\0&-1"
        r"\end{smallmatrix}\bigr)$, $\rho(T_3)=\bigl(\begin{smallmatrix}"
        r"-1&1\\-4&0\end{smallmatrix}\bigr)$ ($x^2+x+4$) and "
        r"$\rho(T_4)=\bigl(\begin{smallmatrix}-1&-1\\-4&0\end{smallmatrix}"
        r"\bigr)$ ($x^2+x-4$, eigenvalues $\frac{-1\pm\sqrt{17}}2$). In "
        r"the second, $R=\rho(T_3)=\bigl(\begin{smallmatrix}1/2&1\\-7/4&"
        r"-3/2\end{smallmatrix}\bigr)$ has order 3 — \emph{literally the "
        r"same matrix} as the order-3 generator in the F056A pair's "
        r"second block — and $\rho(T_{11})=\bigl(\begin{smallmatrix}-1&2\\"
        r"7/2&1\end{smallmatrix}\bigr)$ satisfies $\rho(T_{11})^2=8I$.",
    "auto-klein-minimal": r"Same cubic Klein graph; the rank-13 algebra "
        r"$\C^5\oplus M_2\oplus M_2$ is on the $C_7$-layout page. "
        r"\emph{Non-trivial representations: two $2$-dimensional — "
        r"explicit matrices there.}",
    "klein-z2-quotient": r"Quotients break the symmetry that Hecke "
        r"algebras feed on: a $\mathbb{Z}_2$ quotient of the Klein graph "
        r"is no longer vertex-transitive, so there is no single pair "
        r"$(\mathrm{Aut},\mathrm{Stab})$ — the algebra lives upstairs on "
        r"the cover. \emph{Non-trivial representations: not applicable — "
        r"no transitive pair, hence no Hecke algebra downstairs.}",
    "zhou-z2-quotient": r"Like all proper quotients, the Zhou "
        r"$\mathbb{Z}_2$-quotient trades transitivity for topology (here "
        r"in service of CSS codes), so its Hecke story is its cover's. "
        r"What survives downstairs is the voltage-invariant half of the "
        r"spectrum. \emph{Non-trivial representations: not applicable — "
        r"no transitive pair downstairs.}",
    "zhou-quotient-spectral": r"Same Zhou quotient, drawn by its "
        r"eigenvectors — a spectral layout is itself a small Hecke joke: "
        r"the coordinates are spherical functions of whatever symmetry "
        r"remains. \emph{Non-trivial representations: not applicable — "
        r"no transitive pair.}",
}

# Gelfand pages: say explicitly that there is no higher-dimensional irrep.
_GELFAND_NOTE = (
    r" \emph{Non-trivial representations: none — the pair is Gelfand, so "
    r"every irreducible representation of $\mathcal H$ is $1$-dimensional "
    r"(the spherical rows are the complete list).}")
for _stem in ["auto-k4", "auto-k33", "auto-cube", "auto-petersen",
              "auto-heawood", "auto-moebiuskantor", "auto-pappus",
              "auto-desargues", "auto-nauru", "auto-coxeter",
              "auto-tuttecoxeter", "auto-biggssmith", "auto-pg23",
              "auto-gh22", "dmt"]:
    HECKE[_stem] += _GELFAND_NOTE

# Pages whose graph was redrawn in the dmt ring-palette style
# (clayworth/generate_named_graphs.py et al., 2026-06-10) use the new image
# from lean/named_graphs/ instead of the individual_graphs render.
# Paths are repo-root-relative.
IMAGE_OVERRIDES = {
    "auto-k4": "clayworth/named_graphs/k4-F004A.jpg",
    "auto-k33": "clayworth/named_graphs/k33-F006A.jpg",
    "auto-pappus": "clayworth/named_graphs/pappus-F018A.jpg",
    "auto-desargues": "clayworth/named_graphs/desargues-F020A.jpg",
    "auto-nauru": "clayworth/named_graphs/nauru-F024A.jpg",
    "auto-f26a": "clayworth/named_graphs/f26a-F026A.jpg",
    "auto-coxeter": "clayworth/named_graphs/coxeter-F028A.jpg",
    "auto-tuttecoxeter": "clayworth/named_graphs/tuttecoxeter-F030A.jpg",
    "auto-psl2_5_c3": "clayworth/named_graphs/psl2_5_c3.jpg",
    "auto-cr19": "clayworth/named_graphs/cr19.jpg",
    "auto-cr31": "clayworth/named_graphs/cr31.jpg",
    "auto-pg23": "clayworth/named_graphs/pg23.jpg",
    "auto-petersen_voltage": "clayworth/named_graphs/petersen_voltage.jpg",
    "auto-heawood_voltage": "clayworth/named_graphs/heawood_voltage.jpg",
    "auto-mk_voltage": "clayworth/named_graphs/mk_voltage.jpg",
    "auto-cube": "lean/named_graphs/cube-F008A.jpg",
    "auto-petersen": "lean/named_graphs/petersen-F010A.jpg",
    "auto-heawood": "lean/named_graphs/heawood-F014A.jpg",
    "auto-moebiuskantor": "lean/named_graphs/mobius-kantor-F016A.jpg",
    "auto-dodecahedron": "lean/named_graphs/dodecahedron-F020B.jpg",
    "auto-psl2_7_c3": "clayworth/named_graphs/klein-F056A.jpg",
    "auto-psl2_13": "lean/named_graphs/zhou3-F182A.jpg",
    "auto-biggssmith": "clayworth/named_graphs/biggssmith-F102A.jpg",
    "auto-gh22": "lean/named_graphs/langer-63v.jpg",
    "auto-tuttecage12": "lean/named_graphs/tutte-12cage.jpg",
}

PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{xcolor}
\newcommand{\C}{\mathbb{C}}
\newcommand{\Q}{\mathbb{Q}}
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
\texttt{docs/RESEARCH\_PROGRAM.md}. \emph{Hecke} paragraphs summarize the
Hecke algebra $\mathcal H(\mathrm{Aut}\,\Gamma,\mathrm{Stab}(v))$ of each
pair --- spherical function tables and explicit irreducible matrices ---
computed in the \texttt{hecke} repo and written up in
\texttt{hecke-irreps.pdf} (with \texttt{why13.pdf} for Meinhold-72).
\end{minipage}
\vfill
{\small Generated by \texttt{src/make\_catalogue.py} --- regenerate after adding drawings.\par}
\end{titlepage}

\tableofcontents
\newpage
"""

# Final summary page: every non-trivial (higher-dimensional) irreducible
# representation in the catalogue, with its graph and defining relations.
SUMMARY_PAGE = r"""\addtocontents{toc}{\protect\vspace{0.5em}}
\section{Non-trivial representations --- summary}
Every higher-dimensional irreducible representation occurring in this
catalogue. All other pairs are Gelfand: their Hecke algebras are commutative
and every irrep is $1$-dimensional.

\bigskip
\renewcommand{\arraystretch}{1.75}
\begin{center}
\begin{tabular}{@{}lcl@{}}
\hline
\textbf{Graph (pair)} & \textbf{Block} & \textbf{Defining relations} \\
\hline
F26A \quad $(C_{13}{:}C_6)/C_3$ & $M_2\times2$ &
  $T^2=\tfrac{5\mp\sqrt{13}}2I$, $N_+N_-=3$; Galois over $\Q(\sqrt{13})$ \\
Dodecahedron \quad $\mathrm{PSL}(2,5)/C_3$ & $M_2$ &
  $A^2=I$, \ $C^2+C+3I=0$, \ $(\rho(T_6)-I)^2=4I$ \\
\hline
Klein F056A \quad $\mathrm{PSL}(2,7)/C_3$ & $M_2$ &
  $D^2+D+2I=0$, \ $S^2=2I$, \ $SDS^{-1}=\bar D$ \\
 & $M_2$ & $B^2=2I$, \ $R^3=I$ \\
 & $M_3$ & integral, charpoly $x^3+x-6$; \ $0,\tfrac{-1\pm\sqrt{17}}2$ \\
\hline
CR(19) \quad $(C_{19}{:}C_6)/C_3$ & $M_2\times3$ &
  $T^2=N_aI$; cubic period field of $\Q(\zeta_{19})$ \\
CR(31) \quad $(C_{31}{:}C_6)/C_3$ & $M_2\times5$ &
  $T^2=N_aI$; quintic period field of $\Q(\zeta_{31})$ \\
\hline
Zhou-3 F182A \quad $\mathrm{PSL}(2,13)/S_3$ & $M_2\times3$ &
  Galois conjugate over $\Q(\zeta_7)^+$ \ ($c_k\mapsto c_{2k}$) \\
 & $M_3$ & Steinberg: $\rho(T_3)^3=-I$, \ charpoly $x^3+x-26$ \\
 & $M_3$ & $\rho(T_3)^3=5I$ \ (cubic analogue of $B^2=13I$) \\
\hline
Meinhold-72 \quad $G_2(2)/\mathrm{PSL}(3,2)$ & $M_2$ &
  $A^2=I$, \ $B^2=13I$, \ $AB+BA=-4I$ \\
\hline
Klein $\{7,3\}$ \quad $\mathrm{PGL}(2,7)/S_3$ & $M_2$ &
  integral: $x^2+x+4$, $x^2+x-4$ \ ($\tfrac{-1\pm\sqrt{17}}2$) \\
 & $M_2$ & $R^3=I$, \ $T^2=8I$ \\
\hline
\end{tabular}
\end{center}

\bigskip
The reflection generator $\rho(T_2)=\bigl(\begin{smallmatrix}1&0\\0&-1
\end{smallmatrix}\bigr)$ appears in every $2$-dimensional block: it is the
mirror involution of the pair. Explicit matrices are on the individual pages.
\newpage
"""

PAGE = r"""\section{%(title)s}
\begin{center}
\includegraphics[width=0.92\textwidth,height=0.66\textheight,keepaspectratio]{%(image)s}
\end{center}
\vfill
%(description)s

\medskip
{\small\color{gray}\texttt{%(source_tex)s}}
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


# Drawings deliberately left out of the catalogue (Robin, 2026-06-10).
# zhou-z2-quotient re-included 2026-06-12 after the cover-colour redraw;
# auto-meinhold72 superseded by meinhold72-rings (C7 ring layout).
EXCLUDED = {"zhou-quotient-spectral", "auto-meinhold72"}


def main():
    on_disk = {os.path.splitext(os.path.basename(p))[0]
               for p in glob.glob(os.path.join(INDIV, "*.pdf"))}
    on_disk -= EXCLUDED
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
            coloured = f"individual-coloured/{stem}.jpg"
            if stem in IMAGE_OVERRIDES:
                image = IMAGE_OVERRIDES[stem]
                source = IMAGE_OVERRIDES[stem]
                print(f"  using override {image}")
            elif os.path.exists(os.path.join(ROOT, coloured)):
                image = coloured
                source = f"individual-coloured/{stem}.pdf"
                print(f"  using coloured {image}")
            else:
                print(f"  rendering {stem}.jpg")
                pdf_to_jpeg(stem)
                image = f"individual_graphs/{stem}.jpg"
                source = f"individual_graphs/{stem}.pdf"
            if stem in HECKE:
                desc += ("\n\n\\medskip\n\\textbf{Hecke.} " + HECKE[stem])
            parts.append(PAGE % {
                "title": title,
                "image": image,
                "description": desc,
                "source_tex": tex_escape_stem(source),
            })
            n_pages += 1

    parts.append(SUMMARY_PAGE)
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
