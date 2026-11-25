"""
Create trading cards for Carter classes using Hecke algebra orbital graphs.

Each card shows:
- The class name and archetype
- A visualization (orbital graph or combination)
- Key properties (rank, order, Coxeter number)
- Poetic description
"""

import numpy as np
import subprocess

def create_trading_card_tex(
    class_label,
    archetype_name,
    quote,
    image_file,
    rank,
    properties,
    output_name
):
    """Generate a single trading card as a standalone LaTeX document"""

    tex = r"""
\documentclass[tikz,border=0pt]{standalone}
\usepackage{graphicx}
\usepackage{tikz}
\usetikzlibrary{shapes,shadows,calc}

\definecolor{cardback}{RGB}{10,10,25}
\definecolor{cardframe}{RGB}{150,100,200}
\definecolor{titletext}{RGB}{255,220,150}
\definecolor{bodytext}{RGB}{220,220,255}

\begin{document}
\begin{tikzpicture}

% Card background
\fill[cardback] (0,0) rectangle (8,11);

% Frame
\draw[cardframe, line width=2pt, rounded corners=5pt] (0.2,0.2) rectangle (7.8,10.8);
\draw[cardframe, line width=1pt, rounded corners=5pt] (0.4,0.4) rectangle (7.6,10.6);

% Title section
\fill[cardframe!30!cardback] (0.5,9.5) rectangle (7.5,10.5);
\node[titletext, font=\Large\bfseries] at (4,10) {CLASS_LABEL};

% Archetype name
\node[titletext, font=\large\itshape] at (4,9.2) {ARCHETYPE_NAME};

% Image section (centered square)
\node at (4,6) {\includegraphics[width=6cm]{IMAGE_FILE}};

% Properties box
\fill[cardframe!20!cardback] (0.6,3.8) rectangle (7.4,4.8);
\node[bodytext, font=\small, align=left, text width=6.5cm] at (4,4.3) {
\textbf{Rank:} RANK_VALUE \\
PROPERTIES_TEXT
};

% Quote section
\node[bodytext, font=\small\itshape, align=center, text width=6.5cm] at (4,2.5) {
``QUOTE_TEXT''
};

% Footer
\node[bodytext!50, font=\tiny] at (4,0.5) {Carter Class • W(E8) Conjugacy};

\end{tikzpicture}
\end{document}
"""

    # Replace placeholders
    tex = tex.replace('CLASS_LABEL', class_label)
    tex = tex.replace('ARCHETYPE_NAME', archetype_name)
    tex = tex.replace('IMAGE_FILE', image_file)
    tex = tex.replace('RANK_VALUE', str(rank))
    tex = tex.replace('PROPERTIES_TEXT', properties)
    tex = tex.replace('QUOTE_TEXT', quote)

    # Write tex file
    tex_path = f"/Users/miriam/symmetric-graphs/e8/card_{output_name}.tex"
    with open(tex_path, 'w') as f:
        f.write(tex)

    return tex_path


def compile_card(tex_path):
    """Compile a card to PDF"""
    result = subprocess.run(
        ['pdflatex', '-interaction=nonstopmode', tex_path],
        cwd='/Users/miriam/symmetric-graphs/e8',
        capture_output=True
    )
    return result.returncode == 0


# Create the main hybrid class cards

cards = [
    {
        'label': 'E8',
        'archetype': 'The Ouroboros',
        'quote': 'In my end is my beginning',
        'image': 'e8_mandala.pdf',
        'rank': 8,
        'properties': r'\textbf{Vertices:} 240 \\ \textbf{Symmetry:} 30-fold \\ \textbf{Type:} Exceptional Prime',
        'output': 'e8_ouroboros'
    },
    {
        'label': 'E7A1',
        'archetype': 'The Serpent \& the Spinor',
        'quote': 'The wyrm coiled around the axis mundi',
        'image': 'e7_mandala.pdf',
        'rank': 8,
        'properties': r'\textbf{Decomposition:} E7 $\times$ A1 \\ \textbf{Symmetry:} 18-fold (E7) \\ \textbf{Type:} Exceptional Hybrid',
        'output': 'e7a1_serpent'
    },
    {
        'label': 'E6A2',
        'archetype': 'The Cubic Triangle',
        'quote': 'Where form meets force',
        'image': 'e6_schlafli.pdf',
        'rank': 8,
        'properties': r'\textbf{Decomposition:} E6 $\times$ A2 \\ \textbf{Vertices:} 27 (E6) \\ \textbf{Type:} Exceptional + SU(3)',
        'output': 'e6a2_cubic'
    },
    {
        'label': 'D8',
        'archetype': 'The Octocube',
        'quote': 'The shadow cast by the snake',
        'image': 'd7_classical.pdf',  # Using D7 as proxy
        'rank': 8,
        'properties': r'\textbf{Type:} Classical (D-series) \\ \textbf{Symmetry:} 14-fold \\ \textbf{Vertices:} $2^7 = 128$',
        'output': 'd8_octocube'
    },
    {
        'label': 'A1$^8$',
        'archetype': 'The Octave',
        'quote': 'Eight bells ringing in different towers',
        'image': 'e8_ip0.pdf',  # The bloom - maximally dispersed
        'rank': 8,
        'properties': r'\textbf{Decomposition:} A1$^8$ \\ \textbf{Type:} Eight independent spins \\ \textbf{Character:} Most fragmentary',
        'output': 'a1_8_octave'
    },
    {
        'label': 'G2E6',
        'archetype': 'The Hexagonal Cubic',
        'quote': 'Six petals enfolding the trinity',
        'image': 'g2_roots.pdf',
        'rank': 8,
        'properties': r'\textbf{Decomposition:} G2 $\times$ E6 \\ \textbf{Type:} Hexagon + Cubic \\ \textbf{Rank:} 2 + 6',
        'output': 'g2e6_hexcubic'
    },
]

if __name__ == "__main__":
    print("Creating trading cards...\n")

    for card in cards:
        print(f"Creating card: {card['label']} - {card['archetype']}")
        tex_path = create_trading_card_tex(
            card['label'],
            card['archetype'],
            card['quote'],
            card['image'],
            card['rank'],
            card['properties'],
            card['output']
        )

        if compile_card(tex_path):
            print(f"  ✓ Compiled: card_{card['output']}.pdf")
        else:
            print(f"  ✗ Failed to compile")

    print("\n✓ All trading cards created!")
