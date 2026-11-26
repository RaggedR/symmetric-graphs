# Dynkin Diagram Upgrade: Complete!

## Issues Fixed

### Issue 1: Text Labels Instead of Diagrams
**Before**: Named classes like A2, E8 showed text "A2", "E8"
**After**: Proper Dynkin diagrams rendered in TikZ

Examples:
- **A2**: Now shows o---o (2 nodes connected)
- **E8**: Now shows the full E8 Dynkin diagram (7 nodes + branch)
- **F4**: Now shows F4 diagram with proper structure
- **G2**: Now shows G2 diagram with triple bond

### Issue 2: Composite Types Showed Only One Diagram
**Before**: A2A2A2 showed only single A2 diagram
**After**: Shows 3 disjoint A2 diagrams stacked vertically

Examples:
- **A2A2A2**: 3 disjoint A2 diagrams
- **D4A3A1**: D4 (branched) + A3 (chain) + A1 (single node)
- **A2A2A1A1A1A1**: 2 × A2 + 4 × A1 (all disjoint)
- **E7A1**: E7 diagram + isolated A1
- **F4A2A2**: F4 + 2 disjoint A2

## Implementation

### Files Created

**`dynkin_diagrams.py`**
- `parse_levi_type()`: Parses composite types like "D4A3A1" → [("D", 4), ("A", 3), ("A", 1)]
- `single_dynkin_diagram()`: Generates TikZ for individual Dynkin diagrams
- `dynkin_diagram_tikz()`: Main function - handles composite types with vertical stacking

**`regenerate_improved_with_dynkin.py`**
- Regenerates all 60 Carter images with proper Dynkin diagrams
- Uses improved 1-to-1 mapping from assignment algorithm
- Integrates new `dynkin_diagrams.py` module

### Technical Details

**Diagram Types Supported:**
- **A_n**: Simple chain of n nodes
- **D_n**: Branching structure (fork at the end)
- **E_6, E_7, E_8**: Exceptional diagrams with branches
- **F_4**: Diagram with double bond
- **G_2**: Diagram with triple bond

**Composite Handling:**
- Components parsed via regex: `([A-Z])(\d+)`
- Stacked vertically with 2.0 unit spacing
- Each component rendered independently
- Example: "D4A3A1" → D4 at y=0, A3 at y=-2, A1 at y=-4

## Results

**Total Images**: 112 Carter class PDFs
**Upgraded**: 60 with perfect 1-to-1 + proper Dynkin diagrams
**Success Rate**: 100% (60/60 compiled successfully)

**Visual Quality:**
- All simple types show proper diagrams ✓
- All composite types show multiple disjoint diagrams ✓
- Diagrams positioned in top-left corner ✓
- White/gray color scheme matches graph background ✓

## Examples Verified

**Simple Types:**
```
A2.tex   → o---o
E8.tex   → o---o---o---o---o---o---o
                       |
                       o
```

**Composite Types:**
```
D4A3A1.tex →  D4: o---o
                       / \
                      o   o
              A3: o---o---o
              A1: o

A2A2A1A1A1A1.tex →  A2: o---o
                    A2: o---o
                    A1: o
                    A1: o
                    A1: o
                    A1: o
```

## Manifest

Updated manifest: `carter_improved_manifest_v2.txt`

Contains:
- Label → Graph file → PDF mapping
- 60 entries with perfect 1-to-1
- Composite diagram annotations

---

**Generated**: 2026-11-26
**Status**: Complete ✓
**Next**: Extend to remaining 52 classes when more orbitals available
