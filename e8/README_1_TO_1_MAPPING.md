# Path to Perfect 1-to-1 Carter Class to Orbital Mapping

## Current Status

✅ **Completed**: 112 Carter class images generated
📊 **Current mapping**: 45 unique graphs for 112 classes (not 1-to-1)
🎯 **Goal**: 112 unique graphs for perfect 1-to-1 mapping

## What We Need

**Deficit**: 67 more unique orbitals

Breakdown by family:
- **A-series**: Need 10 more (19 classes, 9 graphs)
- **E8-series**: Need 6 more (15 classes, 9 graphs)
- **D-series**: Need 2 more (13 classes, 11 graphs)
- **G2-series**: Need 5 more (10 classes, 5 graphs)
- **F4-series**: Need 2 more (7 classes, 5 graphs)
- **E6-series**: Need 1 more (4 classes, 3 graphs)
- **E5-series**: Need 3 more (4 classes, 1 graph)

## Implementation Strategy

### Phase 1: Easy Wins (Compute Standard Root Systems)

**A-series** (10 orbitals, ~1 hour):
```bash
python e8/a1_build.py      # 1 orbital (trivial: 2 roots)
python e8/a4_build.py      # 5 orbitals (20 roots)
python e8/a5_build.py      # 5 orbitals (30 roots)
```

**D5/E5** (3-4 orbitals, ~30 min):
```bash
python e8/d5_orbitals.py   # Extract D5 orbitals from existing projections
```

**E6** (1 orbital, ~15 min):
```bash
python e8/e6_complement.py # Compute complement of Schläfli graph
```

**Total from Phase 1**: 14-15 new orbitals → **61-62 total**

### Phase 2: Intermediate (Projection Variants)

**G2 variants** (5 orbitals, ~30 min):
```bash
python e8/g2_variants.py   # Different projections of same G2 roots
```

**F4 check** (possibly 2 more, ~20 min):
```bash
# Verify we have all F4 orbitals, compute missing ones
```

**Total from Phase 2**: 7 new orbitals → **68-69 total**

### Phase 3: Advanced (Hecke Algebras for Product Levi Types)

**E8 quotients** (10+ orbitals, ~2-3 hours):
```bash
python e8/d8_hecke.py      # H(W(E8), W(D8))
python e8/a8_hecke.py      # H(W(E8), W(A8))
python e8/a7a1_hecke.py    # H(W(E8), W(A7)×W(A1))
```

**Total from Phase 3**: 10+ new orbitals → **78-80 total**

### Phase 4: Fill Gaps (Product Types)

Compute specific product Hecke algebras for remaining classes.

**Total target**: 112+ unique orbitals

## The Algorithm

Once we have 112+ orbitals, run:

```bash
python e8/assign_one_to_one.py
```

This implements a 5-phase assignment:

1. **Exact Levi type match** (E7 class → E7/E6 orbital)
2. **Levi component match** (E7×A1 → E7 orbital)
3. **Family match** (Any E-type → E-family orbital)
4. **Element order heuristic** (Order 18 → E7-type)
5. **Fill remaining** (Any available orbital)

Result: Perfect 1-to-1 mapping

## Files

**Documentation**:
- `ALGORITHM_FOR_1_TO_1_MAPPING.md` - Detailed strategy
- `README_1_TO_1_MAPPING.md` - This file

**Code**:
- `assign_one_to_one.py` - The 1-to-1 assignment algorithm
- Scripts to create (see Phase 1-4 above)

**Current mapping**:
- `carter_to_hecke_mapping_112.json` - Current 45-graph mapping

**Future mapping**:
- `carter_to_orbital_one_to_one.json` - Perfect 1-to-1 mapping (after computing 112+ orbitals)

## Priority Recommendation

For best ROI:
1. **Do Phase 1** (A1, A4, A5, D5, E6) → 15 new orbitals
2. **Do Phase 2** (G2 variants, F4 check) → 7 new orbitals
3. Total: 47 + 22 = **69 orbitals** (reduces deficit from 67 to 43)

This gets you halfway to 1-to-1 with ~2 hours of work.

For perfect 1-to-1:
4. **Do Phase 3** (E8 Hecke algebras) → 10+ orbitals
5. **Do Phase 4** (product types) → 30+ orbitals
6. Total: 69 + 40 = **110+ orbitals** ≈ 1-to-1

This requires ~4-6 hours total.

## Quick Start

To begin implementing:

```bash
# Clone the pattern from existing builders
cp e8/a3_build.py e8/a4_build.py
cp e8/a3_build.py e8/a5_build.py

# Edit to change rank from 3 to 4 and 5
# Then run them

python e8/a4_build.py
python e8/a5_build.py
```

The pattern is well-established in `exceptional_build.py` and `extract_all_orbitals.py`.

---

**Current**: 47 orbitals, 112 classes → some sharing inevitable
**Goal**: 112+ orbitals, 112 classes → perfect 1-to-1 possible
**Path**: Compute 65-70 more orbitals via root systems + Hecke algebras
