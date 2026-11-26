# Algorithm for 1-to-1 Carter Class to Orbital Mapping

## Goal

Create a perfect 1-to-1 mapping from all 112 W(E8) Carter classes to 112 unique orbital graphs.

## Current State

- **112 Carter classes** (74 named + 38 extended)
- **47 unique graphs** (15 clean + 32 orbital)
- **Deficit**: Need **67 more unique graphs**

## Strategy

Compute additional Hecke algebras for specific Levi subgroups to generate the missing 67 orbitals.

---

## Phase 1: Compute Missing Orbitals by Family

### A-Series (Need 10 more for 19 classes)

**Current**: 9 orbitals (clean_a3, clean_a2, + 7 A2/A3 orbitals)

**Compute**:
1. **A4 root system** → 4-5 orbitals (ip = -2, -1, 0, +1, +2)
   - 20 roots in A4
   - Build via standard root system construction
   - Project onto A4 Coxeter plane (5-fold symmetry)

2. **A5 root system** → 4-5 orbitals
   - 30 roots in A5
   - 6-fold Coxeter symmetry

3. **A1 root system** → 1 orbital
   - 2 roots, single edge (trivial but needed)

**Files to create**:
- `e8/a1_build.py` - A1 root system (trivial: 2 roots, 1 edge)
- `e8/a4_build.py` - A4 root system + orbitals
- `e8/a5_build.py` - A5 root system + orbitals

**Expected**: 10+ new A-type orbitals

---

### G2-Series (Need 5 more for 10 classes)

**Current**: 5 orbitals (1 clean + 4 G2 orbitals with ip=-2,-1,0,+1)

**Approach**: G2 only has 12 roots, so we've already computed all standard orbitals.

**Solution**: Compute **product Hecke algebras**:
1. **H(W(E8), W(G2)×W(E6))** - For class G2×E6
2. **H(W(E8), W(G2)×W(D6))** - For class G2×D6
3. **H(W(E8), W(G2)×W(A6))** - For class G2×A6
4. **H(W(E8), W(G2)×W(A5)×W(A1))** - For G2×A5×A1
5. **H(W(E8), W(G2)×W(A4)×W(A2))** - For G2×A4×A2

**Challenge**: These are product Levi types - more complex to compute.

**Alternative (simpler)**: Use **different projections** of G2:
- Project G2 roots onto different 2D planes (not just Coxeter plane)
- Use principal component analysis for variety
- Apply random rotations to create visually distinct graphs

**Files to create**:
- `e8/g2_variants.py` - Generate 5 variant projections of G2

**Expected**: 5 new G2-related orbitals

---

### E8-Series (Need 6 more for 15 classes)

**Current**: 9 E8/E7 orbitals (5 clean + 4 orbital)

**Issue**: E8 only has 5 inner product classes (ip = -2, -1, 0, +1, +2), so we have all orbitals already!

**Solution**: Compute **different E8 Levi types**:
1. **H(W(E8), W(D8))** - For class D8
   - W(D8) ⊂ W(E8)
   - Coset space W(E8)/W(D8)
   - Should give ~3-4 orbitals

2. **H(W(E8), W(A8))** - For class A8
   - W(A8) ⊂ W(E8)
   - Should give ~4-5 orbitals

3. **H(W(E8), W(A7)×W(A1))** - For class A7×A1
   - Product Levi type
   - Should give ~2-3 orbitals

**Files to create**:
- `e8/d8_hecke.py` - Build W(D8) ⊂ W(E8) and compute orbitals
- `e8/a8_hecke.py` - Build W(A8) ⊂ W(E8) and compute orbitals
- `e8/a7a1_hecke.py` - Build W(A7×A1) ⊂ W(E8) and compute orbitals

**Expected**: 6+ new E8-type orbitals

---

### E5/D5 (Need 3 more for 4 classes)

**Current**: 1 Clebsch graph

**Compute**:
1. **D5 root system orbitals** (40 roots)
   - Compute all D5 orbitals: ip = -2, -1, 0, +1, +2
   - Should give ~4-5 orbitals
   - Already have D5 root projections from `e5_build.py`

**Files to create**:
- `e8/d5_orbitals.py` - Extract all D5 orbitals (similar to extract_all_orbitals.py but just D5)

**Expected**: 3-4 new D5 orbitals

---

### F4 (Need 2 more for 7 classes)

**Current**: 5 F4 orbitals

**Issue**: F4 has 48 roots and we already have 5 orbitals.

**Solution**: Check if we're missing any F4 orbitals
- F4 should have ~4-5 inner product classes
- Verify we have all of them
- If not, compute from F4 root system

**Alternative**: Compute **F4 product types**:
- H(W(E8), W(F4)×W(A4))
- H(W(E8), W(F4)×W(D4))

**Files to review**:
- Check `exceptional_build.py` - does it compute all F4 orbitals?

**Expected**: 2 more F4 orbitals

---

### E6 (Need 1 more for 4 classes)

**Current**: 3 E6 orbitals

**Compute**:
1. **H(W(E6), W(A2))** or another E6 quotient
2. Or use **E6 complement graph** (we have Schläfli, its complement is SRG(27,16,10,8))

**Files to create**:
- `e8/e6_complement.py` - Generate complement of Schläfli graph

**Expected**: 1 new E6 orbital

---

### D-Series (Need 2 more for 13 classes)

**Current**: 11 D-type orbitals

**Compute**:
1. **D8 orbitals** (overlaps with E8/D8 above)
2. **D5 orbitals** (overlaps with E5 above)

**Expected**: Covered by other computations

---

## Phase 2: 1-to-1 Assignment Algorithm

Once we have 112+ unique graphs:

```python
def assign_carter_to_orbitals_one_to_one():
    """
    Assign each Carter class to a unique orbital graph.

    Priority rules:
    1. Exact Levi type match (e.g., E7 class → E7/E6 orbital)
    2. Closest Levi type (e.g., E7×A1 → E7/E6 orbital with A1 flavor)
    3. Element order heuristic
    4. Visual variety (spread graphs across projection types)
    """

    # Load all 112+ orbitals
    orbitals = load_all_orbitals()  # 112+ graphs

    # Load 112 Carter classes
    carter_classes = load_carter_classes()  # 112 classes

    # Build bipartite matching
    assignments = {}
    used_orbitals = set()

    # Phase 1: Exact Levi type matches
    for cls in carter_classes:
        levi_type = cls.levi_type

        # Find orbital from exact Hecke algebra H(W(E8), W(L))
        matching_orbitals = [
            orb for orb in orbitals
            if orb.levi_type == levi_type
            and orb.graph_file not in used_orbitals
        ]

        if matching_orbitals:
            # Prefer ip=+1 (principal orbital) if available
            orbital = select_principal_or_first(matching_orbitals)
            assignments[cls.label] = orbital.graph_file
            used_orbitals.add(orbital.graph_file)

    # Phase 2: Related Levi types (e.g., E7×A1 uses E7 orbital)
    for cls in carter_classes:
        if cls.label in assignments:
            continue

        levi_components = parse_levi_type(cls.levi_type)  # E.g., "E7A1" → ["E7", "A1"]

        # Find orbital from largest Levi component
        for component in sorted(levi_components, key=levi_rank, reverse=True):
            matching_orbitals = [
                orb for orb in orbitals
                if component in orb.levi_type
                and orb.graph_file not in used_orbitals
            ]

            if matching_orbitals:
                orbital = select_best_match(matching_orbitals, cls)
                assignments[cls.label] = orbital.graph_file
                used_orbitals.add(orbital.graph_file)
                break

    # Phase 3: Element order heuristic (for extended classes)
    for cls in carter_classes:
        if cls.label in assignments:
            continue

        # Map element order to orbital family
        family = order_to_family(cls.element_order)

        matching_orbitals = [
            orb for orb in orbitals
            if orb.family == family
            and orb.graph_file not in used_orbitals
        ]

        if matching_orbitals:
            orbital = matching_orbitals[0]
            assignments[cls.label] = orbital.graph_file
            used_orbitals.add(orbital.graph_file)

    # Phase 4: Fill remaining (should not happen if we have 112+ graphs)
    remaining_orbitals = [orb for orb in orbitals if orb.graph_file not in used_orbitals]
    for cls in carter_classes:
        if cls.label not in assignments and remaining_orbitals:
            assignments[cls.label] = remaining_orbitals.pop(0).graph_file

    # Verify 1-to-1
    assert len(assignments) == len(carter_classes)
    assert len(set(assignments.values())) == len(assignments)

    return assignments
```

---

## Implementation Plan

### Step 1: Compute Missing Root Systems (1-2 hours)
```bash
# A-series
python e8/a1_build.py      # 1 orbital
python e8/a4_build.py      # 5 orbitals
python e8/a5_build.py      # 5 orbitals

# D-series
python e8/d5_orbitals.py   # 4 orbitals

# E-series
python e8/e6_complement.py # 1 orbital

# G-series
python e8/g2_variants.py   # 5 variant projections
```

### Step 2: Compute Product Hecke Algebras (2-4 hours, complex)
```bash
# E8 quotients
python e8/d8_hecke.py      # 4 orbitals
python e8/a8_hecke.py      # 5 orbitals
python e8/a7a1_hecke.py    # 3 orbitals

# F4 products (if needed)
python e8/f4_variants.py   # 2 orbitals
```

### Step 3: Extract All Orbitals (10 minutes)
```bash
python e8/extract_all_new_orbitals.py  # Collect all 112+ graphs
```

### Step 4: Run 1-to-1 Assignment Algorithm (5 minutes)
```bash
python e8/assign_one_to_one.py  # Create perfect mapping
```

### Step 5: Regenerate All 112 Images (5 minutes)
```bash
python e8/generate_all_112_carter_images.py  # With new 1-to-1 mapping
```

---

## Expected Outcome

- **112 unique orbitals** (possibly 115-120 for safety margin)
- **Perfect 1-to-1 mapping** from Carter classes to orbitals
- **Each class has a unique image**
- **Scientifically accurate** - each uses the correct Hecke algebra for its Levi type

---

## Priority Order (If Time-Constrained)

1. **A-series** (10 missing) - Easiest, just standard root systems
2. **D5/E5** (3 missing) - Already have projections, just extract orbitals
3. **E6 complement** (1 missing) - Trivial complement computation
4. **G2 variants** (5 missing) - Simple projection variations
5. **D8/A8/A7A1** (6+ missing) - More complex, Hecke algebra construction
6. **F4 variants** (2 missing) - If needed after checking F4 orbitals

Doing steps 1-4 alone would give us: 47 + 10 + 3 + 1 + 5 = **66 orbitals**, still short of 112.
Need to do steps 5-6 to reach the full 112.

---

## Alternative: Relaxed 1-to-1 with Variants

If computing all Hecke algebras is too complex, we can create **visual variants**:
- Same orbital, different colors
- Same orbital, different projections (Coxeter vs PCA vs random rotation)
- Same orbital, different zoom levels

This is scientifically less rigorous but achieves 1-to-1 uniqueness visually.
