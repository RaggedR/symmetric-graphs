# Handover: Completing the Carter Class Images (74 → 112)

**For:** Robin
**From:** Miriam + Claude
**Date:** 2026-11-26

## What We've Done

We've generated **74 Carter class images** for W(E8), each showing:
- A clean orbital graph from the appropriate Hecke algebra
- Dynkin diagram(s) for the Levi type in the top left corner
- No captions (clean, ready for use)

All files are in `e8/carter_class_*.pdf`.

## The Gap: 74 → 112 Classes

W(E8) has **112 conjugacy classes**, but we've only generated images for **74 Levi types**. The issue is that **some Levi types have multiple conjugacy classes**.

For example:
- The Levi type E7×A1 might correspond to 2-3 different conjugacy classes
- These are distinguished by how the central element embeds
- Carter's classification gives the full partition data

## What You Need to Do

### Step 1: Get the Complete Carter Labels

**Option A: Use SageMath (started but incomplete)**

We have `get_all_carter_classes.sage` which enumerates all 112 conjugacy classes:

```bash
cd e8
sage get_all_carter_classes.sage
```

This creates `carter_classes_complete.json` with all 112 classes. The issue: **we don't yet have Carter labels** (Levi types) for each. We only have:
- Element order
- Class size
- Reduced word

**What you need:** Extract the actual Carter labels from SageMath. Try:

```sage
W = WeylGroup("E8")
cc = W.conjugacy_classes_representatives()
for w in cc:
    # Get Carter label somehow...
    # Sage might have: w.carter_label() or similar
```

**Option B: Use GAP with CHEVIE package**

GAP with the CHEVIE package has complete Carter data for all exceptional groups:

```gap
# This is the ideal approach if you have CHEVIE installed
LoadPackage("chevie");
W := CoxeterGroup("E", 8);
cc := ConjugacyClasses(W);

for c in cc do
    carter := CarterLabel(c);  # or similar function
    Print(carter, "\n");
od;
```

**Option C: Use published tables**

Carter's original paper has the complete list:
- Carter, R.W. "Conjugacy classes in the Weyl group" Compositio Math. 25 (1972)
- Table for E8 lists all 112 classes with their labels

You could manually transcribe or find a digital version.

### Step 2: Map All 112 Classes to Hecke Algebras

Once you have the complete Carter labels, extend `carter_to_hecke_mapping.json`:

```python
# In map_carter_to_hecke.py
carter_to_hecke = {
    # ... existing 74 mappings ...

    # New classes (example - these are illustrative)
    "E7A1_v2": ("E8/E7", "clean_e8_mandala.pdf", "ip+1"),
    "D6A2_v2": ("E7/E6", "d6_graph.pdf", "roots"),
    # ... etc for all 112
}
```

**Key question:** For classes with the same Levi type, **which orbital graph to use?**

Possible strategies:
1. **Use the same graph** (simplest - just duplicate)
2. **Use different orbitals** from the same Hecke algebra (e.g., E7A1_v1 uses mandala, E7A1_v2 uses weave)
3. **Compute the specific Hecke algebra** for each conjugacy class representative

### Step 3: Compute Additional Hecke Algebras (Optional but Better)

For completeness, you might want to compute Hecke algebras for Levi types we're currently approximating:

**Currently missing:**
- H(W(E8), W(D8)) - we use D7 as proxy
- H(W(E8), W(E6)×W(A2)) - we use just E6
- Various product types

**How to compute a new Hecke algebra:**

```python
# Template: e8/e8e6_hecke.g
# 1. Build W(E8) permutation representation (we have this in /tmp/e8_perms.g)
# 2. Find the Levi subgroup W(L)
# 3. Compute coset space W(E8)/W(L)
# 4. Find orbital graphs by inner products
# 5. Project onto Coxeter plane and draw
```

Example for D8:
```gap
Read("/tmp/e8_perms.g");
G := Group(gens);  # W(E8)

# W(D8) is the stabilizer of... (need to identify D8 embedding)
# Then compute orbitals, etc.
```

### Step 4: Generate All 112 Images

Once you have the mapping:

```bash
cd e8
python3 map_carter_to_hecke.py          # Update mapping
python3 generate_all_carter_images.py    # Generate all 112 PDFs
```

## Files You'll Need

**Data:**
- `carter_classes_complete.json` - All 112 classes from SageMath (partial Carter labels)
- `carter_to_hecke_mapping.json` - Current mapping (74 classes)

**Scripts:**
- `get_all_carter_classes.sage` - SageMath: enumerate classes
- `map_carter_to_hecke.py` - Create the mapping
- `generate_all_carter_images.py` - Generate all images
- `create_clean_orbitals.py` - Make clean orbital graphs (no captions)

**Orbital graphs (clean, no captions):**
- `clean_e8_mandala.pdf`, `clean_e8_weave.pdf`, `clean_e8_bloom.pdf`, etc.
- Plus: `e7_mandala.pdf`, `e6_schlafli.pdf`, `d*_graph.pdf`, `f4_roots.pdf`, `g2_roots.pdf`, `a*_roots.pdf`

**Documentation:**
- `CARTER_CLASSES_README.md` - Overview of the 74 classes
- `README.md` - Hecke algebra theory
- `docs/SYMMETRIC_EMBEDDING.md` - Drawing algorithm

## Technical Deep Dive

### Carter Class Theory

A conjugacy class in W(E8) is determined by:
1. **Levi type**: The Weyl group W(L) of the centralizer
2. **Partition data**: How the eigenvalues split

Example: E7×A1 might have multiple classes because:
- The A1 factor can embed in different ways relative to E7
- This is captured by partition data beyond just the Levi type

### Our Mapping Strategy

```
Carter class with Levi L
    ↓
H(W(E8), W(L)) - Hecke algebra on coset space
    ↓
Select principal orbital (ip=+1 usually)
    ↓
Project onto Coxeter plane
    ↓
Draw with Dynkin diagram overlay
```

**Approximations we made:**
- Same Levi type → same orbital graph (loses distinction between multiple classes)
- Missing Levi subgroups → use nearest descent (D8 → D7, etc.)

### Where to Get Help

1. **SageMath documentation**: Search for "Weyl group conjugacy classes" and "Carter classes"
2. **GAP CHEVIE**: https://webusers.imj-prg.fr/~jean.michel/chevie/chevie.html
3. **Carter's paper**: The original classification (1972)
4. **Ask on MathOverflow**: "How to get Carter labels for all 112 W(E8) conjugacy classes in Sage/GAP?"

## Quick Win Strategy

If you need **all 112 images quickly** without full accuracy:

1. Take the 112 conjugacy classes from SageMath
2. For each, guess the Levi type from element order and class size
3. Map to the nearest of our 74 existing types
4. Generate 112 images (with duplicates)

This gives you the full set, even if some assignments are approximate. You can refine later.

## Questions?

The key bottleneck is **getting Carter labels for all 112 classes**. Everything else is mechanical - we have the pipeline working.

Good luck!

— Miriam & Claude
