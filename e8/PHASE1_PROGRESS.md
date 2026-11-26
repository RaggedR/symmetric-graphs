# Phase 1 Progress: 50% to 1-to-1 Mapping

## ✅ Completed: A-Series Root Systems

Successfully built and extracted orbitals from:
- **A1**: 2 roots → 1 orbital
- **A4**: 20 roots → 4 orbitals  
- **A5**: 30 roots → 4 orbitals

**Total added**: 9 new orbitals

## 📊 Current Status

- **Before**: 47 orbitals (15 clean + 32 generated)
- **After**: 56 orbitals (15 clean + 41 generated)
- **Progress**: 50.0% (56/112 orbitals)
- **Deficit**: 56 more needed for perfect 1-to-1

## 📁 New Files Created

Orbitals:
- `orbital_a1_ipm2.pdf` (1 edge)
- `orbital_a4_ipm2.pdf` (10 edges)
- `orbital_a4_ipm1.pdf` (60 edges)
- `orbital_a4_ipp0.pdf` (60 edges)
- `orbital_a4_ipp1.pdf` (60 edges)
- `orbital_a5_ipm2.pdf` (15 edges)
- `orbital_a5_ipm1.pdf` (120 edges)
- `orbital_a5_ipp0.pdf` (180 edges)
- `orbital_a5_ipp1.pdf` (120 edges)

Code:
- `a_series_build.py` - Builder for A1, A4, A5 root systems

## 🎯 Next Steps (Phase 2)

To reach 70+ orbitals (~62% progress):

1. **D5 orbitals** (~30 min)
   - Extract from existing `/tmp/d5_*.npy` data
   - Expected: 4 new orbitals
   - Would bring us to: 60 orbitals

2. **E6 complement** (~15 min)
   - Compute complement of Schläfli graph
   - Expected: 1 new orbital
   - Would bring us to: 61 orbitals

3. **G2 variants** (~30 min)
   - Generate projection variants
   - Expected: 5 new orbitals
   - Would bring us to: 66 orbitals

**Phase 2 total**: ~1.25 hours → 66 orbitals (59% progress)

## 🚀 For Perfect 1-to-1 (Phase 3)

Would need an additional 46 orbitals via:
- E8 Hecke algebras (D8, A8, A7×A1)
- Product Levi types
- More complex computations

**Estimated time**: 3-4 additional hours

## 💡 Current Impact

With 56 orbitals for 112 classes:
- Average: Each graph used by ~2 classes
- Much better than before (45 graphs → ~2.5 classes/graph)
- Many classes now have unique graphs
- Remaining duplicates mostly in extended classes

---
Generated: 2026-11-26
Orbital count: 56/112 (50%)
Next milestone: 70 orbitals (62%)
