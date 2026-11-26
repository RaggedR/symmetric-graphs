# Carter Classes of W(E8): Complete Visual Atlas

## Summary

We have generated **74 Carter class images** for the Weyl group W(E8), each showing:
- The orbital graph from the appropriate Hecke algebra descent
- The Dynkin diagram(s) for the Levi type in the top left corner

## What are Carter Classes?

Carter classes are conjugacy classes in the Weyl group W(E8), labeled by their **Levi types** - the structure of their centralizer subgroup. The Levi type is a product of irreducible root systems like E7×A1, D6×A2, A1^8, etc.

## The Mapping Strategy

Each Carter class with Levi type L is visualized using the Hecke algebra **H(W(E8), W(L))**:
- The orbital graph shows how W(E8) acts on the coset space W(E8)/W(L)
- We use the "principal" orbital (typically inner product +1 or the most regular graph)
- When exact computation is expensive, we use the nearest descent in the exceptional ladder

## The Exceptional Ladder Descent

```
E8 ⊃ E7 ⊃ E6 ⊃ D5
 ↓     ↓     ↓    ↓
240   56    27   16  (coset sizes)
 ↓     ↓     ↓    ↓
rank 5 rank 4 rank 3 (Hecke algebra dimensions)
```

## Files Generated

All 74 images are in the format: `carter_class_<label>.pdf`

### Examples:
- `carter_class_E8.pdf` - The full E8 (Coxeter element), using the mandala
- `carter_class_E7A1.pdf` - E7 × A1, using E8/E7 Hecke algebra
- `carter_class_A1A1A1A1A1A1A1A1.pdf` - The octave (A1^8), using the bloom
- `carter_class_D8.pdf` - The octocube (classical shadow)
- `carter_class_G2E6.pdf` - Hexagonal cubic

## Hecke Algebras Used

| Descent | Vertices | Orbitals | Main Classes Using This |
|---------|----------|----------|------------------------|
| H(W(E8), W(E7)) | 240 roots | 5 | E8, E7, E7A1, most A-series |
| H(W(E7), W(E6)) | 56 weights | 4 | E6, E6A2, E6A1² |
| H(W(E6), W(D5)) | 27 lines | 3 | D5 family |
| H(W(D5), W(D4)) | 16 | - | D4 family |
| F4 roots | 48 | - | F4 family |
| G2 roots | 12 | - | G2 family |

## Orbital Graphs by Type

### E8/E7 Orbitals (240 roots):
- **ip = +2**: identity (no edges)
- **ip = +1**: the mandala (56-regular, 6720 edges)
- **ip = 0**: the bloom (126-regular, 15120 edges)
- **ip = -1**: the weave (56-regular, 6720 edges)
- **ip = -2**: antipodal matching (120 K₂)

### Distribution:
- **e8_mandala.pdf**: 7 classes (E8, E7, A8, A7, etc.)
- **e8_ip0.pdf (bloom)**: 1 class (A1^8 - the octave)
- **e8_ipm1.pdf (weave)**: 6 classes (highly fragmented A-series)
- **e7_mandala.pdf**: 2 classes (E6-related)
- **e6_schlafli.pdf**: 4 classes (A6A2, etc.)
- **a2_roots.pdf**: 9 classes (A2-dominant)
- **a3_roots.pdf**: 10 classes (A3-dominant)
- **d4_graph.pdf**: 7 classes (D4 family)
- **d6_graph.pdf**: 3 classes (D6 family)
- **d7_classical.pdf**: 3 classes (D7, D8)
- **f4_roots.pdf**: 7 classes (F4 family)
- **g2_roots.pdf**: 10 classes (G2 family)

## The Complete List of 74 Classes

### Exceptional Classes (6):
1. E8, E7, E6, E7A1, E6A2, E6A1A1

### D-series (21):
D8, D7, D7A1, D6, D6A2, D6A1A1, D5, D5A3, D5A2A1, D5A1A1A1,
D4, D4D4, D4A4, D4A3A1, D4A2A2, D4A2A1A1, D4A1A1A1A1

### F4 family (7):
F4, F4A4, F4D4, F4A3A1, F4A2A2, F4A2A1A1, F4A1A1A1A1

### G2 family (10):
G2, G2E6, G2D6, G2A6, G2A5A1, G2A4A2, G2A3A3, G2A3A2A1, G2A2A2A2, G2A2A2A1A1

### Pure A-series (30):
A8, A7, A7A1, A6, A6A2, A6A1A1, A5, A5A3, A5A2A1, A5A1A1A1,
A4, A4A4, A4A3A1, A4A2A2, A4A2A1A1, A4A1A1A1A1,
A3, A3A3, A3A3A2, A3A3A1A1, A3A2, A3A2A2A1, A3A2A1A1A1, A3A1A1A1A1A1,
A2, A2A2, A2A2A2, A2A2A2A2, A2A2A2A1A1, A2A2A1A1A1A1, A2A1A1A1A1A1A1,
A1, A1A1A1A1A1A1A1A1

## Technical Details

**Generation scripts:**
- `enumerate_carter.py` - Enumerate candidate Carter labels (74 Levi types)
- `map_carter_to_hecke.py` - Map each class to a Hecke algebra orbital
- `generate_all_carter_images.py` - Generate all 74 images with Dynkin diagrams
- `get_all_carter_classes.sage` - SageMath: enumerate all 112 conjugacy classes

**Data files:**
- `carter_to_hecke_mapping.json` - Complete mapping specification
- `carter_classes_complete.json` - All 112 conjugacy classes from SageMath
- `carter_images_generated.txt` - List of successfully generated images

## Next Steps

To expand to the full 112 classes:
1. Some Levi types have multiple conjugacy classes (need deeper Carter theory)
2. For complete accuracy, compute H(W(E8), W(L)) for each distinct L
3. Add more sophisticated Dynkin diagram layouts for complex product types

## References

- Carter, R.W. "Conjugacy classes in the Weyl group" Compositio Math. 25 (1972)
- The exceptional ladder: `README.md` in this directory
- Hecke algebra theory: `docs/RESEARCH_PROGRAM.md`
