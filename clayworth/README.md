# CSS Surface Codes in Lean 4

Formalisation of CSS quantum error-correcting codes from surface tilings,
following [Breuckmann & Terhal (arXiv:1506.04029)](https://arxiv.org/abs/1506.04029).

## Main results (sorry-free)

- **∂² = 0**: the chain complex condition for any `CellularSurface`
- **rank(d₁) = V − 1**: via kernel characterisation of connected graphs over F₂
- **rank(d₂) = F − 1**: via dual kernel characterisation with two-sides condition
- **k = 2g**: number of logical qubits equals the first Betti number

## Concrete instances

| Graph | Foster | Sabidussi | g | CSS code | Lean |
|-------|--------|-----------|---|----------|------|
| Heawood | F014A | Sab(G₄₂, C₃) | 1 | [[21, 2]] | `Heawood.lean` ✓ |
| Möbius-Kantor | F016A | Sab(GL(2,3), C₃) | 2 | [[24, 4]] | `Bolza.lean` ✓ |
| Klein | F056A | Sab(PSL(2,7), C₃) | 3 | [[84, 6]] | `KleinQuartic.lean` |
| Clayworth | — | Sab(G₂(2), C₃) | 505 | [[6048, 1010]] | `Clayworth.lean` |

✓ = all axioms verified by `decide` (zero sorry)

## Building

```bash
lake build
```

Requires Lean 4 and Mathlib (v4.29.1).

## References

- Breuckmann & Terhal, *Constructions and Noise Threshold of Hyperbolic Surface Codes*, [arXiv:1506.04029](https://arxiv.org/abs/1506.04029)
- Langer, *Symmetric Graphs and their Quotients*, [arXiv:1306.4798](https://arxiv.org/abs/1306.4798)
