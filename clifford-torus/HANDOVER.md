# Clifford Torus — Handover Notes

## What we built

Three interactive Three.js visualizations showing the Clifford torus and how a flat square becomes one:

1. **clifford-torus.html** — Static view of the Clifford torus with two families of coordinate circles (u-lines in cyan/blue, v-lines in red/gold). Stereographic projection from S^3 to R^3. Animated green dot traces a dense irrational winding. Side-by-side with the flat square showing the same grid.

2. **homotopy.html** — Animated homotopy from flat square to Clifford torus via linear interpolation in R^4. Slider + play/pause. Pauses at t=0 (flat square), t~0.184 (cusp singularity, 5 seconds, orange label), and t=1 (Clifford torus, with faster spin). Bounces back and forth.

3. **natural-homotopy.html** — Holomorphic family of tori: (cos(alpha)*exp(2πiu), sin(alpha)*exp(2πiv)) on S^3, with alpha going from 0.18 to π/4. Always smooth, no cusp, but tau varies (shown on screen). Robin said this wasn't interesting — "we have no interest in a thin torus." Probably delete or don't link.

All three are live on GitHub Pages at `https://raggedr.github.io/symmetric-graphs/clifford-torus/`.

## The LinkedIn post (draft)

A post about why the Clifford torus formula is just "exponentiate each coordinate" — the exponential as the canonical covering map, the product structure, the category theory angle (natural transformation from additive to multiplicative group functor).

### Current state of the draft

The first part (exponential explanation) is solid and Robin is happy with it. The ending about the cusp singularity went through several rewrites. The cleanest version so far:

> Unfortunately, watching the wrapping happen is harder than describing the result. If you animate a flat square curling up into the Clifford torus by linearly interpolating in R^4, the surface stays flat throughout — the two coordinates live in orthogonal planes and don't interfere with each other.
>
> But at t = 1/(1 + pi*sqrt(2)), roughly 18.4% of the way through, the curve develops a cusp. The linear pull keeping it straight exactly balances the circular curl bending it closed. For one instant, at one point, the derivative vanishes. Then the curl wins, the curve closes into a circle, and the torus is born.

### Things Robin rejected

- "The square lives outside the category where the exponential is natural" — vague, doesn't mean anything clear
- "tau = i" without explanation — jargon for LinkedIn
- "The exponential can't produce a line segment" — backwards (the line segment is the input, the circle is the output)
- Any mention of the "thin torus" / natural homotopy — not interesting, distracts from the story
- LaTeX in conversation — doesn't render in the terminal

## Key mathematical facts

- **The Clifford torus map** is (u,v) -> (exp(2πiu), exp(2πiv))/√2 — a holomorphic embedding C/Z^2 -> C^2.

- **The linear homotopy** F(u,v,t) = (1-t)*(u-0.5, 0, v-0.5, 0) + t/√2*(cos2πu, sin2πu, cos2πv, sin2πv) decomposes as a product of two curves in orthogonal planes of R^4. This means K=0 and tau=i at every intermediate step (where the surface is non-singular).

- **The cusp** at t = 1/(1+π√2) ≈ 0.184 is where the derivative gamma'(1/4) = 0. At u=1/4, cos(2πu)=0 and sin(2πu)=1, giving (1-t) - 2πt/√2 = 0. The curve has a semicubical cusp (y ~ x^(2/3)) at that point. The cusp is unavoidable for any product-of-curves homotopy from an open line segment to a closed circle.

- **The holomorphic family** (a(t)*exp(2πiu), b(t)*exp(2πiv)) with a^2+b^2=1 is always smooth but has tau = i*b/a, which varies unless a=b (in which case nothing changes). You can't have holomorphic + tau preserved + shape changing — pick two.

- **3D vs 4D**: A torus in R^3 CAN have the right conformal structure (tau=i for the right R/r ratio), but it CANNOT be flat (K must be nonzero somewhere on a closed surface in R^3). The Clifford torus in R^4 is both conformally correct and flat.

## Files

```
clifford-torus/
├── clifford-torus.html    — static torus with grid (linked from LinkedIn)
├── homotopy.html           — square-to-torus animation (linked from LinkedIn)
├── natural-homotopy.html   — holomorphic family (probably don't link)
├── README.md               — math writeup
└── HANDOVER.md             — this file
```

## Repo

`~/git/graphviz/` → `git@github.com:RaggedR/symmetric-graphs.git`

GitHub Pages: `https://raggedr.github.io/symmetric-graphs/clifford-torus/`
