# Clifford Torus Visualizations

Interactive visualizations of the Clifford torus and its relationship to the flat square.

## Files

- **clifford-torus.html** — Static view of the Clifford torus (stereographic projection S^3 -> R^3) with two families of coordinate circles. [Live](https://raggedr.github.io/symmetric-graphs/clifford-torus/clifford-torus.html)

- **homotopy.html** — Animated homotopy from flat square to Clifford torus via linear interpolation in R^4. [Live](https://raggedr.github.io/symmetric-graphs/clifford-torus/homotopy.html)

## The homotopy

The map at time t is:

```
F(u,v,t) = (1-t) * (u-0.5, 0, v-0.5, 0)  +  t/sqrt(2) * (cos 2*pi*u, sin 2*pi*u, cos 2*pi*v, sin 2*pi*v)
```

At t=1 this is the holomorphic embedding (u,v) -> (exp(2*pi*i*u), exp(2*pi*i*v)) / sqrt(2) into C^2.

The u-part and v-part live in orthogonal planes of R^4 and decouple completely. The surface at every t is a product of two curves: gamma_t(u) in the (x1,x2)-plane and delta_t(v) in the (x3,x4)-plane. For any such product:

- The induced metric is ds^2 = |gamma'|^2 du^2 + |delta'|^2 dv^2 (no cross terms)
- E depends only on u, G depends only on v, so K = 0 everywhere
- Since gamma_t and delta_t are the same curve, the periods are equal, giving tau = i throughout

The intrinsic geometry is flat and the conformal structure (tau = i, the square lattice) is preserved at every intermediate step.

## The cusp singularity at t = 1/(1 + pi*sqrt(2)) ~ 0.184

At this value of t, the curve gamma_t develops a semicubical cusp at u = 1/4. The derivative gamma_t'(1/4) = 0: the linear pull (1-t) exactly cancels the circular curl (2*pi*t/sqrt(2)) at the steepest point of the cosine.

Near the cusp, the curve looks like y ~ x^(2/3). The metric degenerates at that single point. For all other values of t, the surface is a smooth immersion.

The animation pauses at this moment so you can see the shape of the surface when the singularity occurs. The cusp is an artefact of linear interpolation — a nonlinear homotopy (e.g. interpolating in polar coordinates) could avoid it entirely.

## 4D -> 3D projection

Both visualizations use stereographic projection from the north pole (0,0,0,1) of S^3:

```
(X, Y, Z) = (x1, x2, x3) / (1 - x4)
```

This is conformal (angle-preserving), so the orthogonality of the coordinate grid is visible in 3D. The torus appears asymmetric (fatter on one side) because stereographic projection distorts distances but preserves angles.
