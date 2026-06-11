# Clebsch Surface, E₆, and the Petersen-Schläfli Chain

The 27 lines on a cubic surface, the E₆ root system, and how the Petersen graph
falls out of the Schläfli graph via the Clebsch diagonal surface.

## Contents

| File | Description |
|------|-------------|
| [Clebsch.md](Clebsch.md) | Full mathematical writeup: blow-up picture, 27 lines, the Petersen → T(6) → Schläfli chain, E₆ connection, Hecke algebras of all three SRGs |
| `hecke_clebsch.g` | GAP script: computes the Schläfli graph as Sab(PSp(4,3), (C₂⁴):A₅, HaH), extracts T(6) and the Petersen graph, verifies SRG parameters and Hecke structure constants |
| `petersen-e6.html` | Interactive browser visualization of the Petersen graph inside the E₆ weight structure |

## The chain

```
Petersen             T(6)                Schläfli
SRG(10, 3, 0, 1)    SRG(15, 6, 1, 3)    SRG(27, 10, 1, 5)
Aut = S₅             Aut = S₆             Aut = W(E₆)
  |                    |                    |
W(A₄)               W(A₅)               W(E₆)
```

All three are Gelfand pairs. The Hecke algebras are commutative (≅ C³), and
the eigenvalue-1 eigenspaces nest: 5 ⊂ 9 ⊂ 20.

## Context

This directory was originally `hecke/` before the full Hecke algebra project
moved in. The content here is the geometric starting point — PSp(4,3) acting
on 27 cosets produces the Schläfli graph, which is where the E₆ connection
and the "why 13" question originate. See [`../hecke/`](../hecke/) for the
algebraic continuation.
