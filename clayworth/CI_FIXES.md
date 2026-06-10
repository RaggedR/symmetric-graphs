# Mathlib4 PR CI Fixes — 2026-05-21/22

## Summary

13 open PRs on leanprover-community/mathlib4. 4 had CI failures; all fixed
over multiple rounds. 12/13 now green, #39636 pending final CI run.

## Reviewer Comments

| PR | Reviewer | Comment |
|----|----------|---------|
| #39650 | **vihdzp** | "Surely if Mathlib is to have graphs for specific polytopes, we should first have those polytopes?" |
| #39649 | **wwylele** | Pointed out how to self-add the `LLM-generated` label |
| #39636 | **themathqueen** | Fixed the `depends on:` format in PR description; showed how to self-add labels |

## CI Failures and Fixes

### PR #39650 — dodecahedron, Petersen, 3-cube
**1 issue, 1 fix**

| # | Source | Issue | Fix |
|---|--------|-------|-----|
| 1 | `mk_all --check` | `Archive.lean` import out of alphabetical order (`NamedGraphs` before `MiuLanguage`) | Move import after `MiuLanguage` entries |

### PR #39653 — cellular surfaces, CSS quantum codes
**2 issues, 2 fixes (across 2 rounds)**

| # | Source | Issue | Fix |
|---|--------|-------|-----|
| 1 | `docBlame` linter | Missing docstring on `CellularSurface.nextPerm` | Add `/-- ... -/` doc comment |
| 2 | `mk_all --check` | `Mathlib.lean` import ordering (`CellularSurface` after `Colex`) | Move before `Colex` (alphabetical) |

### PR #39654 — CellularSurface instances (genus 1, 2, 3, 505)
**5 issues, 3 fix rounds**

| # | Source | Issue | Fix |
|---|--------|-------|-----|
| 1 | `docBlame` linter | Missing docstring on `nextPerm` (shared file with #39653) | Cherry-pick from #39653 |
| 2 | `linter.style.header` | Doc content inside copyright block, not in `/-! ... -/` module docstring (4 files) | Restructure: copyright block minimal, `/-! ... -/` after imports |
| 3 | `linter.style.header` | Double spaces from `  - ` markdown lists | Replace with `* ` lists |
| 4 | `linter.style.longLine` | Lines > 100 chars in data arrays (Klein: 9, Clayworth: 1802) | Wrap array and `++` concatenation lines |
| 5 | `linter.style.longFile` | ClayworthSurface.lean is 4048 lines (limit 1500) | `set_option linter.style.longFile 4200` |

Also: `Mathlib.lean` ordering (same as #39653, fixed separately on this branch).

### PR #39636 — ascending Pochhammer is of binomial type
**5 issues across 5 rounds — the `IsShiftEquivariant` saga**

| # | Source | Issue | Fix attempted | Result |
|---|--------|-------|---------------|--------|
| 1 | Build | `introN` failed: can't see through `IsShiftEquivariant` across modules | `unfold IsShiftEquivariant` | FAIL — `unfold` can't match applied form |
| 2 | Build | `unfold` failed on applied `IsShiftEquivariant backwardDiff` | `delta IsShiftEquivariant` | FAIL — `delta` doesn't recognize it |
| 3 | Build | `delta` failed: "did not delta reduce" | `simp only [IsShiftEquivariant]` | FAIL — "not an exposed body" |
| 4 | Build | `simp only` failed; term-mode `fun a =>` also fails | `@[expose] public section` in DeltaOperator.lean | Fixed the root cause |
| 5 | `simpNF` linter | `backwardDiff_X` and `backwardDiff_C` `@[simp]` now redundant (expose made `simp` too powerful) | Remove `@[simp]` attributes | Pending CI |

**Root cause**: Lean 4.30's `module` keyword makes definitions unexposed
by default. `@[expose] public section` is required for cross-module
definitional transparency.

## Other Changes

- Merged symmetric-graphs PRs #1, #2, #3 to main
- Created `lean/named_graphs/` with stable image URLs
- Added visualization links to 6 mathlib Archive files
- Renamed Zhou -> Zhou-3 everywhere, regenerated image
- Added Zhou-6 image and links
