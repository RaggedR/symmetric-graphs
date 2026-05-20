/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under MIT license.

# Binary Linear Codes over F₂

A binary linear code is the kernel of a parity check matrix over F₂ = ZMod 2.
This file defines the basic structure and connects code dimension to matrix rank
via the rank-nullity theorem.
-/

import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.LinearAlgebra.Dimension.RankNullity
import Mathlib.InformationTheory.Hamming

open Matrix

/-- A binary linear code of length `n` is specified by a parity check matrix
    `H : F₂^m → F₂^n`. The code is `C = ker(H) ⊆ F₂^n`. -/
structure BinaryLinearCode (n : ℕ) where
  /-- Number of parity checks (rows of H) -/
  m : ℕ
  /-- The parity check matrix H ∈ F₂^{m × n} -/
  parityCheck : Matrix (Fin m) (Fin n) (ZMod 2)

namespace BinaryLinearCode

variable {n : ℕ} (C : BinaryLinearCode n)

/-- The dimension of a binary linear code: `k = n - rank(H)`.
    By rank-nullity, this equals `dim(ker H)`. -/
noncomputable def dim : ℕ := n - C.parityCheck.rank

/-- The code as a submodule of F₂^n: the kernel of the parity check matrix. -/
noncomputable def codeSpace : Submodule (ZMod 2) (Fin n → ZMod 2) :=
  LinearMap.ker C.parityCheck.mulVecLin

/-- A codeword is an element of the code space. -/
def IsCodeword (v : Fin n → ZMod 2) : Prop :=
  C.parityCheck *ᵥ v = 0

theorem isCodeword_iff_mem_codeSpace (v : Fin n → ZMod 2) :
    C.IsCodeword v ↔ v ∈ C.codeSpace := by
  simp [IsCodeword, codeSpace, LinearMap.mem_ker]

end BinaryLinearCode
