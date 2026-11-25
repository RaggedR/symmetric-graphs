# The Hecke algebra H(W(E8), W(E6)) for the skipped rung W(E8)/W(E6).
# W(E6) = stabiliser of the flag (r0, w0) = two roots at 60 degrees.
# The coset space is the 13440 ordered pairs of E8 roots at inner product +1.
#
# dim(Hecke) = rank of the permutation action = number of suborbits.
# commutative  <=>  the permutation character is multiplicity-free.
# irreps of the Hecke algebra are M_{m_i}(C) blocks; dims = the multiplicities m_i.

Read("/tmp/e8_perms.g");;
G := Group(gens);;
Print("|W(E8)| = ", Size(G), "  (expect 696729600)\n");

H7 := Stabilizer(G, r0);;
Print("|Stab(r0)| = ", Size(H7), "  (expect 2903040 = |W(E7)|)\n");
H6 := Stabilizer(G, [r0, w0], OnTuples);;
Print("|Stab(r0,w0)| = ", Size(H6), "  (expect 51840 = |W(E6)|)\n");
Print("index |G:H6| = ", Index(G, H6), "  (expect 13440)\n");

# coset action, rank and subdegrees
hom := FactorCosetAction(G, H6);;
P := Image(hom);;
Print("coset action degree = ", NrMovedPoints(P), "\n");
stab1 := Stabilizer(P, 1);;
suborbits := Orbits(stab1, [1..Index(G,H6)]);;
sub := SortedList(List(suborbits, Length));;
Print("RANK (= dim of Hecke algebra = #suborbits) = ", Length(suborbits), "\n");
Print("subdegrees (sorted) = ", sub, "\n");
# number of self-paired (symmetric) suborbits: suborbit O is self-paired iff
# g in O implies g^-1-partner in O; count via pairing on a base point
Print("sum of subdegrees = ", Sum(sub), " (should be 13440)\n");

# ---- multiplicities (commutativity + irrep dims) ----
Print("\n-- decomposing the permutation character --\n");
pc := PermutationCharacter(G, H6);;
Print("pc(1) = degree = ", pc[1], "\n");
tbl := CharacterTable(G);;
irr := Irr(tbl);;
mult := List(irr, chi -> ScalarProduct(tbl, pc, chi));;
nz := Filtered(mult, m -> m > 0);;
Print("multiplicities (nonzero) = ", SortedList(nz), "\n");
Print("number of distinct irreducible constituents = ", Length(nz), "\n");
Print("sum m_i          = ", Sum(nz), "\n");
Print("sum m_i^2 (=dim) = ", Sum(nz, m -> m^2), "\n");
Print("MULTIPLICITY-FREE (Hecke commutative)? ", ForAll(nz, m -> m = 1), "\n");
Print("irrep dimensions of the Hecke algebra (the m_i): ", SortedList(nz), "\n");
Print("dims of paired W(E8)-modules (chi(1) for constituents): ",
      SortedList(List(Filtered([1..Length(irr)], i -> mult[i] > 0), i -> irr[i][1])), "\n");
