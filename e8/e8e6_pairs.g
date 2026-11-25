# Follow-up: the (W(E8)-module dimension, multiplicity) pairs for H(W(E8),W(E6)).
Read("/tmp/e8_perms.g");;
G := Group(gens);;
H6 := Stabilizer(G, [r0, w0], OnTuples);;
pc := PermutationCharacter(G, H6);;
tbl := CharacterTable(G);;
irr := Irr(tbl);;
pairs := [];;
for i in [1..Length(irr)] do
  m := ScalarProduct(tbl, pc, irr[i]);
  if m > 0 then Add(pairs, [irr[i][1], m]); fi;
od;
Sort(pairs);
Print("(module dim, multiplicity = Hecke irrep dim):\n");
for p in pairs do Print("  ", p[1], "  x", p[2], "\n"); od;
# Wedderburn block counts
blocks := Collected(List(pairs, p -> p[2]));;
Print("Wedderburn: dim m appears #blocks:\n", blocks, "\n");
Print("dim check: ", Sum(pairs, p -> p[2]^2), " = 79\n");
