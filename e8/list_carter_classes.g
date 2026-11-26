# List all 112 conjugacy classes of W(E8) with their Carter labels
# (Levi types computed from centralizers)

# Load E8 permutation representation
Read("/tmp/e8_perms.g");;
G := Group(gens);;
Print("|W(E8)| = ", Size(G), " (expect 696729600)\n\n");

# Get all conjugacy classes
Print("Computing conjugacy classes of W(E8)...\n");
Print("(This may take a while)\n\n");
cc := ConjugacyClasses(G);;

Print("Found ", Length(cc), " conjugacy classes (expect 112)\n\n");

# For each class, print basic information
Print("Class | Size | Order | Carter Type (TODO)\n");
Print("------|------|-------|-------------------\n");

for i in [1..Length(cc)] do
    rep := Representative(cc[i]);;
    class_size := Size(cc[i]);;
    elt_order := Order(rep);;

    # TODO: Compute Carter label from centralizer structure
    # For now, just print the basic data

    Print(i, " | ", class_size, " | ", elt_order, " | ?\n");
od;

Print("\n");
Print("To get Carter labels, we need to:\n");
Print("1. Compute the centralizer C_G(g) for each class rep g\n");
Print("2. Identify C_G(g) as a parabolic subgroup W(L)\n");
Print("3. Extract the Levi type L\n");
Print("\nThis requires deeper reflection group machinery...\n");
