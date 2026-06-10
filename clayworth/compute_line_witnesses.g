# Compute BFS witness words for G₂(2) line-action transitivity
# Same approach as G2Action.lean: BFS from vertex 0, record generator words

gen1_lines := PermList([5, 4, 6, 12, 10, 11, 18, 16, 17, 7, 8, 9, 2, 1, 3,
  19, 20, 21, 15, 13, 14, 31, 37, 52, 53, 55, 54, 58, 25, 59, 62, 24, 63, 56,
  34, 57, 60, 28, 61, 32, 29, 44, 46, 22, 40, 23, 50, 35, 48, 38, 41, 33, 36,
  49, 47, 27, 51, 26, 42, 39, 43, 30, 45]);;  # 1-indexed

gen2_lines := PermList([16, 26, 27, 37, 46, 47, 28, 61, 60, 39, 17, 30, 18, 15,
  9, 38, 41, 62, 29, 49, 54, 12, 31, 43, 58, 35, 6, 3, 2, 1, 45, 44, 24, 21,
  57, 63, 52, 22, 40, 10, 4, 36, 33, 42, 13, 53, 7, 56, 8, 51, 14, 11, 34, 5,
  32, 20, 19, 59, 25, 50, 48, 23, 55]);;  # 1-indexed

# Verify orders
Print("gen1 order: ", Order(gen1_lines), "\n");;
Print("gen2 order: ", Order(gen2_lines), "\n");;

# BFS from vertex 1 (0-indexed: vertex 0)
# Alphabet: 0=gen1, 1=gen2, 2=gen1^-1, 3=gen2^-1
gens := [gen1_lines, gen2_lines, gen1_lines^-1, gen2_lines^-1];;
genNames := ["0", "1", "2", "3"];;

visited := [1];;  # 1-indexed
words := rec();;
words.(1) := "[]";
queue := [1];;

while Length(queue) > 0 do
  v := queue[1];
  Remove(queue, 1);
  for gi in [1..4] do
    w := v^gens[gi];
    if not w in visited then
      Add(visited, w);
      Add(queue, w);
      # Word for w = word for v ++ [gi-1]
      if words.(v) = "[]" then
        words.(w) := Concatenation("[", genNames[gi], "]");
      else
        oldWord := words.(v);
        # Remove trailing ]
        inner := oldWord{[2..Length(oldWord)-1]};
        words.(w) := Concatenation("[", inner, ",", genNames[gi], "]");
      fi;
    fi;
  od;
od;

Print("Reached ", Length(visited), " vertices\n");;

# Output as Lean Array
Print("\nprivate def dualLangerWitnessWordData : Array (List (Fin 4)) := #[\n");;
for i in [1..63] do
  Print("  ", words.(i));
  if i < 63 then Print(","); fi;
  Print("\n");
od;
Print("]\n");;
