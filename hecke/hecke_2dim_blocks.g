# hecke_2dim_blocks.g — When does H(PSL(2,q), C3) have a 2-dimensional block?
#
# A 2-dimensional Wedderburn block is an irreducible constituent chi of
# pi = Ind_{C3}^G 1 with multiplicity m_chi = 2.  By Frobenius reciprocity,
#
#     m_chi = < Ind_{C3} 1, chi > = (1/3)( chi(1) + chi(g) + chi(g^2) ),
#
# where g is an element of order 3.  So m_chi = dim of the C3-fixed subspace
# of chi.  This script prints, for each q, the location of the order-3 class
# (split torus if q=1 mod 3, non-split if q=2 mod 3, unipotent if q=3^k) and,
# for every irrep, [degree, chi(g), m_chi], flagging the 2-dimensional blocks.
#
# Result (see hecke-irreps.tex, "When does a 2-dimensional block exist?"):
# a 2-dim block occurs for an odd prime power q >= 5 iff q in {5, 7, 9}.
#
# Usage: gap -q hecke_2dim_blocks.g 2>&1 | tee hecke_2dim_blocks.out

analyse := function(q)
  local G, tbl, irr, ccl, c3, g2pos, i, chi, deg, val, m, rows, qmod3, loc, any;
  G := PSL(2,q);;
  tbl := CharacterTable(G);;
  irr := Irr(tbl);;
  ccl := ConjugacyClasses(tbl);;
  c3 := First([1..Length(ccl)], j -> Order(Representative(ccl[j]))=3);;
  if c3 = fail then Print("q=",q,": no order-3 element\n\n"); return; fi;
  # locate the class of g^2 (= g^{-1}); equals c3 iff the class is real
  g2pos := PowerMap(tbl, 2)[c3];;
  qmod3 := q mod 3;;
  if   qmod3 = 1 then loc := "split torus   (3 | q-1)";
  elif qmod3 = 2 then loc := "non-split torus (3 | q+1)";
  else                loc := "unipotent     (q = 3^k)"; fi;
  Print("================ q = ", q,
        "   (q mod 3 = ", qmod3, ",  g in ", loc, ") ================\n");
  Print("  deg   chi(g)   m_chi = dim fixed\n");
  rows := [];;
  for i in [1..Length(irr)] do
    chi := irr[i];;
    deg := chi[1];;
    val := chi[c3];;
    m   := (deg + chi[c3] + chi[g2pos]) / 3;;   # (chi(1)+chi(g)+chi(g^2))/3
    Add(rows, [deg, val, m]);
  od;
  Sort(rows, function(a,b) return a[1] < b[1]; end);;
  any := false;;
  for i in [1..Length(rows)] do
    Print("  ", String(rows[i][1],4), "  ", String(rows[i][2],6),
          "    ", rows[i][3]);
    if rows[i][3] = 2 and rows[i][1] > 1 then
      Print("   <== 2-DIM BLOCK"); any := true;
    fi;
    Print("\n");
  od;
  if any then Print("  ** q=",q," HAS a 2-dim block **\n");
         else Print("  (no 2-dim block)\n"); fi;
  Print("\n");
end;;

# finite check (q <= 21, i.e. min nontrivial degree (q-1)/2 <= 10) plus
# square-q corroboration q = 25, 49.
for q in [5, 7, 9, 11, 13, 17, 19, 23, 25, 29, 31, 49] do analyse(q); od;
Print("done\n");
