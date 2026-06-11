# hecke_catalog_full.g — COMPLETE data for hecke-catalog.tex:
#   for EVERY non-trivial irrep of every catalogued Hecke algebra,
#   (1) generators and complete relations (min polys + rewriting rules),
#   (2) explicit matrices for ALL orbital generators T_1..T_r,
#   (3) every rho(T_i) expanded as a word/linear combination in the
#       chosen block generators,
#   plus the table of all 1-dimensional characters.
# Emits LaTeX fragments to catalog-data/<key>.tex; console log = audit trail.
#
# Method: isotypic restriction (hecke_irreps3.g), exact cyclotomics,
# symbol dictionaries per case for LaTeX conversion.
#
# Usage: gap -q --quitonbreak hecke_catalog_full.g < /dev/null \
#            > hecke_catalog_full.out 2>&1

LoadPackage("grape");;

# ───────────────────────── shared machinery ─────────────────────────

suborbitdata := function(im, n)
    local st, orbs, mats, prs, o, A, L, p;
    st := Stabilizer(im, 1);;
    orbs := ShallowCopy(Orbits(st, [1..n]));;
    SortBy(orbs, Length);;
    mats := [];;  prs := [];;
    for o in orbs do
        A := NullMat(n, n);;
        L := Orbit(im, [1, o[1]], OnPairs);;
        for p in L do A[p[1]][p[2]] := 1; od;
        Add(mats, A);;  Add(prs, L);;
    od;
    return [mats, prs];
end;;

spin := function(v, mats)
    local V, i, M, w, changed;
    V := [v];;
    changed := true;;
    while changed do
        changed := false;;
        for i in [1..Length(V)] do
            for M in mats do
                w := V[i] * M;;
                if SolutionMat(V, w) = fail then
                    Add(V, w);; changed := true;;
                fi;
            od;
        od;
    od;
    return V;
end;;

rootsfan := function(f)
    local c, b0, b1, d;
    c := CoefficientsOfUnivariatePolynomial(f);;
    if Length(c) = 2 then return [ -c[1]/c[2] ]; fi;
    if Length(c) = 3 then
        b1 := c[2]/c[3];;  b0 := c[1]/c[3];;
        d := b1^2 - 4*b0;;
        return [ (-b1 + Sqrt(d))/2, (-b1 - Sqrt(d))/2 ];
    fi;
    return [];
end;;

# ───────────────────────── LaTeX machinery ─────────────────────────

# positive rational -> tex
qtex := function(q)
    if IsInt(q) then return String(q); fi;
    return Concatenation("\\tfrac{", String(NumeratorRat(q)), "}{",
                         String(DenominatorRat(q)), "}");
end;;

# is the string a single term (no internal +/-)?
issingle := function(s)
    local i;
    if Length(s) = 0 then return false; fi;
    if Position(s, '+') <> fail then return false; fi;
    for i in [2..Length(s)] do
        if s[i] = '-' then return false; fi;
    od;
    return true;
end;;

# rational linear combination of symbol tex strings ("" = 1)
combtex := function(coeffs, texs)
    local s, j, q, aq, sgn, piece;
    s := "";;
    for j in [1..Length(coeffs)] do
        q := coeffs[j];;
        if q = 0 then continue; fi;
        if q < 0 then sgn := "-";; aq := -q;; else sgn := "+";; aq := q;; fi;
        if texs[j] = "" then
            piece := qtex(aq);;
        elif aq = 1 then
            piece := texs[j];;
        else
            piece := Concatenation(qtex(aq), texs[j]);;
        fi;
        if s = "" then
            if sgn = "-" then s := Concatenation("-", piece);;
            else s := piece;; fi;
        else
            s := Concatenation(s, sgn, piece);;
        fi;
    od;
    if s = "" then return "0"; fi;
    return s;
end;;

# global warning counter for fallbacks
FALLBACKS := [];;

# cyclotomic -> tex via the case's symbol dictionary
enttex := function(x, case)
    local N, best, sub, rows, sol, nz, full;
    if x = 0 then return "0"; fi;
    N := Conductor(Concatenation(case.symvals, [x]));;
    best := fail;;
    for sub in case.subsets do
        rows := List(case.symvals{sub}, s -> CoeffsCyc(s, N));;
        sol := SolutionMat(rows, CoeffsCyc(x, N));;
        if sol <> fail and ForAll(sol, IsRat) then
            nz := Number(sol, c -> c <> 0);;
            if best = fail or nz < best.nz then
                best := rec(nz := nz, sol := sol, sub := sub);;
            fi;
        fi;
    od;
    if best = fail then
        Add(FALLBACKS, [case.key, x]);;
        return Concatenation("\\mathtt{", String(x), "}");
    fi;
    return combtex(best.sol, case.symtexs{best.sub});
end;;

# matrix -> pmatrix
mattex := function(M, case)
    local rows;
    rows := List(M, row ->
        JoinStringsWithSeparator(List(row, x -> enttex(x, case)), " & "));;
    return Concatenation("\\begin{pmatrix} ",
        JoinStringsWithSeparator(rows, " \\\\ "), " \\end{pmatrix}");
end;;

# cyclotomic-coefficient combination of word strings ("I" = identity word)
gencomb := function(coeffs, wtexs, case)
    local s, j, c, cs, sgn, piece;
    s := "";;
    for j in [1..Length(coeffs)] do
        c := coeffs[j];;
        if c = 0 then continue; fi;
        cs := enttex(c, case);;
        sgn := "+";;
        if cs = "1" then
            piece := wtexs[j];;
        elif cs = "-1" then
            sgn := "-";;  piece := wtexs[j];;
        elif issingle(cs) then
            if cs[1] = '-' then
                sgn := "-";;  cs := cs{[2..Length(cs)]};;
            fi;
            if wtexs[j] = "I" then
                piece := Concatenation(cs, "I");;
            else
                piece := Concatenation(cs, "\\,", wtexs[j]);;
            fi;
        else
            piece := Concatenation("\\bigl(", cs, "\\bigr)", wtexs[j]);;
        fi;
        if s = "" then
            if sgn = "-" then s := Concatenation("-", piece);;
            else s := piece;; fi;
        else
            s := Concatenation(s, sgn, piece);;
        fi;
    od;
    if s = "" then return "0"; fi;
    return s;
end;;

# minimal polynomial coefficients c_0..c_k (monic) of a matrix
minpolycoeffs := function(M)
    local pows, k, sol, cs;
    pows := [ Flat(M^0) ];;
    k := 1;;
    while true do
        sol := SolutionMat(pows, Flat(M^k));;
        if sol <> fail then
            cs := List([1..k], i -> -sol[i]);;
            Add(cs, 1);;
            return cs;
        fi;
        Add(pows, Flat(M^k));;
        k := k + 1;;
    od;
end;;

# minimal polynomial relation "X^k + ... = 0" in tex (delegates to gencomb,
# powers of the variable as "words", leading term first)
polytex := function(cs, var, case)
    local texs, j;
    texs := ["I"];;
    for j in [1..Length(cs) - 1] do
        if j = 1 then Add(texs, var);;
        else Add(texs, Concatenation(var, "^", String(j)));; fi;
    od;
    return Concatenation(gencomb(Reversed(cs), Reversed(texs), case),
                         " = 0");
end;;

# word (list of generator indices) -> tex via letters
wordtex := function(w, letters)
    local s, g;
    if Length(w) = 0 then return "I"; fi;
    s := "";;
    for g in w do s := Concatenation(s, letters[g]);; od;
    return s;
end;;

# --- rendered-width estimation (to break overlong align rows) ---

# read a {...} brace group starting at s[i] = '{'; return [contents, next i]
bracegroup := function(s, i)
    local d, j;
    d := 0;;  j := i;;
    repeat
        if s[j] = '{' then d := d + 1;
        elif s[j] = '}' then d := d - 1; fi;
        j := j + 1;;
    until d = 0;
    return [s{[i+1..j-2]}, j];
end;;

# approximate rendered width (in math-char units) of a simple math string;
# \tfrac stacks (max of num/den), \sqrt adds 1, symbol commands count 1
effw := function(s)
    local w, i, n, j, cmd, g1, g2;
    w := 0;;  i := 1;;  n := Length(s);;
    while i <= n do
        if s[i] = '\\' then
            j := i + 1;;
            while j <= n and (s[j] in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") do
                j := j + 1;;
            od;
            cmd := s{[i+1..j-1]};;
            i := j;;
            if cmd = "" then
                i := i + 1;;          # \, \; etc: zero-width, skip symbol
            elif cmd = "tfrac" then
                g1 := bracegroup(s, i);;
                g2 := bracegroup(s, g1[2]);;
                w := w + Maximum(effw(g1[1]), effw(g2[1])) + 1;;
                i := g2[2];;
            elif cmd = "sqrt" then
                g1 := bracegroup(s, i);;
                w := w + effw(g1[1]) + 1;;
                i := g1[2];;
            elif cmd = "bigl" or cmd = "bigr" then
                w := w + 0;;          # delimiter counted when ( ) reached
            else
                w := w + 1;;          # \eta, \alpha, \omega, ...: one symbol
            fi;
        elif s[i] in "{}_^ " then
            i := i + 1;;
        else
            w := w + 1;;  i := i + 1;;
        fi;
    od;
    return w;
end;;

# pmatrix width: sum of column widths (max entry per column) + spacing/parens
matw := function(M, case)
    local ncols, j, wsum;
    ncols := Length(M[1]);;
    wsum := 2 + 2 * (ncols - 1);;
    for j in [1..ncols] do
        wsum := wsum + Maximum(List([1..Length(M)],
            i -> effw(enttex(M[i][j], case))));;
    od;
    return wsum;
end;;

# total row budget in math-char units before we split matrix and word
ROWCAP := 58;;

# BFS monomial closure of a generator set, capped at msq dims
closure := function(gens, msq)
    local words, mats, flats, idx, g, M2;
    if Length(gens) = 0 then
        return rec(words := [[]], mats := [], flats := []);
    fi;
    words := [ [] ];;
    mats := [ gens[1]^0 ];;
    flats := [ Flat(gens[1]^0) ];;
    idx := 1;;
    while idx <= Length(mats) and Length(mats) < msq do
        for g in [1..Length(gens)] do
            M2 := mats[idx] * gens[g];;
            if SolutionMat(flats, Flat(M2)) = fail then
                Add(words, Concatenation(words[idx], [g]));;
                Add(mats, M2);;
                Add(flats, Flat(M2));;
                if Length(mats) = msq then break; fi;
            fi;
        od;
        idx := idx + 1;;
    od;
    return rec(words := words, mats := mats, flats := flats);
end;;

GENLETTERS := ["A", "B", "C", "D", "E"];;

# ───────────────────────── main per-case driver ─────────────────────────

docase := function(case)
    local im, n, tbl, irr, ccl, csums, C, M, g, x, sd, mats, prs, r, lens,
          pivals, pi, onedims, blocks, idx, chi, m, P, W, d, B, vals, a,
          cp, f, lam, E, v, V, rho, rhos, evecs, mu, i, j, c2, S, w,
          tex, gidx, gmats, cl, cdim, test, letters, sols, sol, wt, mt, out,
          mp, relstrs, nw, pos, prod, bi, fname, cols, onetex, k, ws;
    Print("############ ", case.name, " ############\n");;
    im := case.build();;
    n := case.n;;
    tbl := CharacterTable(im);;
    irr := Irr(tbl);;
    ccl := ConjugacyClasses(tbl);;
    csums := [];;
    for C in ccl do
        M := NullMat(n, n);;
        for g in Elements(C) do
            for x in [1..n] do M[x][x^g] := M[x][x^g] + 1; od;
        od;
        Add(csums, M);;
    od;
    sd := suborbitdata(im, n);;
    mats := sd[1];;  prs := sd[2];;
    r := Length(mats);;
    lens := List(mats, M -> Sum(M[1]));;
    Print("rank ", r, ",  suborbit lengths ", lens, "\n");;
    pivals := List(ccl, C -> Number([1..n], x -> x^Representative(C) = x));;
    pi := ClassFunction(tbl, pivals);;
    onedims := [];;
    blocks := [];;
    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m = 0 then continue; fi;
        P := (chi[1] / Size(im)) * Sum([1..Length(ccl)],
                 i -> ComplexConjugate(chi[i]) * csums[i]);;
        if m = 1 then
            vals := List(prs, L -> Sum(L, p -> P[p[2]][p[1]]) / chi[1]);;
            Add(onedims, rec(idx := idx, deg := chi[1], vals := vals));;
            Print("1-dim on chi_", idx, " (deg ", chi[1], "): ",
                  vals, "\n");;
            continue;;
        fi;
        if m > case.mmax then
            Print("chi_", idx, " mult ", m, " — SKIPPED (over mmax)\n");;
            continue;;
        fi;
        Print(">>> ", m, "-dim block on chi_", idx, " (deg ", chi[1],
              ")\n");;
        W := BaseMat(P);;
        d := Length(W);;
        B := List(mats, M -> List(W, w -> SolutionMat(W, w * M)));;
        V := fail;;
        for a in [2..r] do
            cp := CharacteristicPolynomial(B[a]);;
            for f in Set(Factors(cp)) do
                for lam in rootsfan(f) do
                    E := NullspaceMat(B[a] - lam * IdentityMat(d));;
                    if Length(E) <> chi[1] then continue; fi;
                    v := E[1];;
                    V := spin(v, B);;
                    if Length(V) = m then break; else V := fail; fi;
                od;
                if V <> fail then break; fi;
            od;
            if V <> fail then break; fi;
        od;
        if V = fail then
            Print("  !! NO MODULE FOUND — needs attention\n");;
            continue;;
        fi;
        rho := List(V, w -> SolutionMat(V, w * B[a]));;
        cp := CharacteristicPolynomial(rho);;
        evecs := [];;
        for f in Set(Factors(cp)) do
            for mu in rootsfan(f) do
                Append(evecs, NullspaceMat(rho - mu * IdentityMat(m)));;
            od;
        od;
        if Length(evecs) = m then
            V := List(evecs, e -> e * V);;
        fi;
        if m = 2 then
            for i in [2..r] do
                rho := List(V, w -> SolutionMat(V, w * B[i]));;
                if rho[1][2] <> 0 and rho[2][1] <> 0 then
                    c2 := rho[1][2] / rho[2][1];;
                    if IsRat(c2) and c2 > 0 then
                        V[2] := Sqrt(c2) * V[2];;
                    fi;
                    break;;
                fi;
            od;
        fi;
        rhos := List(B, M -> List(V, w -> SolutionMat(V, w * M)));;
        # sanity: algebra dim and rho(sum T_i) = 0
        S := [];;
        for i in [1..r] do
            w := Flat(rhos[i]);;
            if Length(S) = 0 or SolutionMat(S, w) = fail then Add(S, w); fi;
            for j in [1..r] do
                w := Flat(rhos[i] * rhos[j]);;
                if SolutionMat(S, w) = fail then Add(S, w); fi;
            od;
        od;
        if Length(S) <> m^2 then
            Print("  !! algebra dim ", Length(S), " <> ", m^2, "\n");;
        fi;
        if Sum(rhos) <> 0 * rhos[1] then
            Print("  !! rho(sum T_i) <> 0\n");;
        fi;
        # greedy minimal generating set
        gidx := [];;  cdim := 1;;
        for i in [2..r] do
            if cdim = m^2 then break; fi;
            test := closure(rhos{Concatenation(gidx, [i])}, m^2);;
            if Length(test.mats) > cdim then
                Add(gidx, i);;
                cdim := Length(test.mats);;
                cl := test;;
            fi;
        od;
        if cdim <> m^2 then
            Print("  !! generators do not reach M_", m, "\n");;
            continue;;
        fi;
        gmats := rhos{gidx};;
        letters := GENLETTERS{[1..Length(gidx)]};;
        Print("  generators: ", List([1..Length(gidx)],
              j -> Concatenation(letters[j], "=T_",
                                 String(gidx[j]))), "\n");;
        ws := List(cl.words, w -> wordtex(w, letters));;
        Print("  monomial basis: ", ws, "\n");;
        # expansions of all T_i
        sols := [];;
        for i in [1..r] do
            sol := SolutionMat(cl.flats, Flat(rhos[i]));;
            if sol = fail then Print("  !! T_", i, " not in span\n");; fi;
            Add(sols, sol);;
            # verify exactly
            if sol <> fail and
               Sum([1..m^2], k -> sol[k] * cl.mats[k]) <> rhos[i] then
                Print("  !! expansion of T_", i, " WRONG\n");;
            fi;
        od;
        # relations: min polys + rewriting rules
        relstrs := [];;
        for j in [1..Length(gidx)] do
            mp := minpolycoeffs(gmats[j]);;
            Add(relstrs, polytex(mp, letters[j], case));;
        od;
        for bi in [1..Length(cl.words)] do
            for j in [1..Length(gidx)] do
                nw := Concatenation(cl.words[bi], [j]);;
                if nw in cl.words then continue; fi;
                if Length(nw) = 2 and nw[1] = nw[2] then continue; fi;
                prod := cl.mats[bi] * gmats[j];;
                sol := SolutionMat(cl.flats, Flat(prod));;
                Add(relstrs, Concatenation(wordtex(nw, letters), " = ",
                    gencomb(sol, ws, case)));;
            od;
        od;
        Print("  relations: ", Length(relstrs), " emitted\n");;
        Add(blocks, rec(idx := idx, deg := chi[1], m := m,
            rhos := rhos, gidx := gidx, letters := letters,
            cl := cl, ws := ws, sols := sols, relstrs := relstrs));;
    od;
    # ───── emit LaTeX ─────
    tex := "% auto-generated by hecke_catalog_full.g — do not edit\n";;
    # 1-dim table
    if Length(onedims) > 0 then
        Append(tex, "\\paragraph{One-dimensional characters.}\n\\[\n");;
        Append(tex, "\\begin{array}{c|");;
        Append(tex, Concatenation(List(onedims, o -> "c")));;
        Append(tex, "}\n");;
        Append(tex, Concatenation(List(onedims, o ->
            Concatenation(" & \\chi_{", String(o.idx), "}\\,(",
                          String(o.deg), ")"))));;
        Append(tex, " \\\\\\hline\n");;
        for i in [1..r] do
            Append(tex, Concatenation("T_{", String(i), "}"));;
            for j in [1..Length(onedims)] do
                Append(tex, Concatenation(" & ",
                    enttex(onedims[j].vals[i], case)));;
            od;
            Append(tex, " \\\\\n");;
        od;
        Append(tex, "\\end{array}\n\\]\n\n");;
    fi;
    # blocks
    for bi in [1..Length(blocks)] do
        x := blocks[bi];;
        Append(tex, Concatenation(
            "\\subsubsection*{Block ", String(bi), ": $M_",
            String(x.m), "(\\C)$ on $\\chi_{", String(x.idx),
            "}$ (degree ", String(x.deg), ")}\n"));;
        # generators
        Append(tex, "\\noindent\\textbf{Generators.}\n\\[\n");;
        Append(tex, JoinStringsWithSeparator(
            List([1..Length(x.gidx)], j -> Concatenation(
                x.letters[j], " = \\rho(T_{", String(x.gidx[j]), "}) = ",
                mattex(x.rhos[x.gidx[j]], case))),
            ", \\qquad\n"));;
        Append(tex, " .\n\\]\n");;
        # relations
        Append(tex, Concatenation(
            "\\textbf{Complete relations} (monomial basis $",
            JoinStringsWithSeparator(x.ws, ",\\;"), "$ of $M_",
            String(x.m), "(\\C)$):\n\\begin{gather*}\n"));;
        Append(tex, JoinStringsWithSeparator(x.relstrs, " ,\\\\\n"));;
        Append(tex, " .\n\\end{gather*}\n");;
        # all orbital images
        Append(tex,
            "\\textbf{All orbital images} ($\\rho(T_i)$, valency $k_i$, ");;
        Append(tex, "word in the generators):\n\\begin{align*}\n");;
        for i in [1..r] do
            mt := mattex(x.rhos[i], case);;
            wt := gencomb(x.sols[i], x.ws, case);;
            # 20 units overhead: "\rho(T_{nn}) =", " = ", "(k_{nn} = v)"
            if matw(x.rhos[i], case) + effw(wt) + 20 > ROWCAP then
                # too wide: word on an aligned continuation row
                Append(tex, Concatenation(
                    "\\rho(T_{", String(i), "}) &= ", mt,
                    " & (k_{", String(i), "} &= ", String(lens[i]),
                    ") \\\\\n &= ", wt, " &&"));;
            else
                Append(tex, Concatenation(
                    "\\rho(T_{", String(i), "}) &= ", mt,
                    " = ", wt,
                    " & (k_{", String(i), "} &= ", String(lens[i]), ")"));;
            fi;
            if i < r then Append(tex, " \\\\\n"); fi;
        od;
        Append(tex, " .\n\\end{align*}\n\n");;
    od;
    fname := Concatenation("/Users/robin/git/hecke/catalog-data/",
                           case.key, ".tex");;
    out := OutputTextFile(fname, false);;
    SetPrintFormattingStatus(out, false);;   # no line-wrap backslashes
    PrintTo(out, tex);;
    CloseStream(out);;
    Print("  -> wrote ", fname, "\n\n");;
end;;

# ───────────────────────── case constructions ─────────────────────────

al  := (-1 + ER(-7)) / 2;;
alb := (-1 - ER(-7)) / 2;;
eta19 := a -> Sum([1,7,8,11,12,18], s -> E(19)^(a*s));;
eta31 := a -> Sum([1,5,6,25,26,30], s -> E(31)^(a*s));;
cc := k -> E(7)^k + E(7)^(-k);;

frobpair := function(p, mult)
    local trans, m, F, H;
    trans := PermList(List([0..p-1], x -> ((x+1) mod p) + 1));;
    m := PermList(List([0..p-1], x -> ((mult*x) mod p) + 1));;
    F := Group(trans, m);;
    H := Subgroup(F, [m^2]);;
    return Image(FactorCosetAction(F, H));
end;;

psl2c3 := function(q)
    local G, H;
    G := PSL(2, q);;
    H := Group(Representative(First(ConjugacyClasses(G),
            c -> Order(Representative(c)) = 3)));;
    return Image(FactorCosetAction(G, H));
end;;

Read("/Users/robin/git/hecke/klein_adj.g");;

kleinaut := function(adj)
    local n, gamma;
    n := Length(adj);;
    gamma := Graph(Group(()), [1..n], OnPoints,
                   function(x, y) return y in adj[x]; end, true);;
    return AutGroupGraph(gamma);
end;;

f26aaut := function()
    local code, full, m, edges, i, j, dir, gamma;
    code := [-7, 7];;
    full := Concatenation(List([1..13], r -> code));;
    m := Length(full);;
    edges := [];;
    for i in [1..m] do
        Add(edges, [i, (i mod m) + 1]);;
        j := ((i - 1 + full[i]) mod m) + 1;;
        if i < j then Add(edges, [i, j]); else Add(edges, [j, i]); fi;
    od;
    edges := Set(edges);;
    dir := Concatenation(edges, List(edges, e -> [e[2], e[1]]));;
    gamma := EdgeOrbitsGraph(Group(()), dir, 26);;
    return AutGroupGraph(gamma);
end;;

zhou3 := function()
    local G, emb;
    G := PSL(2, 13);;
    emb := First(IsomorphicSubgroups(G, SymmetricGroup(3)), e -> true);;
    return Image(FactorCosetAction(G, Image(emb)));
end;;

meinhold := function()
    local G63, K;
    G63 := PrimitiveGroup(63, 2);;
    K := Group(
      ( 1,33)( 2,61)( 3,57)( 4,27)( 5,47)( 6,36)( 7,60)( 8,32)( 9,59)(10,54)
        (12,37)(13,55)(14,29)(16,48)(17,20)(18,26)(19,31)(21,40)(23,56)(24,63)
        (25,39)(28,51)(30,41)(35,44)(38,58)(42,46)(43,62)(49,53),
      ( 1,12)( 2,56)( 3,55)( 4,51)( 5,59)( 6,61)( 7,39)( 8,34)( 9,57)(10,54)
        (11,22)(13,47)(14,38)(15,31)(16,36)(18,46)(19,52)(20,30)(23,48)(24,53)
        (25,27)(26,37)(28,60)(32,45)(33,42)(35,49)(43,63)(44,62),
      ( 2, 4)( 3,27)( 5,28)( 6, 9)( 7,36)( 8,42)(10,20)(11,14)(13,48)(15,44)
        (16,39)(17,41)(18,46)(19,24)(21,50)(22,38)(23,47)(25,55)(26,32)(30,54)
        (31,62)(33,34)(37,45)(43,63)(51,56)(52,53)(57,61)(59,60),
      ( 1,38)( 2,51)( 3,20)( 4,43)( 5,41)( 6,35)( 7,53)( 8,40)(10,31)(11,50)
        (12,37)(13,55)(14,26)(15,52)(16,25)(17,57)(18,29)(19,54)(21,32)(23,24)
        (27,62)(28,61)(30,47)(33,58)(36,44)(39,48)(49,60)(56,63) );;
    return Image(FactorCosetAction(G63, K));
end;;

cases := [
rec(key := "dodeca", name := "PSL(2,5)/C3 (n=20)", n := 20, mmax := 3,
    build := function() return psl2c3(5); end,
    symvals := [1, ER(5)], symtexs := ["", "\\sqrt5"],
    subsets := [[1,2]]),
rec(key := "f26a", name := "F26A (C13:C6)/C3 (n=26)", n := 26, mmax := 3,
    build := f26aaut,
    symvals := [1, ER(13)], symtexs := ["", "\\sqrt{13}"],
    subsets := [[1,2]]),
rec(key := "cr19", name := "CR(19) (C19:C6)/C3 (n=38)", n := 38, mmax := 3,
    build := function() return frobpair(19, 8); end,
    symvals := [1, eta19(1), eta19(2), eta19(4)],
    symtexs := ["", "\\eta_1", "\\eta_2", "\\eta_4"],
    subsets := [[1,2,3],[1,2,4],[1,3,4]]),
rec(key := "cr31", name := "CR(31) (C31:C6)/C3 (n=62)", n := 62, mmax := 3,
    build := function() return frobpair(31, 26); end,
    symvals := [1, eta31(1), eta31(2), eta31(3), eta31(4), eta31(8)],
    symtexs := ["", "\\eta_1", "\\eta_2", "\\eta_3", "\\eta_4", "\\eta_8"],
    subsets := [[1,2,3,4,5],[1,2,3,4,6],[1,2,3,5,6],
                [1,2,4,5,6],[1,3,4,5,6]]),
rec(key := "f056a", name := "F056A PSL(2,7)/C3 (n=56)", n := 56, mmax := 3,
    build := function() return psl2c3(7); end,
    symvals := [1, al, alb], symtexs := ["", "\\alpha", "\\bar\\alpha"],
    subsets := [[1,2],[1,3]]),
rec(key := "klein37", name := "Klein {3,7} graph (n=24)", n := 24, mmax := 3,
    build := function() return kleinaut(klein37_adj); end,
    symvals := [1, E(3), E(3)^2, ER(7), E(3)*ER(7), E(3)^2*ER(7),
                al, alb],
    symtexs := ["", "\\omega", "\\bar\\omega", "\\sqrt7",
                "\\omega\\sqrt7", "\\bar\\omega\\sqrt7",
                "\\alpha", "\\bar\\alpha"],
    subsets := [[1,2,4,5],[1,2,4,6],[1,3,4,5],[1,3,4,6],
                [1,7],[1,8]]),
rec(key := "klein73", name := "Klein {7,3} graph PGL(2,7)/S3 (n=56)",
    n := 56, mmax := 3,
    build := function() return kleinaut(klein73_adj); end,
    symvals := [1, ER(17), ER(2), al, alb],
    symtexs := ["", "\\sqrt{17}", "\\sqrt2", "\\alpha", "\\bar\\alpha"],
    subsets := [[1,2],[1,3],[1,4],[1,5]]),
rec(key := "zhou3", name := "Zhou-3 PSL(2,13)/S3 (n=182)", n := 182,
    mmax := 3,
    build := zhou3,
    symvals := [1, cc(1), cc(2), cc(3), ER(13)],
    symtexs := ["", "c_1", "c_2", "c_3", "\\sqrt{13}"],
    subsets := [[1,2,3],[1,2,4],[1,3,4],[1,5]]),
rec(key := "meinhold", name := "Meinhold-72 G2(2)/PSL(3,2) (n=72)",
    n := 72, mmax := 2,
    build := meinhold,
    symvals := [1, ER(-3)], symtexs := ["", "\\sqrt{-3}"],
    subsets := [[1,2]]),
];;

for case in cases do
    docase(case);;
od;

Print("fallback warnings: ", Length(FALLBACKS), "\n");;
for w in FALLBACKS do Print("  ", w[1], ": ", w[2], "\n"); od;
Print("════════ DONE ════════\n");;
