# Usage: gap -q -b enumerate.g with N, IDX, ORI, OUT set via -c "N:=9;; IDX:=36;; ORI:=false;; OUT:=\"file\";;"
# ORI=false: torsion-free subgroups of index IDX in the extended triangle group (NEC surface groups,
#            Fuchsian ones of that index in Delta^pm cannot occur by parity).
# ORI=true : torsion-free subgroups of index IDX in the Fuchsian triangle group Delta(2,3,N) = <x,y|x^2,y^3,(xy)^N>,
#            then transported to Delta^pm (x=ab, y=bc) and the action on Delta^pm/K is exported (2*IDX points).
F := FreeGroup("a","b","c");;
a := F.1;; b := F.2;; c := F.3;;
G := F / [a^2, b^2, c^2, (a*b)^2, (b*c)^3, (c*a)^N];;
a := G.1;; b := G.2;; c := G.3;;
excl := [a, b, c, a*b, b*c];;
Append(excl, List([1..N-1], j -> (c*a)^j));;
t0 := Runtime();;
if ORI = false then
  subs := LowIndexSubgroupsFpGroup(G, TrivialSubgroup(G), IDX, List(excl, UnderlyingElement));;
  subs := Filtered(subs, K -> Index(G, K) = IDX);;
  # torsion-free double check: no excluded element fixes a coset
  subs := Filtered(subs, function(K) local hom; hom := FactorCosetAction(G, K);
      return ForAll(excl, w -> ForAll([1..IDX], p -> p^Image(hom, w) <> p)); end);;
else
  F2 := FreeGroup("x","y");;
  x := F2.1;; y := F2.2;;
  Gf := F2 / [x^2, y^3, (x*y)^N];;
  x := Gf.1;; y := Gf.2;;
  excl2 := [x, y];;
  Append(excl2, List([1..N-1], j -> (x*y)^j));;
  subs0 := LowIndexSubgroupsFpGroup(Gf, TrivialSubgroup(Gf), IDX, List(excl2, UnderlyingElement));;
  subs0 := Filtered(subs0, K -> Index(Gf, K) = IDX);;
  subs0 := Filtered(subs0, function(K) local hom; hom := FactorCosetAction(Gf, K);
      return ForAll(excl2, w -> ForAll([1..IDX], p -> p^Image(hom, w) <> p)); end);;
  subs := List(subs0, K -> Subgroup(G, List(GeneratorsOfGroup(K), g -> MappedWord(UnderlyingElement(g), [F2.1, F2.2], [a*b, b*c]))));;
fi;
Print("N=", N, " IDX=", IDX, " ORI=", ORI, " classes=", Length(subs), " time=", Runtime()-t0, "\n");
out := OutputTextFile(OUT, false);;
SetPrintFormattingStatus(out, false);
for K in subs do
  hom := FactorCosetAction(G, K);
  d := Index(G, K);
  AppendTo(out, "{\"N\":", N, ",\"ori\":", ORI, ",\"d\":", d,
     ",\"a\":", ListPerm(Image(hom, a), d), ",\"b\":", ListPerm(Image(hom, b), d), ",\"c\":", ListPerm(Image(hom, c), d),
     ",\"gens\":", List(GeneratorsOfGroup(K), g -> String(g)), "}\n");
od;
CloseStream(out);
QUIT;
