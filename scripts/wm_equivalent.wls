ClearAll["Global`*"];

g[w_, t_, signed_] := 
 With[{n = Length@w}, 
  Graph@Join[
    Flatten@Table[
      pr[t, i] \[DirectedEdge] r[t, i, s], {i, n}, {s, {-1, 1}}], 
    Flatten@Table[
      c[t, j, s] \[DirectedEdge] pc[t, j], {j, n}, {s, {-1, 1}}], 
    Flatten@Table[
      If[w[[i, j]] == 0, Nothing, 
       r[t, i, s] \[DirectedEdge] c[t, j, s w[[i, j]]]], {i, n}, {j, 
       n}, {s, {-1, 1}}], 
    If[signed, {}, 
     Join[Table[rm[t, i] \[DirectedEdge] r[t, i, 1], {i, n}], 
      Table[c[t, j, 1] \[DirectedEdge] cm[t, j], {j, n}]]]]]
WMEquivalent[a_, b_, signed_ : True] := 
 Module[{n = Length@a, f, ri, ci, p, q},
  f = FindGraphIsomorphism[g[a, 0, signed], g[b, 1, signed], 1];
  If[f === {}, Return["DISTINCT"]];
  f = First@f;
  ri = Table[f[r[0, i, 1]] /. r[1, j_, s_] :> {j, s}, {i, n}];
  ci = Table[f[c[0, i, 1]] /. c[1, j_, s_] :> {j, s}, {i, n}];
  p = SparseArray[Table[{ri[[i, 1]], i} -> ri[[i, 2]], {i, n}], {n, n}];
  q = SparseArray[Table[{ci[[i, 1]], i} -> ci[[i, 2]], {i, n}], {n, n}];
  {p, q}]
