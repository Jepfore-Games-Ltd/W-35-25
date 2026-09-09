ClearAll["Global`*"];
swizzled[v_] := Table[v[[Mod[i + j - 2, 5] + 1]], {i, 5}, {j, 5}];
border[core_, edge_, corner_] := ArrayFlatten[{
    {core, ConstantArray[edge, {5, 1}]},
    {ConstantArray[edge, {1, 5}], {{corner}}}
}];
a = {1, 1, -1, 1, -1};
b = {-1, 1, 1, 0, 0};
c = {1, -1, 0, 0, 1};
a0 = border[swizzled[a], 0, 1];
a1 = border[swizzled[b], 1, 1];
a2 = border[swizzled[c], -1, 1];
f0 = {{0, 1, -1, -1, 1, 0}, 
      {1, 0, 1, -1, -1, 0},
      {-1, 1, 0, 1, -1, 0},
      {-1, -1, 1, 0, 1, 0},
      {1, -1, -1, 1, 0, 0}};
f = f0[[All, {1, 2, 4, 5, 3, 6}]];
S = RotateRight[IdentityMatrix[5]];
moving = KroneckerProduct[a0, IdentityMatrix[5]] +
         KroneckerProduct[a1, S + Transpose[S]] +
         KroneckerProduct[a2, MatrixPower[S, 2] + Transpose[MatrixPower[S, 2]]];
rolling = KroneckerProduct[f, ConstantArray[1, {1, 5}]];
w35 = ArrayFlatten[{
    {moving, Transpose[rolling]},
    {rolling, ConstantArray[1, {5, 5}]}
}];

On[Assert];
Assert[SymmetricMatrixQ[w35]];
Assert[Total[Abs@w35] == ConstantArray[25, 35]];
Assert[w35 . Transpose[w35] == 25 IdentityMatrix[35]];
MatrixPlot[w35, Mesh -> {6, 6}, MeshStyle -> Black]
