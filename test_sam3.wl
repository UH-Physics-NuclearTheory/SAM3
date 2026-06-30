(* ===================================================================
   Self-contained tests for the SAM-3.0 Mathematica notebook (SAM-3.0.nb).
   Wolfram counterpart of cpp/test_sam3.cpp.

     wolframscript -file test_sam3.wl

   Part A: KappaCE vs analytic ground truth (Gaussian, exact binomial,
           two observables, multi-charge Schur), evaluated with BOTH notebook
           methods ("ColoredPartitions" and "CoefficientMatching").
   Part B: the two independent methods must agree on a generic non-Gaussian,
           multi-charge input.
   Exits with code 1 if any check fails.
   =================================================================== *)

scriptDir = Which[
  $InputFileName =!= "" && $InputFileName =!= Null,        DirectoryName[$InputFileName],
  TrueQ[$Notebooks] && Quiet[NotebookDirectory[]] =!= $Failed, NotebookDirectory[],
  True,                                                    Directory[]];
nbPath = FileNameJoin[{scriptDir, "SAM-3.0.nb"}];
If[! FileExistsQ[nbPath], Print["cannot find ", nbPath]; Exit[2]];

nb = Get[nbPath];
codecells = Cases[nb, Cell[BoxData[b_], "Code", ___] :> b, Infinity];
Quiet[Scan[ToExpression, codecells]];
If[! ListQ[Options[KappaCE]], Print["FAILED to load KappaCE from ", nbPath]; Exit[1]];

relTol = 1.*^-9; absTol = 1.*^-9; fails = 0;
close[g_, e_] := Abs[g - e] <= absTol + relTol*Max[1, Abs[g], Abs[e]];
check[label_, got_, exp_] := If[close[N@got, N@exp],
  Print["  PASS  ", label, " = ", N@got],
  Print["  FAIL  ", label, "  got=", N@got, "  expected=", N@exp,
        "  diff=", N[got - exp]]; fails++];
chargeLabels[s_] := Take[{"B", "Q", "S", "C", "D", "E"}, s];
ce[idx_, s_, meth_] :=
  KappaCE[idx, chargeLabels[s], "Form" -> "Raw", "Simplify" -> Identity, "Method" -> meth];
(* check a cumulant against an analytic value via BOTH notebook methods *)
checkBoth[label_, idx_, s_, exp_] := (
  check[label <> " [ColoredPartitions]",   ce[idx, s, "ColoredPartitions"],   exp];
  check[label <> " [CoefficientMatching]", ce[idx, s, "CoefficientMatching"], exp];
);

(* ---------- Part A: analytic ground truth (mirrors cpp/test_sam3.cpp) ---------- *)
Print["== Part A: KappaCE vs analytic ground truth (both methods) =="];

ClearAll[kGCE]; kGCE[{2},{0}]=10; kGCE[{1},{1}]=3; kGCE[{0},{2}]=5; kGCE[_List,_List]:=0;
ResetKappaCache[];
checkBoth["Gaussian C2", {2}, 1, 41/5];           (* 10 - 9/5 *)

ClearAll[kGCE]; mu=100; al=3/10; kGCE[{n_},{m_}] := If[n>0, al*mu, If[m>0, mu, 0]];
ResetKappaCache[];
checkBoth["Binomial C1", {1}, 1, 30];
checkBoth["Binomial C2", {2}, 1, 21];
checkBoth["Binomial C3", {3}, 1, 42/5];
checkBoth["Binomial C4", {4}, 1, -546/100];

ClearAll[kGCE];
kGCE[{2,0},{0}]=10; kGCE[{0,2},{0}]=20; kGCE[{1,1},{0}]=4;
kGCE[{1,0},{1}]=3; kGCE[{0,1},{1}]=-2; kGCE[{0,0},{2}]=5; kGCE[_List,_List]:=0;
ResetKappaCache[];
checkBoth["TwoObs C20", {2,0}, 1, 41/5];
checkBoth["TwoObs C02", {0,2}, 1, 96/5];
checkBoth["TwoObs C11", {1,1}, 1, 26/5];

ClearAll[kGCE];
kGCE[{2},{0,0}]=12; kGCE[{1},{1,0}]=3; kGCE[{1},{0,1}]=-2;
kGCE[{0},{2,0}]=5; kGCE[{0},{1,1}]=1; kGCE[{0},{0,2}]=4; kGCE[_List,_List]:=0;
ResetKappaCache[];
With[{KK = {{5,1},{1,4}}, CC = {3,-2}},
  checkBoth["MultiCharge(s=2) C2", {2}, 2, 12 - CC.Inverse[KK].CC]];

(* ---------- Part B: the two methods must agree (generic, non-Gaussian) ---------- *)
Print["== Part B: ColoredPartitions vs CoefficientMatching (generic s=2, order 3) =="];
ClearAll[kGCE];
kGCE[{1},{0,0}]=0;   kGCE[{2},{0,0}]=8;    kGCE[{3},{0,0}]=2;
kGCE[{1},{1,0}]=3/2; kGCE[{1},{0,1}]=-1;   kGCE[{2},{1,0}]=7/10;
kGCE[{2},{0,1}]=-3/10; kGCE[{1},{2,0}]=1/10; kGCE[{1},{0,2}]=-1/5; kGCE[{1},{1,1}]=3/20;
kGCE[{0},{2,0}]=5;   kGCE[{0},{0,2}]=4;    kGCE[{0},{1,1}]=1;
kGCE[{0},{3,0}]=1/2; kGCE[{0},{0,3}]=-2/5; kGCE[{0},{2,1}]=1/5; kGCE[{0},{1,2}]=3/10;
kGCE[_List,_List]:=0;
ResetKappaCache[];
Do[
  check["Generic s=2 C" <> ToString[k] <> " (methods agree)",
        ce[{k}, 2, "ColoredPartitions"], ce[{k}, 2, "CoefficientMatching"]],
  {k, 1, 3}];

Print[""];
If[fails == 0,
  Print["All SAM-3.0 Mathematica tests passed."]; Exit[0],
  Print[fails, " test(s) FAILED."]; Exit[1]];
