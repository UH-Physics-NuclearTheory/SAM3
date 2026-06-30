(* ===================================================================
   Cross-check the C++ implementation (cpp/sam3.h) against the Mathematica
   notebook (SAM-3.0.nb) on identical random inputs.

     wolframscript -file crosscheck_cpp.wl <bridge1.m> [bridge2.m ...]

   Each bridge*.m is produced by cpp/bridge_sam3 and defines dim, chg, nmax,
   the kGCE downvalues, and ceCpp = {{n, value}, ...}. This script loads the
   notebook's KappaCE, evaluates it on the same kGCE, and compares to ceCpp.
   Exits with code 1 if any cumulant disagrees.
   =================================================================== *)

scriptDir = Which[
  $InputFileName =!= "" && $InputFileName =!= Null,        DirectoryName[$InputFileName],
  TrueQ[$Notebooks] && Quiet[NotebookDirectory[]] =!= $Failed, NotebookDirectory[],
  True,                                                    Directory[]];
nbPath = FileNameJoin[{scriptDir, "SAM-3.0.nb"}];
bridgeFiles = Rest[$ScriptCommandLine];
If[bridgeFiles === {}, Print["usage: ... <bridge.m> [bridge.m ...]"]; Exit[2]];
If[! FileExistsQ[nbPath], Print["cannot find ", nbPath]; Exit[2]];

nb = Get[nbPath];
codecells = Cases[nb, Cell[BoxData[b_], "Code", ___] :> b, Infinity];
Quiet[Scan[ToExpression, codecells]];
If[! ListQ[Options[KappaCE]], Print["FAILED to load KappaCE from ", nbPath]; Exit[1]];

relTol = 1.*^-9; absTol = 1.*^-9; fails = 0;
close[g_, e_] := Abs[g - e] <= absTol + relTol*Max[1, Abs[g], Abs[e]];
chargeLabels[s_] := Take[{"B", "Q", "S", "C", "D", "E"}, s];

Do[
  Get[bf];
  ResetKappaCache[];
  Print[FileNameTake[bf], "  (d=", dim, " s=", chg, " Nmax=", nmax,
        ", ", Length[ceCpp], " cumulants)"];
  maxDiff = 0.;
  Do[
    With[{idx = pair[[1]], cpp = pair[[2]]},
      mma = N@KappaCE[idx, chargeLabels[chg], "Form" -> "Raw", "Simplify" -> Identity];
      d2 = Abs[mma - cpp];
      maxDiff = Max[maxDiff, d2];
      If[! close[mma, cpp],
        Print["    FAIL  C", idx, "  mma=", mma, "  cpp=", cpp, "  diff=", d2];
        fails++]],
    {pair, ceCpp}];
  Print["    max |MMA - C++| = ", maxDiff,
        If[close[maxDiff, 0.], "   OK", "   <-- MISMATCH"]],
  {bf, bridgeFiles}];

Print[""];
If[fails == 0,
  Print["C++ and Mathematica agree on all cumulants."]; Exit[0],
  Print[fails, " cumulant(s) DISAGREE."]; Exit[1]];
