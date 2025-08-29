(* ::Package:: *)

f="tocheck_P7P8.in"
(*Print[$CommandLine];
init = ToExpression[$CommandLine[[4]]];*)
init=0
While [init<2673,
init=init+1; (*start at 1, end at 2674*)
final= init;
(*Print[Head @ init];*)
Print[init];
Print[final];
timec = 600;
outfname = StringJoin["./results/", "checked_", ToString[init], "_", ToString[final], ".out"];
(*tempdir1="//"
tempdir=StringJoin[tempdir1,"tmp_",ToString[init],"_",ToString[final]];
CreateDirectory[tempdir];
$TemporaryDirectory=tempdir;*)
outf = OpenWrite[outfname];
WriteString[outf,StringJoin["#aborting after ", ToString[timec]],"\n"];
file = OpenRead[f];
line = ReadLine[file];
i = 1;
While[Not[SameQ[line, EndOfFile]] && i <= final,
 If[i >= init && i <= final,
  WriteString["stdout", i, ";"];
  coefexpr = StringSplit[line, ";"];
  expr = ToExpression[coefexpr[[2]]];
  WriteString[outf, i, ";"];
  WriteString[outf, coefexpr[[2]], ";"];
  out1 = With[{assum = alpha > 1 && 0 < beta < 1 && 0 < omega < 1 && 0 < gamma <= 1}, 
    TimeConstrained[Simplify[Reduce[assum \[Implies] expr >= 0, Reals], assum], timec]];
  WriteString[outf, out1, "\n"];
  Close[outf];, 1];
 line = ReadLine[file];
 i = i + 1;
 ]
Close[file];
]



