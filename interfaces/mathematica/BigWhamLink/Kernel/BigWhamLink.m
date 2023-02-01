(* ::Package:: *)

(* Wolfram Language Package *)

(* Created by Brice Lecampion 4 Dec 2020 *)
(* Copyright (c) EPFL (Ecole Polytechnique Federale de Lausanne) , Switzerland,*)
(* Geo-Energy Laboratory, 2020.  All rights reserved. *)


BeginPackage["BigWhamLink`"]

Hmat::usage = "Hmat is a symbol to which BigWhamLink messages are associated."

(* Privately load LTemplate. Note the leading ` character!! *)
(* Get["`LTemplate`LTemplatePrivate`"] *)
Get["@LTemplateLibPath@/LTemplate.m"]
Get["@LTemplateLibPath@/LTemplatePrivate.m"]

(* ConfigureLTemplate[] must be called at this point. *)
ConfigureLTemplate[
	(* You should also supply a symbol from the MyApp` context (called Hmat here) to
       associate LTemplate's standard ::error, ::warning, etc. messages with. *)
	"MessageSymbol" -> Hmat, 
	(* If lazy loading is enabled, functions are loaded only on first use.
	   This improves package loading performance, but it is not convenient
	   during development and debugging. *)
	"LazyLoading" -> False
]

Print[LClassContext[]];
AppendTo[$ContextPath,LClassContext[]];

(* Exported symbols added here with SymbolName::usage *)

$BigWhamVersion::usage = "$BigWhamVersion gives the version number of the BigWham library ";
GetFilesDate::usage = "Get dates of BigWham library and BigWhamLink library object files "

toHMatExpr::usage = "toHMatExpr[pts, conn,kernelname, properties, maxleafsize,eta,eps_aca] create the HmatExpr object \
                      holding all the hierarchical matrix information and data";
testHmatExpr::usage="testHmatExpr[id] return True or False depending if id is a HmatExpr or not ";
PlotHpattern::usage = "Plotting the hierarchical matrix pattern";
GetSize::usage = "GetSize[id] Return the hierarchical matrix size" ;
GetCollocationPoints::usage = "GetCollocationPoints[id] - return the coordinates of the collocation points"
GetCompressionRatio::usage = "GetCompressionRatio[id] output the compression ratio of the hierarchical matrix"
GetFullBlocks::usage = "GetFullBlocks[id] return the full blocks as a sparse matrix" ;
Hdot::usage = "Hdot[id,y] - return the hierarchical dot product H(id).y  (taking care of permutation)" ;
(*HdotRec::usage = "HdotRec[id,y] - return the hierarchical dot product H(id).y  (taking care of permutation) using recursion" ;*)
GetPermutation::usage = "GetPermutation[id] returns vector of permutation";
GetKernel::usage = "GetKernel[id] returns the name (as string) of the BEM kernel used for this hierarchical matrix"
GetProblemDimension::usage = "GetProblemDimension[id] returns the dimension of the vectorized BEM problem"

(*PlotHpattern2::usage = "Plotting the hierarchical matrix pattern";*)
GetPattern::usage = "GetPattern[id] return the pattern of the Hmatrix  "

ComputeStresses::usage = "ComputeStresses[id,sol,obsPts,nPts,prop,coor,conn,areDDglobal] - return the stress tensor for a list of obs pts"
ComputeDisplacements::usage = "ComputeDisplacements[id,sol,obsPts,nPts,prop,coor,conn,areDDglobal] - return the displacement vector for a list of obs pts"

Begin["`Private`"]
(* Implementation of the package *)

(* Helper function: abort package loading and leave a clean $ContextPath behind *)
packageAbort[] := (End[]; EndPackage[]; Abort[])

(*-----------------------------------------------*)

$BigWhamVersion  = "1.2.2";

$packageDirectory  = DirectoryName[$InputFileName];
$libraryDirectory  = FileNameJoin[{$packageDirectory, "../LibraryResources", $SystemID}];
$sourceDirectory   = FileNameJoin[{$packageDirectory, "../LibraryResources", "Source"}];

$buildSettingsFile = FileNameJoin[{$packageDirectory, "BuildSettings.m"}];

$buildSettings = None;
If[FileExistsQ[$buildSettingsFile], Get[$buildSettingsFile] ]

(* loading the tbb manually - seems to fix things on some OS*)
LibraryLoad[$TBBLIBPATH];

(* Add $libraryDirectory to $LibraryPath in case the package 
   is not installed in $UserBaseDirectory/Applications. 
 *)
If[Not@MemberQ[$LibraryPath, $libraryDirectory],
  AppendTo[$LibraryPath, $libraryDirectory]
];

(*-----------------------------------------------*)
(* BigWham lib onject file *)
$BigwhamLib=FileNameJoin[{$BigWhamDirectory,"/libBigWham.a"}];

(***** The library template specification*****)
(* template Specification *)
template = LTemplate["HmatExpressions",
   {
    (* Functions that are not attached to any HmatExpr (i.e. 
    could also be free) will be members of this class. 
    There will be only once instance of this class. *)
    LClass["Manager",
     {
      LFun["releaseHMatExpr", {Integer}, "Void"]
      }
     ],
    (* The matrix class. *)
    LClass["HMatExpr",
     {
      LFun["set", {{Real, 2, "Constant"}, {Integer, 2, "Constant"}, 
        "UTF8String", {Real, 1, "Constant"}, Integer, Real, Real}, 
       "Void"],
      LFun["isBuilt", {}, "Boolean"],
      LFun["getKernel", {}, "UTF8String"],
      LFun["getPermutation", {}, {Integer, 1}],
      LFun["getCollocationPoints", {}, {Real, 2}],
      LFun["getCompressionRatio", {}, Real],
      LFun["getHpattern", {}, {Integer, 2}],
(*      LFun["getHpattern2", {}, {Integer, 2}],*)
      LFun["hdot", {{Real, 1, "Constant"}}, {Real, 1}],
(*      LFun["hdotRec", {{Real, 1, "Constant"}}, {Real, 1}],*)
      LFun["getSize", {}, {Integer, 1}],
      LFun["getProblemDimension", {}, Integer],
      LFun["getFullBlocks", {}, LType[SparseArray, Real]]
(*)      ,
      LFun["computeStresses", {{Real, 1, "Constant"},{Real, 2, "Constant"},
        Integer,{Real, 1, "Constant"},{Real, 2, "Constant"},{Integer, 2, "Constant"}, 
        "Boolean"}, {Real, 2}],
        LFun["computeDisplacements", {{Real, 1, "Constant"},{Real, 2, "Constant"},
          Integer,{Real, 1, "Constant"},{Real, 2, "Constant"},{Integer, 2, "Constant"},
          "Boolean"}, {Real, 2}] *)
       (*[sol,obsPts,nPts,prop,coor,conn,areDDglobal];*)
      }
     ]
    }
   ];

(***** Compilation, loading and initialization *****)
(* include directory files *)
incdir = {
   $MKLROOT <> "/include/",
   $TBBROOT <> "/include/",
   "@CMAKE_SOURCE_DIR@",
   FileNameJoin[{"@CMAKE_SOURCE_DIR@", "/src/"}],
    FileNameJoin[{"@CMAKE_SOURCE_DIR@", "/il/"}]
   };

(* lib files directory *)
libdir = {
	$BigWhamDirectory,
  $MKLROOT <> "/lib/", $TBBROOT <> "/lib/"
  }

$linkerOptions={"-ltbb", "-ldl", "-lpthread", "-lm"};

$compilerOptions={"-std=c++17","-m64", "-DIL_MKL", "-DIL_BLAS"};

Switch[$SystemID,"MacOSX-x86-64",
$extraObjectFiles={$BigwhamLib,
    $MKLROOT <> "/lib/libmkl_intel_lp64.a",
    $MKLROOT <> "/lib/libmkl_core.a",
    $MKLROOT <> "/lib/libmkl_sequential.a"}
,"Linux-x86-64",
       $extraObjectFiles={$BigwhamLib,
    $MKLROOT <> "/lib/intel64/libmkl_intel_lp64.a",
    $MKLROOT <> "/lib/intel64/libmkl_core.a",
    $MKLROOT <> "/lib/intel64/libmkl_sequential.a"}
]
 
compileBigWhamLink[] := 
	Module[{res},
		      (* create directory for binary if it doesn't exist yet *)
      	If[Not@DirectoryQ[$libraryDirectory],
        	CreateDirectory[$libraryDirectory]
      	];
		 
		SetDirectory[$sourceDirectory];
		res = CompileTemplate[template, "IncludeDirectories" -> incdir ,
  		"LibraryDirectories" -> libdir , 
  		"ExtraObjectFiles" -> $extraObjectFiles,
   		"CompileOptions" -> $compilerOptions,
  		"LinkerOptions" -> $linkerOptions,
   		"ShellOutputFunction" -> Print, "TargetDirectory" -> $libraryDirectory];

		ResetDirectory[];
		res
	];
	
Switch[$SystemID, "Linux-x86-64",	
$HmatExpLib=FileNameJoin[{$libraryDirectory,"HmatExpressions.so"}],
"MacOSX-x86-64",
$HmatExpLib=FileNameJoin[{$libraryDirectory,"HmatExpressions.dylib"}]
]

If[FileExistsQ[$HmatExpLib],
	LoadTemplate[template];,
	compileBigWhamLink[];
	LoadTemplate[template];
]

(*------------------------------------------*)
(*------ Interface Functions ---------------*)

GetFilesDate[]:={$HmatExpLib->FileDate[$HmatExpLib],$BigwhamLib->FileDate[$BigwhamLib]};


toHMatExpr[coor_?(MatrixQ[#, NumericQ] &), 
  conn_?(MatrixQ[#, IntegerQ] &), Kernelname_?StringQ, 
  prop_?(VectorQ[#, NumericQ] &), maxleaf_?IntegerQ , eta_?NumericQ, 
  eps_?NumericQ] := Block[{vec = Make[HMatExpr]},
  vec@"set"[coor, conn-1, Kernelname, prop, maxleaf, eta, eps];  (* switch to 0 -index base for connectivity *)
  vec
  ]

testHmatExpr[h_] := (Length@Cases[LExpressionList[BigWhamLink`LTemplate`Classes`HMatExpr], h] == 1)

Options[PlotHpattern] = {Colors -> {Red, Green}};
PlotHpattern[id_?(testHmatExpr), OptionsPattern[]] :=  Module[{tt, all, col = OptionValue[Colors]}, tt = id@"getHpattern"[];
   all = {EdgeForm[Thin], If[#[[5]] == 0, col[[1]], col[[2]]], 
       Rectangle[{#[[1]], -#[[2]]}, {#[[3]], -#[[4]]}]} & /@ tt;
   Graphics[all]];

(*Options[PlotHpattern2] = {Colors -> {Red, Green}};
PlotHpattern2[id_?(testHmatExpr), OptionsPattern[]] :=  Module[{tt, all, col = OptionValue[Colors]}, tt = id@"getHpattern2"[];
all = {EdgeForm[Thin], If[#[[5]] == 0, col[[1]], col[[2]]],
  Rectangle[{#[[1]], -#[[2]]}, {#[[3]], -#[[4]]}]} & /@ tt;
Graphics[all]];
*)
GetPattern[id_?(testHmatExpr), OptionsPattern[]] :=id@"getHpattern"[];

GetCompressionRatio[id_?(testHmatExpr)]:=id@"getCompressionRatio"[];

GetPermutation[id_?(testHmatExpr)]:=id@"getPermutation"[] + 1; (* 1 index based for mathematica *)

GetSize[id_?(testHmatExpr)]:=id@"getSize"[];

GetKernel[id_?(testHmatExpr)]=id@"getKernel"[];

(* get problem vector dimension  *)
GetProblemDimension[id_?(testHmatExpr)]:=id@"getProblemDimension"[];

GetCollocationPoints[id_?(testHmatExpr)]:=id@"getCollocationPoints"[];
(*)Module[{permut,xycolPermut,dim,xycol},
  permut=id@"getPermutation"[] + 1;
  xycolPermut=id@"getCollocationPoints"[];
  dim=GetProblemDimension[id];
  xycol=ConstantArray[0.,{Length[xycolPermut],dim}];
  xycol[[permut]]=xycolPermut;
  xycol
];*) (* return in the original state *)


getPermuttedDOF[id_?(testHmatExpr)]:=Module[{dim,permut},
  dim=GetProblemDimension[id];
  permut=id@"getPermutation"[] + 1;
  Switch[dim,1,permut,
    2, Flatten[{dim*(# - 1) , dim*(# - 1) + 1} & /@ permut ] + 1,
    3,Flatten[{dim*(# - 1) , dim*(# - 1) + 1, dim*(# - 1) + 2} & /@ permut ] + 1,_,0]
];

GetFullBlocks[id_?(testHmatExpr)]:=id@"getFullBlocks"[];
(*    Module[{pdof,fullb},*)
(*  pdof=getPermuttedDOF[id];*)
(*  fullb=SparseArray[{{1,1}-> 0.},GetSize[id]];*)
(*  fullb[[pdof,pdof]]=id@"getFullBlocks"[];*)
(*  fullb*)
(*]; *)(* return in the original state *)

(*Hdot::error = "";*)
Hdot::nnarg = "The argument `1` has not the proper length"
Hdot[id_?(testHmatExpr),x_?(VectorQ[#, NumericQ] &)]:= If[Length[x]==(id@"getSize"[])[[2]], id@"hdot"[x],
Message[Hdot::nnarg,x];Abort[]];

(*HdotRec::nnarg = "The argument `1` has not the proper length"
HdotRec[id_?(testHmatExpr),x_?(VectorQ[#, NumericQ] &)]:= If[Length[x]==(id@"getSize"[])[[2]], id@"hdotRec"[x],
  Message[HdotRec::nnarg,x];Abort[]];*)


ComputeStresses[id_?(testHmatExpr), sol_?(VectorQ[#, NumericQ] &), 
obsPts_?(MatrixQ[#, NumericQ] &), nPts_?IntegerQ, prop_?(VectorQ[#, NumericQ] &),
coor_?(MatrixQ[#, NumericQ] &), conn_?(MatrixQ[#, IntegerQ] &),areDDglobal_?BooleanQ]:=id@"computeStresses"[sol,obsPts,nPts,prop,coor,conn-1,areDDglobal]; (* switch to 0 -index base for connectivity *)

ComputeDisplacements[id_?(testHmatExpr), sol_?(VectorQ[#, NumericQ] &),
             obsPts_?(MatrixQ[#, NumericQ] &), nPts_?IntegerQ, prop_?(VectorQ[#, NumericQ] &),
             coor_?(MatrixQ[#, NumericQ] &), conn_?(MatrixQ[#, IntegerQ] &),areDDglobal_?BooleanQ]:=id@"computeDisplacements"[sol,obsPts,nPts,prop,coor,conn-1,areDDglobal]; (* switch to 0 -index base for connectivity *)


End[]

EndPackage[]

