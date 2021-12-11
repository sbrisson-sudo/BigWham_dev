BeginTestSection["GriffithCrack-2D"]

BeginTestSection["S3DP0"]

VerificationTest[(* 1 *)
	CompoundExpression[Set[Nn, 70], Set[pts, Table[List[Plus[-1, Times[2.`, Times[i, Power[Nn, -1]]]], 0.`], List[i, 0, Nn]]], Set[conn, Table[List[i, Plus[i, 1]], List[i, 1, Plus[Length[pts], -1]]]], Set[YME, 1.`], Set[Nu, 0.2`], Set[inst, toHMatExpr[pts, conn, "S3DP0", List[YME, Nu, 1000.`], 16, 30.`, 0.001`]], Set[xy, GetCollocationPoints[inst]], Set[ww, Map[Function[Times[Times[4, Power[Times[YME, Power[Plus[1, Times[-1, Power[Nu, 2]]], -1]], -1]], Sqrt[Plus[1, Times[-1, Power[Part[Slot[1], 1], 2]]]]]], xy]], Set[xx, ConstantArray[0.`, Times[2, Length[xy]]]], Set[Part[xx, Span[2, -1, 2]], ww], Set[f, Hdot[inst, xx]], Set[test, Less[Plus[1, Times[-1, Median[Part[f, Span[2, -1, 2]]]]], 0.0013`]], Unset[inst], test]
	,
	True	
	,
	TestID-> "S3DP0-Griffith-Hdot"
]

VerificationTest[(* 2 *)
	CompoundExpression[Set[Nn, 70], Set[pts, Table[List[Plus[-1, Times[2.`, Times[i, Power[Nn, -1]]]], 0.`], List[i, 0, Nn]]], Set[conn, Table[List[i, Plus[i, 1]], List[i, 1, Plus[Length[pts], -1]]]], Set[YME, 1.`], Set[Nu, 0.2`], Set[inst, toHMatExpr[pts, conn, "S3DP0", List[YME, Nu, 1000.`], 16, 30.`, 0.0001`]], Set[xy, GetCollocationPoints[inst]], Set[ww, Map[Function[Times[Times[4, Power[Times[YME, Power[Plus[1, Times[-1, Power[Nu, 2]]], -1]], -1]], Sqrt[Plus[1, Times[-1, Power[Part[Slot[1], 1], 2]]]]]], xy]], Set[feb, GetFullBlocks[inst]], Set[F, ConstantArray[0.`, Times[2, Length[xy]]]], Set[Part[F, Span[2, -1, 2]], 1.`], SetDelayed[fdot[Pattern[x, Blank[]]], Hdot[inst, x]], Set[dataP, SparseArray`SparseMatrixILU[feb, Rule[Method, "ILUT"], Rule["FillIn", 5], Rule["Tolerance", Power[10, -5.`]]]], Set[prec, Function[x, SparseArray`SparseMatrixApplyILU[dataP, x]]], Set[List[sol, stats], Reap[LinearSolve[fdot, F, Rule[Method, List["Krylov", List[Rule["Method", "BiCGSTAB"], Rule["Preconditioner", prec], Rule["PreconditionerSide", Right], Rule["Tolerance", Power[10, -7.`]], Rule["MaxIterations", 1000], Rule["ResidualNormFunction", Function[Sow[Norm[Slot[1], 2]]]]]]]]]], Set[w, Part[sol, Span[2, -1, 2]]], Set[relErrr, ReplaceAll[Times[Abs[Plus[w, Times[-1, ww]]], Power[ww, -1]], List[Rule[ComplexInfinity, 0]]]], Set[test, Less[Median[relErrr], Power[10, -2.`]]], Unset[inst], test]
	,
	True	
	,
	TestID->"S3DP0-Griffith-HIterativeSolve"
]

EndTestSection[]

BeginTestSection["2DP1"]

VerificationTest[(* 3 *)
	CompoundExpression[Set[Nn, 70], Set[pts, Table[List[Plus[-1, Times[2.`, Times[i, Power[Nn, -1]]]], 0.`], List[i, 0, Nn]]], Set[conn, Table[List[i, Plus[i, 1]], List[i, 1, Plus[Length[pts], -1]]]], Set[YME, 1], Set[Nu, 0.2`], Set[inst, toHMatExpr[pts, conn, "2DP1", List[YME, Nu], 16, 30.`, 0.0001`]], Set[xy, GetCollocationPoints[inst]], Set[xyn, Map[Function[Part[pts, Slot[1]]], Flatten[Map[Function[Part[conn, Slot[1]]], Range[Nn]]]]], Set[ww, Map[Function[Times[Times[4, Power[Times[YME, Power[Plus[1, Times[-1, Power[Nu, 2]]], -1]], -1]], Sqrt[Plus[1, Times[-1, Power[Part[Slot[1], 1], 2]]]]]], xyn]], Set[xx, ConstantArray[0.`, Times[2, Length[xyn]]]], Set[Part[xx, Span[2, -1, 2]], ww], Set[f, Hdot[inst, xx]], Set[test, Less[Median[Abs[Plus[1, Times[-1, Part[f, Span[2, -1, 2]]]]]], 0.004`]]]
	,
	True	
	,
	TestID->"2DP1-Griffith-Hdot"
]

VerificationTest[(* 4 *)
	CompoundExpression[Set[Nn, 70], Set[pts, Table[List[Plus[-1, Times[2.`, Times[i, Power[Nn, -1]]]], 0.`], List[i, 0, Nn]]], Set[conn, Table[List[i, Plus[i, 1]], List[i, 1, Plus[Length[pts], -1]]]], Set[YME, 1], Set[Nu, 0.2`], Set[inst, toHMatExpr[pts, conn, "2DP1", List[YME, Nu], 16, 30.`, 0.0001`]], Set[xyn, Map[Function[Part[pts, Slot[1]]], Flatten[Map[Function[Part[conn, Slot[1]]], Range[Nn]]]]], Set[ww, Map[Function[Times[Times[4, Power[Times[YME, Power[Plus[1, Times[-1, Power[Nu, 2]]], -1]], -1]], Sqrt[Plus[1, Times[-1, Power[Part[Slot[1], 1], 2]]]]]], xyn]], Set[feb, GetFullBlocks[inst]], Set[F, ConstantArray[0.`, Part[GetSize[inst], 1]]], Set[Part[F, Span[2, -1, 2]], 1.`], SetDelayed[fdot[Pattern[x, Blank[]]], Hdot[inst, x]], Set[dataP, SparseArray`SparseMatrixILU[feb, Rule[Method, "ILUT"], Rule["FillIn", 5], Rule["Tolerance", Power[10, -5.`]]]], Set[prec, Function[x, SparseArray`SparseMatrixApplyILU[dataP, x]]], Set[List[sol, stats], Reap[LinearSolve[fdot, F, Rule[Method, List["Krylov", List[Rule["Method", "BiCGSTAB"], Rule["Preconditioner", prec], Rule["PreconditionerSide", Right], Rule["Tolerance", Power[10, -7.`]], Rule["MaxIterations", 1000], Rule["ResidualNormFunction", Function[Sow[Norm[Slot[1], 2]]]]]]]]]], Set[relErrr, Times[Abs[Plus[Part[sol, Span[4, -3, 2]], Times[-1, Part[ww, Span[2, -2]]]]], Power[Part[ww, Span[2, -2]], -1]]], Less[Median[relErrr], 0.00224`]]
	,
	True	
	,
	TestID->"2DP1-Griffith-HIterativeSolve"
]

EndTestSection[]

EndTestSection[]
