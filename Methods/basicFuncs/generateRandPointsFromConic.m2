restart
needsPackage "NumericalAlgebraicGeometry"
generateRandPointsFromConic = (M) -> (
    xDeno = M_(1,1);
    yDeno = M_(0,0);
    if xDeno > yDeno then(
	R := CC[y];
	k := -M_(1,2);
	b := sqrt yDeno;
	per := random(1,100)/100;
	x := k+per*b;
    	uT := matrix{{x,y,1}};
	u := transpose uT;
	eqn := flatten entries(uT*M*u);
	sol := flatten solveSystem(eqn);
	return transpose(matrix{{x,first(flatten(entries(matrix(sol_0)))),1}});
    )else(
	F := CC[xNew];
	hNew := -M_(0,2);
	bNew := sqrt xDeno;
	perNew := random(1,100)/100;
	yNew := hNew+perNew*bNew;
	uTNew := matrix{{xNew,yNew,1}};
	uNew := transpose uTNew;
	eqnNew := flatten entries(uTNew*M*uNew);
	solNew := flatten solveSystem(eqnNew);
	return transpose(matrix{{first(flatten(entries(matrix(solNew_0)))),yNew,1}});
	);
    )

generateRandRealPointsFromConic = (M) -> (
    xDeno = M_(1,1);
    yDeno = M_(0,0);
    if xDeno > yDeno then(
	R := QQ[y];
	k := -M_(1,2);
	b := realPart(sqrt(yDeno));
	per := random(1,100)/100;
	x := k+per*b;
    	uT := matrix{{x_QQ,y,1}};
	u := transpose uT;
	eqn := flatten entries(uT*M*u);
	sol := flatten solveSystem(eqn);
	return transpose(matrix{{x,first(flatten(entries(matrix(sol_0)))),1}});
    )else(
	F := QQ[xNew];
	hNew := -M_(0,2);
	bNew := realPart(sqrt(xDeno));
	perNew := random(1,100)/100;
	yNew := hNew+perNew*bNew;
	uTNew := matrix{{xNew,yNew_QQ,1}};
	uNew := transpose uTNew;
	eqnNew := flatten entries(uTNew*M*uNew);
	solNew := flatten solveSystem(eqnNew);
	return transpose(matrix{{first(flatten(entries(matrix(solNew_0)))),yNew,1}});
	);
    )
