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

generateRandIntPointsFromConic = (M) -> (
    xDeno = M_(1,1);
    yDeno = M_(0,0);
    if xDeno > yDeno then(
	R := ZZ[y];
	k := sub(-M_(1,2),ZZ);
	b := round realPart(sqrt(yDeno));
	x := random(k-b,k+b);
    	uT := matrix{{x,y,1}};
	u := transpose uT;
	eqn := flatten entries(uT*sub(M,ZZ)*u);
	sol := flatten solveSystem(eqn);
	return transpose(matrix{{x,first(flatten(entries(matrix(sol_0)))),1}});
    )else(
	F := ZZ[xNew];
	hNew := sub(-M_(0,2),ZZ);
	bNew := round realPart(sqrt(xDeno));
	yNew := random(hNew-bNew,hNew+bNew);
	uTNew := matrix{{xNew,yNew,1}};
	uNew := transpose uTNew;
	eqnNew := flatten entries(uTNew*sub(M,ZZ)*uNew);
	solNew := flatten solveSystem(eqnNew);
	return transpose(matrix{{first(flatten(entries(matrix(solNew_0)))),yNew,1}});
	);
    )
