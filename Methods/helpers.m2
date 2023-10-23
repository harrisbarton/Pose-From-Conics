restart
needsPackage "NumericalAlgebraicGeometry";

------------------------------ Small functions that will be used ---------------------------------
-- Cross Product
cross = method(Vector,Vector) := (u,v) -> (
    f := new MutableList from {3:0};
    f#0 = u_(1)*v_(2) - u_(2)*v_(1);
    f#1 = u_(2)*v_(0) - u_(0)*v_(2);
    f#2 = u_(0)*v_(1) - u_(1)*v_(0);
    out := vector(toList f);
    return out
    )

-- Define a plane through 3 points
plane = method(Vector, Vector, Vector) := (u,v,w) -> (
    coordOneMat = matrix{{u_1, v_1, w_1}, {u_2, v_2, w_2}, {u_3, v_3, w_3}};
    coordTwoMat = matrix{{u_0, v_0, w_0}, {u_2, v_2, w_2}, {u_3, v_3, w_3}};
    coordThreeMat = matrix{{u_0, v_0, w_0}, {u_1, v_1, w_1}, {u_3, v_3, w_3}};
    coordFourMat = matrix{{u_0, v_0, w_0}, {u_1, v_1, w_1}, {u_2, v_2, w_2}};
    return (det coordOneMat, (-1)*(det coordTwoMat), det coordThreeMat, (-1)*(det coordFourMat))
    )
        
-- Pointwise Multiplication
ptwsProd = method(Vector,Vector) := (u,v) -> (
    return diagonalMatrix(entries(u))*v
    )

-- Dot Product
dot = method(Vector,Vector) := (u,v) -> (
    return sum(entries(diagonalMatrix(entries(u))*v))
    )

-- Scalar Multiplication
scalarProd = method(RingElement,Vector) := (c,v) -> (
    n := #entries(v);
    u := vector({n:c});
    return ptwsProd(u,v)
    )

-- Pointwise divide
ptwsDiv = method(Vector,Vector) := (f,g) ->(
    n = #entries(f);
    h = new MutableList from {n:0};
    for i from 0 to n-1 do(
	    h#i = f_(i)/g_(i);
	  );
    return h;
)

-- Return the diagonal entries of a square matrix
diagEntries = method(Matrix) := (M) -> (
    n = rank source M;
    f := new MutableList from {n:0};
    for i from 0 to n-1 do (
	f#i = M_(i,i);
	);
    return toList f;
    )

reconstructConic = (u,v) -> (
    R := QQ;
    x := sub(u,R);
    y := sub(v,R);
    ones := sub(transpose matrix{{-1,-1,-1,-1,-1}},R);
    --conicCoeffs = transpose matrix{{a,b,c,d,e}};
    xsq := ptwsProd(x,x);
    xy := ptwsProd(x,y);
    ysq := ptwsProd(y,y);
    sourceMat := transpose matrix{entries xsq, entries xy, entries ysq, entries x, entries y};
    sol = flatten entries (inverse(sourceMat)*ones);
    --polys = entries (sourceMat*conicCoeffs-ones);
    --sol := entries solve(sourceMat,ones);
    a := sol_0;
    b := sol_1;
    c := sol_2;
    d := sol_3;
    e := sol_4;
    output = matrix{{a,b/2,d/2},{b/2,c,e/2},{d/2,e/2,1}}
)

-- Turn points in R^3 into points in P^2
zoomIn = (pts) -> (
    n := rank source pts;
    rowOf1s := mutableMatrix(ZZ,1,n);
    for i from 0 to n-1 do (rowOf1s_(0,i) = 1);
    rowOfOnes := matrix rowOf1s;
    ptsT := transpose pts;
    c := ptsT_2;
    firstRow := ptwsDiv(ptsT_0,c);
    secondRow := ptwsDiv(ptsT_1,c);
    ptsZoomIn := matrix{toList firstRow, toList secondRow} || rowOfOnes
    )

-- Finding Tangency Points
findTangencyPts = (Cw) -> (
    R := QQ[x,y];
    Aw := sub(Cw,R);
    u := transpose matrix{{x,y,1}};
    uT := transpose u;
    u0 := transpose matrix{{x,y}};
    u0T := transpose u0;
    -- Solve for the tangency points
    qw := first flatten entries(uT*Aw*u);
    Jacw := matrix(diff(u0,qw));
    JacEq := first flatten entries(u0T*Jacw);
    solw1 := flatten entries matrix first solveSystem({qw,JacEq});
    tangencyPt1 := transpose matrix{{solw1_0,solw1_1,1}};
    solw2 := flatten entries matrix ((solveSystem({qw,JacEq}))_1);
    tangencyPt2 := transpose matrix{{solw2_0,solw2_1,1}};
    return tangencyPt1 | tangencyPt2;
    )

-- Generate Random Points from Conics (with complex numbers)
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

-- generate Random Points from Conics (only integers)
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

end
