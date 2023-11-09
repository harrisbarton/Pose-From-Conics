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

-- Pseudoinverse
pseudoInverse = method(Matrix) := (M) -> (
    a := M_(0,0); b := M_(0,1); c := M_(1,1); d := M_(0,2); e := M_(1,2); f := M_(2,2);
    e11 = det matrix{{c,e},{e,f}}; 
    e12 = det matrix{{b,e},{d,f}};
    e13 = det matrix{{b,c},{d,e}};
    e21 = det matrix{{b,d},{e,f}};
    e22 = det matrix{{a,d},{d,f}};
    e23 = det matrix{{a,b},{d,e}};
    e31 = det matrix{{b,d},{c,e}};
    e32 = det matrix{{a,d},{b,e}};
    e33 = det matrix{{a,b},{b,c}};
    return matrix{{e11,e12,e13},{e21,e22,e23},{e31,e32,e33}};
    
)
---------------------------------------------------------------------------------------------------

-- Reconstruct conic matrix given five points
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

-- Turn R^3 points into P^2 points
zoomIn = (pts) -> (
    n := rank source pts;
    rowOf1s := mutableMatrix(ZZ,1,n);
    for i from 0 to n-1 do (rowOf1s_(0,i) = 1);
    rowOfOnes := matrix rowOf1s;
    ptsT := transpose pts;
    c := ptsT_2;
    firstRow := ptwsDiv(ptsT_0,c);
    secondRow := ptwsDiv(ptsT_1,c);
    ptsZoomedIn := matrix{toList firstRow, toList secondRow} || rowOfOnes
    )