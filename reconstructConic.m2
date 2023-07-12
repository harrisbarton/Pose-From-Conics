restart
needsPackage "NumericalAlgebraicGeometry";
load "smallFuncs.m2";

--R = CC;
--UimPart = random(R^3,R^5);
--Uim = UimPart||matrix{{1,1,1,1,1}};

--x = (random(R^5,R^1))_0;
--y = (random(R^5,R^1))_0;

reconstructConic = (x,y) -> (
    R := CC;
    x := sub(x,R);
    y := sub(y,R);
    ones := sub(transpose matrix{{-1,-1,-1,-1,-1}},R);
    --conicCoeffs = transpose matrix{{a,b,c,d,e}};
    xsq := ptwsProd(x,x);
    xy := ptwsProd(x,y);
    ysq := ptwsProd(y,y);
    sourceMat := transpose matrix{entries xsq, entries xy, entries ysq, entries x, entries y};
    --sol = inverse(sourceMat)*ones
    --polys = entries (sourceMat*conicCoeffs-ones);
    sol := solve(sourceMat,ones)
    a := sol_0;
    b := sol_1;
    c := sol_2;
    d := sol_3;
    e := sol_4;
    output = matrix{{a,b/2,c/2},{b/2,c,e/2},{d/2,e/2,1}}
)
