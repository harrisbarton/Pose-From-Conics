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
    ones := sub(transpose matrix{{-1,-1,-1,-1,-1}},R);
    --conicCoeffs = transpose matrix{{a,b,c,d,e}};
    xsq := ptwsProd(x,x);
    xy := ptwsProd(x,y);
    ysq := ptwsProd(y,y);
    sourceMat := transpose matrix{entries xsq, entries xy, entries ysq, entries x, entries y};
    --sol = inverse(sourceMat)*ones
    --polys = entries (sourceMat*conicCoeffs-ones);
    sol := solve(sourceMat,ones)
)
