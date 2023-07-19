restart
needsPackage "NumericalAlgebraicGeometry"

R = frac(QQ[a,b,k,c1,c2,c3]);
Q = matrix{{1/(a^2),0,0},{0,1/(b^2),0},{0,0,-1}};
C = matrix{{1,0,c1},{0,1,c2},{0,0,c3}};
tC = transpose C;
B = inverse(tC)*Q*inverse(C);
F = R[lambda];
charPoly = det(k*B-lambda*id_(R^3));
end

restart
(l1,l2,l3,a,b) = (1,2,3,4,5)
R = ZZ/7772777[k,c1,c2,c3,d];
load "oneConicMaSetup.m2"
dim I
degree I
d = 6
F = ZZ/7772777[k,c1,c2,c3];
load "oneConicMaSetup.m2"
dim I
degree I

needsPackage "NumericalAlgebraicGeometry"
F2 = CC[k,c1,c2,c3];
sols = solveSystem flatten(entries(gens(sub(I,F2))))
netList sols
end
