restart 
load "construct_matrix.m2"
load "createPointConstraints.m2"
load "generateRandPointsFromConic.m2"
load "findTangencyPts.m2"

R = QQ[r1,r2,r3,l1,l2,l3];
-- This parametrization has a common denominator which is (r1^2+r2^2+r3^2+1)
rotationCayley = matrix{{-r1^2-r2^2+r3^2+1, -2*r2*r3-2*r1, 2*r1*r3-2*r2}, 
    {-2*r2*r3 + 2*r1, -r1^2+r2^2-r3^2+1, -2*r1*r2-2*r3},
     { 2*r1*r3 + 2*r2, -2*r1*r2 + 2*r3,r1^2-r2^2-r3^2+1}}
translation = transpose(matrix{{l1,l2,l3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros

-------------------Constraints from points---------------------------------
Aim = sub(Cim, R);
m1 = sub(n1, R);
m2 = sub(n2, R);
m3 = sub(n3, R);
-- Five pts
matrixOfPoints = sub(randMatrix, R);
-- Six pts
-*
sixthPt = sub(generateRandRealPointsFromConic(Cw),R);
Minterm = sub(randMatrix, R);
matrixOfPoints = Minterm | submatrix'(sixthPt,{2},);
*-
-- One Pt
--matrixOfPoints = sub(submatrix'(randMatrix, ,{1,2,3}),R);
output = flatten createPointConstraints(Aim, m1, m2, m3, matrixOfPoints);
J = ideal output;
numgens J

-----------------Constraints from tangency points--------------------------
findTangencyPts Cw
