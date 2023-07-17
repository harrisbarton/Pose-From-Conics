restart
load "reconstructConic.m2";
load "smallFuncs.m2";
R = frac(QQ[a,b,c,d,e,f,t1, t2, t3]);
-- arbitrary parameters
(a,b,c,t1,t2,t3) = (1,2,3,4,5,6);
(n1, n2, n3) = (7,5,-3);
-- creating P_4x4 
A = matrix{{0, a, b}, {-a, 0, c}, {-b, -c, 0}};
I = id_(R^3);
rotationCayley =(I + A) * inverse(I - A);
translation = transpose(matrix{{t1,t2,t3}});
zeros = matrix{{0,0,0,1}};
P4by4 = (rotationCayley | translation) || zeros;
-- defining P3 points
ones = matrix{{1,1,1,1,1}};
randMatrix = sub(random(ZZ^2, ZZ^5),R);
randMatrixT = transpose randMatrix;
thirdRow = transpose  (matrix{n1*randMatrixT_0 + n2*randMatrixT_1}  - transpose ones ) * -(1/n3);
pointMatrix =  randMatrix || thirdRow  || ones;
pointMatrixT = transpose pointMatrix;
--reconstruct world conic
Cw = reconstructConic(pointMatrixT_0, pointMatrixT_1);
-- project world points
    -- P3by3 = submatrix'(P4by4, {3}, {2,3}) | (matrix{r3}+translation);
projectedPoints = P4by4*pointMatrix;
projectedPointsT = transpose projectedPoints;
firstCol = ptwsDiv(projectedPointsT_0,projectedPointsT_2);
secondCol = ptwsDiv(projectedPointsT_1, projectedPointsT_2);
projectedPointsLegit = matrix{toList firstCol, toList secondCol} || ones;
projectedPointsLegitT = transpose projectedPointsLegit;
-- reconstruct conic in image plane
Cim = reconstructConic(projectedPointsLegitT_0, projectedPointsLegitT_1);


