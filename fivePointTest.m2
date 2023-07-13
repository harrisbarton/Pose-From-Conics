restart 
load "construct_matrix.m2"
R = QQ[r1,r2,r3,l1,l2,l3];
rotationCayley = matrix{{-r1^2-r2^2+r3^2+1, -2*r2*r3-2*r1, 2*r1*r3-2*r2}, 
    {-2*r2*r3 + 2*r1, -r1^2+r2^2-r3^2+1, -2*r1*r2-2*r3},
     { 2*r1*r3 + 2*r2, -2*r1*r2 + 2*r3,r1^2-r2^2-r3^2+1}}
translation = transpose(matrix{{l1,l2,l3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros

fivePointReconstruction = (Cim, m1, m2, m3, matrixOfPoints) -> (
    matrixOfPointsT = transpose(matrixOfPoints);
    thirdRow = transpose  (matrix{m1*matrixOfPointsT_0 + m2*matrixOfPointsT_1} - transpose ones ) * (1/m3);
    matrixOfPointsRowDrop = submatrix'(matrixOfPoints, {2}, );
    pointMatrix =  matrixOfPointsRowDrop || thirdRow  || ones;
    pointMatrixT = transpose pointMatrix;
    projectedPoints = P4by4*pointMatrix;
    projectedPointsT = transpose projectedPoints;
    firstCol = ptwsDiv(projectedPointsT_0,projectedPointsT_2);
    secondCol = ptwsDiv(projectedPointsT_1, projectedPointsT_2);
    projectedPointsLegit = matrix{toList firstCol, toList secondCol} || ones;
    projectedPointsLegitT = transpose projectedPointsLegit;
    polysol = entries(projectedPointsLegitT*Cim*projectedPointsLegit)
    )
Aim := sub(Cim, R);
m1 := sub(n1, R);
m2 := sub(n2, R);
m3 := sub(n3, R);
matrixOfPoints := sub(randMatrix, R);
output = flatten fivePointReconstruction(Aim, m1, m2, m3, matrixOfPoints)   
output_1 
J = ideal output
gens gb J
