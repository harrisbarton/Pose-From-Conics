restart
load "smallFuncs.m2";
createPointConstraints = (Cim, m1, m2, m3, matrixOfPoints) -> (
    n = rank source matrixOfPoints;
    rowOfOnes = mutableMatrix(ZZ,1,n);
    for i from 0 to n-1 do (rowOfOnes_(0,i) = 1);
    rowOfOnes = matrix rowOfOnes;
    -- Create points in worlds
    matrixOfPointsT = transpose(matrixOfPoints);
    thirdRow = transpose  (matrix{m1*matrixOfPointsT_0 + m2*matrixOfPointsT_1} - transpose rowOfOnes) *(-1/m3);
    matrixOfPointsRowDrop = submatrix'(matrixOfPoints, {2}, );
    pointMatrix =  matrixOfPointsRowDrop || thirdRow  || rowOfOnes;
    pointMatrixT = transpose pointMatrix;
    -- Projection
    projectedPoints = P4by4*pointMatrix;
    projectedPointsT = transpose projectedPoints;
    -- Create points in the image plane
    firstCol = ptwsDiv(projectedPointsT_0,projectedPointsT_2);
    secondCol = ptwsDiv(projectedPointsT_1, projectedPointsT_2);
    projectedPointsLegit = matrix{toList firstCol, toList secondCol} || rowOfOnes;
    projectedPointsLegitT = transpose projectedPointsLegit;
    -- Make constrains
    polysol = diagEntries(projectedPointsLegitT*Cim*projectedPointsLegit)
    )
end
