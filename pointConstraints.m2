restart
fivePointReconstruction = (Cim, m1, m2, m3, matrixOfPoints) -> (
    -- Create points in worlds
    matrixOfPointsT = transpose(matrixOfPoints);
    thirdRow = transpose  (matrix{m1*matrixOfPointsT_0 + m2*matrixOfPointsT_1} - transpose ones ) * -(1/m3);
    matrixOfPointsRowDrop = submatrix'(matrixOfPoints, {2}, );
    pointMatrix =  matrixOfPointsRowDrop || thirdRow  || ones;
    pointMatrixT = transpose pointMatrix;
    -- Projection
    projectedPoints = P4by4*pointMatrix;
    projectedPointsT = transpose projectedPoints;
    -- Create points in the image plane
    firstCol = ptwsDiv(projectedPointsT_0,projectedPointsT_2);
    secondCol = ptwsDiv(projectedPointsT_1, projectedPointsT_2);
    projectedPointsLegit = matrix{toList firstCol, toList secondCol} || ones;
    projectedPointsLegitT = transpose projectedPointsLegit;
    -- Make constrains
    polysol = entries(projectedPointsLegitT*Cim*projectedPointsLegit)
    )
end
