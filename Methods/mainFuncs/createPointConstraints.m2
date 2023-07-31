restart
load "Methods/basicFuncs/smallFuncs.m2";
createPointConstraints = (Cim, P3by3, wPts) -> (
    n := rank source wPts;
    rowOf1s := mutableMatrix(FF,1,n);
    for i from 0 to n-1 do (rowOf1s_(0,i) = 1);
    rowOfOnes := matrix rowOf1s;
    -- Projection
    projectedPoints := P3by3*wPts;
    projectedPointsT := transpose projectedPoints;
    -- Create points in the image plane
    firstCol := ptwsDiv(projectedPointsT_0,projectedPointsT_2);
    secondCol := ptwsDiv(projectedPointsT_1, projectedPointsT_2);
    imPts := matrix{toList firstCol, toList secondCol} || rowOfOnes;
    imPtsT := transpose imPts;
    -- Make constrains
    apply(diagEntries(imPtsT*Cim*imPts),numerator)
    )
end
