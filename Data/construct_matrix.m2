restart
load "Methods/basicFuncs/reconstructConic.m2";
load "Methods/basicFuncs/smallFuncs.m2";
R = frac(QQ[a,b,c,t1, t2, t3]);
-- arbitrary parameters
(a,b,c,t1,t2,t3) = (1,2,3,-1/2,-1/2,-1/2)
(n1, n2, n3) = (7,5,-3)
-- creating P_4x4 
A = matrix{{0, a, b}, {-a, 0, c}, {-b, -c, 0}}
I = id_(R^3)
rotationCayley =(I + A) * inverse(I - A)
translation = transpose(matrix{{t1,t2,t3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros
-- creating P_3x3
zoomOut = matrix{{1,0,0},{0,1,0},{-n1/n3,-n2/n3,1/n3},{0,0,1}}
P3by3 = submatrix'(P4by4*zoomOut,{3},)
-- defining P3 points
ones = matrix{{1,1,1,1,1}}
wPts = sub(random(ZZ^2, ZZ^5),R) || ones
wPtsT = transpose wPts
--reconstruct world conic
Cw = reconstructConic(wPtsT_0, wPtsT_1)
-- project world points
projectedPoints = P3by3*wPts
projectedPointsT = transpose projectedPoints
firstCol = ptwsDiv(projectedPointsT_0,projectedPointsT_2)
secondCol = ptwsDiv(projectedPointsT_1, projectedPointsT_2)
imPts = matrix{toList firstCol, toList secondCol} || ones
imPtsT = transpose imPts
-- reconstruct conic in image plane
Cim = reconstructConic(imPtsT_0, imPtsT_1)
end
