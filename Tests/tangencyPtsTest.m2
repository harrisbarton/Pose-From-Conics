restart
load "Data/construct_matrix.m2"
load "Methods/basicFuncs/findTangencyPts.m2"
load "Methods/basicFuncs/zoomIn.m2"

P3by3T = transpose P3by3
P3by3InvT = transpose(inverse(P3by3))
k = (Cw*pw)_(0,0)/(P3by3T*Cim*pim)_(0,0)
Cw*pw-k*P3by3T*Cim*pim

-*
t = transpose matrix{{t1/t3,t2/t3}};
Tw = findTangencyPts Cw;
Tim = findTangencyPts Cim;
--Tw = Ptsw_0 | Ptsw_1;
TwT = transpose Tw;
--Tim = Ptsim_0 | Ptsw_1;
TimT = transpose Tim;
ones = transpose matrix{{1,1}};
P = sub(P3by3,CC)
-- Constraint
v3 = (transpose(P*Tw))_2
LHS = submatrix'(P*Tw,{2,3},)
--zeros = transpose matrix{{0,0}};
tmat = t|t
RH = submatrix'(Tim,{2},)+tmat
firstR = matrix ptwsProd((transpose(RH))_0,v3)
secondR = matrix ptwsProd((transpose(RH))_1,v3)
RHS = transpose(firstR|secondR)
LHS
RHS-LHS
*-
