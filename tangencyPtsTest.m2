restart
load "construct_matrix.m2"
load "findTangencyPts.m2"
t = transpose matrix{{t1/t3,t2/t3}};
Tw = findTangencyPts Cw;
Tim = findTangencyPts Cim;
--Tw = Ptsw_0 | Ptsw_1;
TwT = transpose Tw;
--Tim = Ptsim_0 | Ptsw_1;
TimT = transpose Tim;
ones = transpose matrix{{1,1}};
P = sub(P4by4,CC)
-- Crating Points in PP^3
thirdRow = transpose  (matrix{n1*TwT_0 + n2*TwT_1}  - ones ) * -(1/n3);
TwRowDrop = submatrix'(Tw,{2},);
realTw = TwRowDrop || thirdRow || transpose(ones)
-- Constraint
v3 = (transpose(P*realTw))_2
LHS = submatrix'(P*realTw,{2,3},)
--zeros = transpose matrix{{0,0}};
tmat = t|t
RH = submatrix'(Tim,{2},)+tmat
firstR = matrix ptwsProd((transpose(RH))_0,v3)
secondR = matrix ptwsProd((transpose(RH))_1,v3)
RHS = transpose(firstR|secondR)
LHS
RHS-LHS
