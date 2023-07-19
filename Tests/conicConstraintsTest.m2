restart
load "Data/construct_matrix.m2"
B = matrix{{1,0,0},{0,1,0},{-n1/n3,-n2/n3,1/n3},{0,0,1}}
--B = matrix{{1,0,0},{0,1,0},{0,0,0},{0,0,1}}
P4by4*B
P3by3 = submatrix'(P4by4*B,{3},)
P3by3T = transpose P3by3
LH = P3by3T*Cim*P3by3
k =LH_(0,0)/Cw_(0,0)
round(sub(LH_(0,0)-k*Cw_(0,0),QQ))
apply(flatten entries sub(LH-k*Cw,QQ),round)
end
