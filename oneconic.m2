restart
load "reconstructConic.m2"
load "smallFuncs.m2"

R = frac(QQ[a,b,c,d,e,f,t1, t2, t3])
-- arbitrary parameters
(a,b,c,t1,t2,t3) = (1,2,3,1/2,1/2,1/2)
(n1, n2, n3) = (7,5,-3)
-- creating P_4x4 
A = matrix{{0, a, b}, {-a, 0, c}, {-b, -c, 0}}
I = id_(R^3);
rotationCayley =(I + A) * inverse(I - A)
translation = transpose(matrix{{t1,t2,t3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros
-- defining P3 points
ones = matrix{{1,1,1,1,1}}
randMatrix = sub(random(ZZ^2, ZZ^5),R)

randMatrix
randMatrixT = transpose randMatrix
thirdRow = transpose  (matrix{n1*randMatrixT_0 + n2*randMatrixT_1}  - transpose ones ) * -(1/n3)
pointMatrix =  randMatrix || thirdRow  || ones
pointMatrixT = transpose pointMatrix
--reconstruct world conic
Cw = reconstructConic(pointMatrixT_0, pointMatrixT_1)
-- project world points
    -- P3by3 = submatrix'(P4by4, {3}, {2,3}) | (matrix{r3}+translation);
projectedPoints = P4by4*pointMatrix
projectedPointsT = transpose projectedPoints
firstCol = ptwsDiv(projectedPointsT_0,projectedPointsT_2)
secondCol = ptwsDiv(projectedPointsT_1, projectedPointsT_2)
projectedPointsLegit = matrix{toList firstCol, toList secondCol} || ones
projectedPointsLegitT = transpose projectedPointsLegit
-- reconstruct conic in image plane
Cim = reconstructConic(projectedPointsLegitT_0, projectedPointsLegitT_1)

v=projectedPointsLegit_0

projectedPointsLegit

transpose(matrix(v))*Cim*v



T = QQ[x,y]
Pt  = sub(P4by4,T)
Cim2 = sub(Cim,T) 
Cw2 = sub(Cw, T) 



u  = matrix({{x},{y},{1}})
u0 = matrix({{x},{y}})


r  = matrix({{x},{y},{1}})
r0 = matrix({{x},{y}})


qw =flatten(entries(transpose(u)*Cw2*u))
qw=qw_0
qim=flatten(entries(transpose(u)*Cim2*u))
qim=qim_0

diff( transpose vars T, qw) 

i0 = u0-matrix({{t1/t3},{t2/t3}})

jacim  = flatten(entries(transpose(i0)*matrix(diff(transpose vars T ,qim))))
jacim = jacim_0
jacw  = flatten(entries(transpose(u0)*matrix(diff(transpose vars T ,qw))))
jacw = jacw_0

Jim  = {qim,jacim}

Jw = {qw,jacw}


solw = solveSystem({qw,jacw})

tt = flatten entries matrix solw_0
vv = matrix {{tt_0}, {tt_1},{1}}
transpose(vv)*Cw*vv

Pt = sub(Pt,CC)

solw

projectedsols={}

for  sol in solw do( 
    sol1 = flatten  entries matrix sol ;
    newt = Pt * matrix {{sol1_0},{sol1_1},{-(1/n3) * (n1*sol1_0 + n2*sol1_1 -1) },{1}};
    newt2 = {newt_(0,0)/newt_(2,0),newt_(1,0)/newt_(2,0)};
    projectedsols = append(projectedsols,newt2);
    )

imagesols = solveSystem({qim,jacim})


