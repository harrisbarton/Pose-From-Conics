restart 
load "Data/construct_matrix.m2"
load "Methods/mainFuncs/createPointConstraints.m2"
load "Methods/mainFuncs/createConicConstraints.m2"
load "Methods/mainFuncs/createTangencyConstraints.m2"
load "Methods/mainFuncs/createPolarityConstraints.m2"
load "Methods/basicFuncs/generateRandPointsFromConic.m2"

-------------------------------Problem Setups-----------------------------
R = QQ[r1,r2,r3,l1,l2,l3,k1,k2];
-- This parametrization has a common denominator which is (r1^2+r2^2+r3^2+1)
-- Here we are just ignoring it since it will be cancelled
rotationCayley = matrix{{-r1^2-r2^2+r3^2+1, -2*r2*r3-2*r1, 2*r1*r3-2*r2}, 
                        {-2*r2*r3 + 2*r1, -r1^2+r2^2-r3^2+1, -2*r1*r2-2*r3},
     			{2*r1*r3 + 2*r2, -2*r1*r2 + 2*r3,r1^2-r2^2-r3^2+1}}
translation = transpose(matrix{{l1,l2,l3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros
P3by3 = submatrix'(P4by4*zoomOut,{3},)
--------------------------------Inputs------------------------------------
-- Two conics
Aw = sub(Cw,R)
Aim = sub(Cim, R)
-- Five pts
randwPts = sub(wPts,R);
-- Projection P_3x3
Proj3by3 = sub(P3by3, R)
pW = sub(pw,R)
pIm = sub(pim,R)
--------------------------Constraints from points-------------------------
ptsConstraints = createPointConstraints(Aim,Proj3by3,randwPts);
J = ideal ptsConstraints;
numgens J
gb J

--------------------------Constraints from one conic-----------------------
conicConstraints = createConicConstraints(Aim,Aw,Proj3by3,k1);
J2 = ideal conicConstraints;
numgens J2
--gens gb J2
--dim J2
-----------------Constraints from polarity---------------------------------
polarityConstraints = createPolarityConstraints(Aw,Aim,pW,pIm,Proj3by3,k2);
J3 = ideal polarityConstraints;
numgens J3

I = ideal flatten{conicConstraints,polarityConstraints};
numgens I
dim I
degree I


-*
-----------------Constraints from tangency points--------------------------
-- Inputs for the function
F = CC[r1,r2,r3,l1,l2,l3]
P = sub(P4by4,F)
t = transpose matrix{{l1,l2}};
cd = sub(r1^2+r2^2+r3^2+1,F);

-- Last argument means the option for mapping points
-- 0 means we map the tangency pts as the given order
-- 1 means we swap the tangency pts and then match them
constraints = createTangencyConstraints(Cw,P,t,l3,0);
constraintsSwap = createTangencyConstraints(Cw,P,t,l3,1);
J2 = ideal constraints;
--dim J2
--degree J2
numgens J2

-- Plugin data to check the system
L = CC;
F = ring constraints_0;
specification = {1,2,3,-1/2,-1/2,-1/2};
f = map(L,F,specification);
-- Apply the ring map
func = (poly) -> (f(poly));
lst = apply(constraints,func)
lst2 = apply(constraintsSwap,func)
*-
