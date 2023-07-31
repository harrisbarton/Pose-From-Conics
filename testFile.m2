restart 
load "Data/construct_matrix.m2"
load "Methods/mainFuncs/createPointConstraints.m2"
load "Methods/mainFuncs/createConicConstraints.m2"
load "Methods/mainFuncs/createTangencyConstraints.m2"
load "Methods/mainFuncs/createPolarityConstraints.m2"
load "Methods/basicFuncs/generateRandPointsFromConic.m2"

R = ZZ/911[r1,r2,r3,l1,l2,l3];
-- This parametrization has a common denominator which is (r1^2+r2^2+r3^2+1)
-- Here we are just ignoring it since it will be cancelled
rotationCayley = matrix{{-r1^2-r2^2+r3^2+1, -2*r2*r3-2*r1, 2*r1*r3-2*r2}, 
                        {-2*r2*r3 + 2*r1, -r1^2+r2^2-r3^2+1, -2*r1*r2-2*r3},
     			{2*r1*r3 + 2*r2, -2*r1*r2 + 2*r3,r1^2-r2^2-r3^2+1}}
translation = transpose(matrix{{l1,l2,l3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros
--P3by3 = submatrix'(P4by4*zoomOut,{3},)

(pw1,pim1,Cw1,Cim1,wPts1,imPts1,zoomOut1) = generateConicData()
(pw2,pim2,Cw2,Cim2,wPts2,imPts2,zoomOut2) = generateConicData()
P3by3First = submatrix'(P4by4*sub(zoomOut1,R),{3},)
P3by3Second = submatrix'(P4by4*sub(zoomOut2,R),{3},)
--------------------------------Inputs------------------------------------
-- Two conics
Aw1 = sub(Cw1,R)
Aw2 = sub(Cw2,R)
Aim1 = sub(Cim1, R)
Aim2 = sub(Cim2, R)
--(3,3) case
randwPts1P3P = submatrix'(sub(wPts1,R),,{3,4});
randwPts2P3P = submatrix'(sub(wPts2,R),,{3,4});
-- (4,2) case
randwPts1P42P = submatrix'(sub(wPts1,R),,{3,4});
randwPts2P42P = submatrix'(sub(wPts2,R),,{3,4});
-- (5,1) case
randwPts1P51P = submatrix'(sub(wPts1,R),,{3,4});
randwPts2P51P = submatrix'(sub(wPts2,R),,{3,4});


-- Point correspondences
pW1 = sub(pw1,R)
pW2 = sub(pw2,R)
pIm1 = sub(pim1,R)
pIm1 = sub(pim2,R)


--------------------------Point constraints in (3,3) case-------------------------
ptsConstraints1P3P = createPointConstraints(Aim1,P3by3First,randwPts1P3P);
ptsConstraints2P3P = createPointConstraints(Aim2,P3by3Second,randwPts2P3P);
JP3P = ideal flatten(ptsConstraints1P3P,ptsConstraints2P3P);
numgens JP3P
time GP3P = groebnerBasis (JP3P,Strategy=>"F4");
dim ideal leadTerm GP3P
degree ideal leadTerm GP3P

dim JP3P
gb JP3P

--------------------------Point constraints in (4,2) case-------------------------
ptsConstraints1P42P = createPointConstraints(Aim1,P3by3First,randwPts1P42P);
ptsConstraints2P42P = createPointConstraints(Aim2,P3by3Second,randwPts2P42P);
JP42P = ideal flatten(ptsConstraints1P42P,ptsConstraints2P42P);
numgens JP42P
time GP42P = groebnerBasis (JP42P,Strategy=>"F4");
dim ideal leadTerm GP42P
degree ideal leadTerm GP42P

dim JP42P
gb JP42P



--------------------------Point constraints in (5,1) case-------------------------
ptsConstraints1P51P = createPointConstraints(Aim1,P3by3First,randwPts1P51P);
ptsConstraints2P51P = createPointConstraints(Aim2,P3by3Second,randwPts2P51P);
JP51P = ideal flatten(ptsConstraints1P51P,ptsConstraints2P51P);
numgens JP51P
time GP51P = groebnerBasis (JP51P,Strategy=>"F4");
dim ideal leadTerm GP51P
degree ideal leadTerm GP51P

dim JP51P
gb JP51P



--------------------------Constraints from two conics----------------------
conicConstraints1 = createConicConstraints(Aim1,Aw1,P3by3First,k1);
conicConstraints2 = createConicConstraints(Aim2,Aw2,P3by3Second,k1);
J2 = ideal flatten(conicConstraints1,conicConstraints2);
numgens J2
time G2 = groebnerBasis (J2,Strategy=>"F4");
dim ideal leadTerm G2
dim J2
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
