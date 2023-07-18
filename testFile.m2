restart 
load "construct_matrix.m2"
load "createPointConstraints.m2"
load "createTangencyConstraints.m2"
load "generateRandPointsFromConic.m2"
load "findTangencyPts.m2"

-------------------------------Problem Setups-----------------------------
R = QQ[r1,r2,r3,l1,l2,l3];
-- This parametrization has a common denominator which is (r1^2+r2^2+r3^2+1)
-- Here we are just ignoring it since it will be cancelled
rotationCayley = matrix{{-r1^2-r2^2+r3^2+1, -2*r2*r3-2*r1, 2*r1*r3-2*r2}, 
                        {-2*r2*r3 + 2*r1, -r1^2+r2^2-r3^2+1, -2*r1*r2-2*r3},
     			{2*r1*r3 + 2*r2, -2*r1*r2 + 2*r3,r1^2-r2^2-r3^2+1}}
translation = transpose(matrix{{l1,l2,l3}})
zeros = matrix{{0,0,0,1}}
P4by4 = (rotationCayley | translation) || zeros

--------------------------Constraints from points-------------------------
-- Inputs
Aim = sub(Cim, R);
m1 = sub(n1, R);
m2 = sub(n2, R);
m3 = sub(n3, R);
-- Five pts
randMatrix = random(R^3,R^5);

output = flatten createPointConstraints(Aim, m1, m2, m3, randMatrix);
J = ideal output;
numgens J

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
