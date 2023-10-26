restart 
load "Methods/degPtsConstraints.m2"

R = ZZ/911[r1,r2,r3,l1,l2,l3];
-- Experiments with different point constraints

degPtsConstraints(R, {4,2})

