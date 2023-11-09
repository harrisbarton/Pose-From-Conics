restart 
load "Methods/degPtsConstraints.m2"

R = ZZ/911[a,b,c,t1,t2,t3];

-- Experiments with different point constraints
degPtsConstraints(R, {4,2})

