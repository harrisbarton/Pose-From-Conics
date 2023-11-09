restart
load "Methods/helpers.m2";

--S = frac(FF[a,b,c,t1, t2, t3]);

generateConicData = () -> (
FF := ZZ/911;
-- arbitrary parameters
--(a,b,c,t1,t2,t3) = (1,2,3,-1/2,-1/2,-1/2)
(a,b,c,t1,t2,t3) := toSequence apply(6,i->random FF);
--(n1, n2, n3) = (7,5,-3)
(n1,n2,n3) := toSequence apply(3,i->random FF);
-- creating P_4x4 
S := matrix{{0, a, b}, {-a, 0, c}, {-b, -c, 0}};
I := id_(FF^3);
RCayley := (I + S) * inverse(I - S);
t := transpose(matrix{{t1,t2,t3}});
zeroOne := matrix{{0,0,0,1}};
P4by4 := (RCayley | t) || zeroOne;
-- creating P_3x3
T := matrix{{1,0,0},{0,1,0},{-n1/n3,-n2/n3,1/n3},{0,0,1}};
P3by3 := submatrix'(P4by4*T,{3},);
-- defining P3 points
ones := matrix{{1,1,1,1,1}};
pw := sub(random(ZZ^2, ZZ^5),FF) || ones;
pwT := transpose pw;
--reconstruct world conic
Cw := reconstructConic(pwT_0, pwT_1);
-- project world points
projectedPoints := P3by3*pw;
pim := zoomIn(projectedPoints);
pimT := transpose pim;
-- reconstruct conic in image plane
Cim := reconstructConic(pimT_0, pimT_1);
-- Take a random point on the plane n1xw+n2yw+n3zw=1
wPts := random(FF^2,FF^1) || matrix{{1}};
imPts := zoomIn(P3by3*pw);
return (wPts,imPts,Cw,Cim,pw,pim,T);
)

--(wPts,imPts,Cw,Cim,pw,pim,T) = generateConicData()
end
