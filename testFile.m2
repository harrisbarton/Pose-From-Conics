restart 
load "Data/construct_matrix.m2"
load "Methods/mainFuncs/createPointConstraints.m2"
load "Methods/mainFuncs/createConicConstraints.m2"
load "Methods/mainFuncs/createTangencyConstraints.m2"
load "Methods/mainFuncs/createPolarityConstraints.m2"
load "Methods/basicFuncs/generateRandPointsFromConic.m2"

R = ZZ/911[r1,r2,r3,l1,l2,l3];

degPtsConstraints=(R, L) -> (
   -- This parametrization has a common denominator which is (r1^2+r2^2+r3^2+1)
   -- Here we are just ignoring it since it will be cancelled    
   rotationCayley := matrix{{-r1^2-r2^2+r3^2+1, -2*r2*r3-2*r1, 2*r1*r3-2*r2}, 
                        {-2*r2*r3 + 2*r1, -r1^2+r2^2-r3^2+1, -2*r1*r2-2*r3},
     			{2*r1*r3 + 2*r2, -2*r1*r2 + 2*r3,r1^2-r2^2-r3^2+1}}; 	
   translation := transpose(matrix{{l1,l2,l3}});
   zeros := matrix{{0,0,0,1}};
   P4by4 := (rotationCayley | translation) || zeros;
   conicData := apply(#L, i -> generateConicData());
   -- (pw,pim,Cw,Cim,wPts,imPts,zoomOut) = generateConicData();
   -- worldConics = apply(#conicData, i -> sub(conicData#i#2,R));
   imageConics := apply(#conicData, i -> sub(conicData#i#3,R));
   P3by3Matrices := apply(#conicData, i -> submatrix'(P4by4*sub(conicData#i#6,R),{3},));
   randwPts := apply(#conicData, i -> submatrix(sub(conicData#i#4,R),,{0..(L#i-1)}));
   ptsConstraints := apply(#imageConics, i -> createPointConstraints(imageConics#i,P3by3Matrices#i,randwPts#i));
   J := ideal flatten ptsConstraints;
   G := ideal leadTerm groebnerBasis(J, Strategy=>"F4");
   (dim G, degree G)
)
degPtsConstraints(R, {4,2})

