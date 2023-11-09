load "Data/construct_matrix.m2"
load "Methods/createPointConstraints.m2"
load "Methods/helpers.m2"

-- See the Experiment section in README for more details about the input
degPtsConstraints=(R, L) -> (
   -- This parametrization has a common denominator which is (a^2+b^2+c^2+1)
   -- Here we are just ignoring it since it will be cancelled    
   RCayley := matrix{{-a^2-b^2+c^2+1, -2*b*c-2*a, 2*a*c-2*b}, 
                           {-2*b*c + 2*a, -a^2+b^2-c^2+1, -2*a*b-2*c},
     			               {2*a*c + 2*b, -2*a*b + 2*c,a^2-b^2-c^2+1}}; 	
   t := transpose(matrix{{t1,t2,t3}});
   zeroOne := matrix{{0,0,0,1}};
   P4by4 := (RCayley | t) || zeroOne;
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