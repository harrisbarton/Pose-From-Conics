restart
load "Methods/helpers.m2"
createTangencyConstraints = (Cw,P,t,l3,option) -> (
    Tw = findTangencyPts Cw;
    Tim = findTangencyPts Cim;
    TwT = transpose Tw;
    Tw = findTangencyPts Cw;
    ones = transpose matrix{{1,1}};    
    -- Crating Points in PP^3
    thirdRow = transpose  (matrix{n1*TwT_0 + n2*TwT_1}  - ones ) * -(1/n3);
    TwRowDrop = submatrix'(Tw,{2},);
    realTw = TwRowDrop || thirdRow || transpose(ones);
    -- Left Hand Side
    v3 = (transpose(P*realTw))_2;
    LHS = cd*l3*submatrix'(P*realTw,{2,3},);
    -- Right Hand Side
    tmat = t|t;
    RH = l3*submatrix'(Tim,{2},)+tmat;
    firstR = matrix ptwsProd((transpose(RH))_0,v3);
    secondR = matrix ptwsProd((transpose(RH))_1,v3);
    if (option == 0) then(
       RHS = cd*transpose(firstR|secondR);
       )else(
       RHSold = cd*transpose(firstR|secondR);
       RHS = matrix(RHSold_1) | matrix(RHSold_0);
       );
    eqns1 = flatten entries(RHS-LHS);
    -- 2 more constraints from maping the origins
    wOrgn = transpose matrix{{0,0,1/n3,1}};
    u3 = first (entries(P*wOrgn))_2;
    LHS2 = l3*submatrix'((P*wOrgn),{2,3},);
    RHS2 = cd*(u3*t);
    eqns2 = flatten entries(RHS2-LHS2);
    return flatten {eqns1,eqns2}
    )
