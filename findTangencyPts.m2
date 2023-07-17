restart
needsPackage "NumericalAlgebraicGeometry";

findTangencyPts = (Cw) -> (
    R := QQ[x,y];
    Pr4by4 := sub(P4by4,R);
    Aw := sub(Cw,R);
    u := transpose matrix{{x,y,1}};
    uT := transpose u;
    u0 := transpose matrix{{x,y}};
    u0T := transpose u0;
    -- Solve for the tangency points
    qw := first flatten entries(uT*Aw*u);
    Jacw := matrix(diff(u0,qw));
    JacEq := first flatten entries(u0T*Jacw);
    solw1 := flatten entries matrix first solveSystem({qw,JacEq});
    tangencyPt1 := transpose matrix{{solw1_0,solw1_1,1}};
    solw2 := flatten entries matrix ((solveSystem({qw,JacEq}))_1);
    tangencyPt2 := transpose matrix{{solw2_0,solw2_1,1}};
    return tangencyPt1 | tangencyPt2;
    )
end
