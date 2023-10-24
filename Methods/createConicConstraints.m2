restart
createConicConstraints = (Cim,Cw,P3by3,k) -> (
    P3by3T := transpose P3by3;
    LHS := P3by3T*Cim*P3by3;
    l := flatten entries (LHS-k*Cw);
    conicConstraints := l_{0,1,2,4,5}
    )
