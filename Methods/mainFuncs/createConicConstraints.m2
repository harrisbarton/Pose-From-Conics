restart
createConicConstraints = (Cim,Cw,P3by3) -> (
    P3by3T := transpose P3by3;
    LHS := P3by3T*Cim*P3by3;
    conicConstraints := flatten entries (LHS-k*Cw)
    )
