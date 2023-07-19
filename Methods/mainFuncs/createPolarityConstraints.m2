restart
createPolarityConstraints = (Cw,Cim,pw,pim,P3by3,k) -> (
    P3by3T := transpose P3by3;
    k = (P3by3T*Cw*pw)_(0,0)/(Cim*pim)_(0,0);
    polarityConstraints := flatten entries(Cw*pw-k*P3by3T*Cim*pim)
    )
end
