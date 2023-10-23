restart
load "Methods/basicFuncs/smallFuncs.m2"

zoomIn = (pts) -> (
    n := rank source pts;
    rowOf1s := mutableMatrix(ZZ,1,n);
    for i from 0 to n-1 do (rowOf1s_(0,i) = 1);
    rowOfOnes := matrix rowOf1s;
    ptsT := transpose pts;
    c := ptsT_2;
    firstRow := ptwsDiv(ptsT_0,c);
    secondRow := ptwsDiv(ptsT_1,c);
    ptsZoomIn := matrix{toList firstRow, toList secondRow} || rowOfOnes
    )
