function mH = uminus(H)
%UMINUS Change the sign of a HODLR Matrix
if is_leafnode(H)
    mH = H;
    mH.F = -H.F;
else
    mH = H;
    mH.A11 = -H.A11;
    mH.A22 = -H.A22;
    mH.U21 = -H.U21;
    mH.U12 = -H.U12;
end
