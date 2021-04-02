function mH = uminus(H)
%UMINUS Change the sign of a HALR Matrix
mH = H;
if is_leafnode(H)
    if H.admissible
        mH.U = -H.U;
    else
        mH.F = -H.F;
    end
else
    mH.A11 = -H.A11;
    mH.A12 = -H.A12;
    mH.A21 = -H.A21;
    mH.A22 = -H.A22;
end
