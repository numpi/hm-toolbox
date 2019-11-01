function H = hodlr2hmatrix(H1)
%HODLR2HMATRIX Convert a HODLR matrix into HMATRIX format. 

H = hmatrix;
H.sz = H1.sz;

if is_leafnode(H1)
    H.F = H1.F;
    H.admissible = false;
else
    H.A11 = hodlr2hmatrix(H1.A11);
    H.A22 = hodlr2hmatrix(H1.A22);
    
    [m1, n1] = size(H.A11);
    [m2, n2] = size(H.A22);
    
    H.A12 = hmatrix;
    H.A12.admissible = true;
    H.A12.sz = [m1 n2];
    H.A12.U = H1.U12;
    H.A12.V = H1.V12;
    
    H.A21 = hmatrix;
    H.A21.admissible = true;
    H.A21.sz = [m2 n1];
    H.A21.U = H1.U21;
    H.A21.V = H1.V21;
end

