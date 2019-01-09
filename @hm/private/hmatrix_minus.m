function H = hmatrix_minus(H1, H2)
%HMATRIX_MINUS Difference of two HODLR matrices

H = H1;

if is_leafnode(H)
    H.F = H1.F - H2.F;
else
    H.A11 = hmatrix_minus(H1.A11, H2.A11);
    H.A22 = hmatrix_minus(H1.A22, H2.A22);
    
    H.U12 = [ H1.U12, -H2.U12 ];
    H.V12 = [ H1.V12, H2.V12 ];
    
    H.U21 = [ H1.U21, -H2.U21 ];
    H.V21 = [ H1.V21, H2.V21 ];
end


end

