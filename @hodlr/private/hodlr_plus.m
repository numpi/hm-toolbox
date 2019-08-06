function H = hodlr_plus(H1, H2, compress)
%HODLRATRIX_PLUS Sum of two HODLR matrices
if ~exist('compress', 'var')
    compress = true;
end
H =  hodlr_plus_ric(H1, H2);

if compress
    H = compress_hodlr(H);
end
end

function H = hodlr_plus_ric(H1, H2)
H = H1;

if is_leafnode(H1)
    H.F = H1.F + H2.F;
else
    H.A11 = hodlr_plus_ric(H1.A11, H2.A11);
    H.A22 = hodlr_plus_ric(H1.A22, H2.A22);
    
    H.U12 = [ H1.U12, H2.U12 ];
    H.V12 = [ H1.V12, H2.V12 ];
    
    H.U21 = [ H1.U21, H2.U21 ];
    H.V21 = [ H1.V21, H2.V21 ];
end
end

