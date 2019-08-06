function H = compress_hodlr(H, nrm)
%COMPRESS_HODLRATRIX Recursive compression of an hodlr object.
%

if ~exist('nrm', 'var')
    nrm = norm(H, 2);
end

if ~is_leafnode(H)
    H.A11 = compress_hodlr(H.A11, nrm);
    H.A22 = compress_hodlr(H.A22, nrm);
    
    [H.U21, H.V21] = compress_factors(H.U21, H.V21, nrm);
    [H.U12, H.V12] = compress_factors(H.U12, H.V12, nrm);
end

end

