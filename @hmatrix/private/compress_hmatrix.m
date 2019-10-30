function H = compress_hmatrix(H, nrm)
%COMPRESS_HODLRATRIX Recursive compression of an hmatrix object.

if ~exist('nrm', 'var')
    nrm = norm(H, 2);
end

if ~is_leafnode(H)
    H.A11 = compress_hmatrix(H.A11, nrm);
    H.A22 = compress_hmatrix(H.A22, nrm);    
    H.A12 = compress_hmatrix(H.A12, nrm);
    H.A21 = compress_hmatrix(H.A21, nrm);    
else
    % No compression needed for non-admissible blocks
    if H.admissible    
        [H.U, H.V] = compress_factors(H.U, H.V, nrm);
    end    
end

end

