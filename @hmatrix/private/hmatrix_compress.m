function H = hmatrix_compress(H, nrm)
%COMPRESS_HODLRATRIX Recursive compression of an hmatrix object.

if ~exist('nrm', 'var')
    nrm = norm(H, hmatrixoption('norm'));
end

if ~is_leafnode(H)
    H.A11 = hmatrix_compress(H.A11, nrm);
    H.A22 = hmatrix_compress(H.A22, nrm);
    H.A12 = hmatrix_compress(H.A12, nrm);
    H.A21 = hmatrix_compress(H.A21, nrm);
else
    % No compression needed for non-admissible blocks
    if H.admissible
        [H.U, H.V] = compress_factors(H.U, H.V, nrm);
    end
end

end

