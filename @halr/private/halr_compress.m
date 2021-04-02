function H = halr_compress(H, nrm)
%COMPRESS_HODLRATRIX Recursive compression of an halr object.

if ~exist('nrm', 'var')
    nrm = norm(H, halroption('norm'));
end

if ~is_leafnode(H)
    H.A11 = halr_compress(H.A11, nrm);
    H.A22 = halr_compress(H.A22, nrm);
    H.A12 = halr_compress(H.A12, nrm);
    H.A21 = halr_compress(H.A21, nrm);
else
    % No compression needed for non-admissible blocks
    if H.admissible
        [H.U, H.V] = compress_factors(H.U, H.V, nrm);
    end
end

end

