function H = hodlr_hadamard_mul(H1, H2, compress);
if ~exist('compress', 'var')
    compress = true;
end
H = hodlr_hadamard_mul_ric(H1, H2);
if compress
	    H = compress_hodlr(H);
end
end

function H = hodlr_hadamard_mul_ric(H1, H2)
H = H1;
if xor(~is_leafnode(H1), ~is_leafnode(H2))
    error('HODLR_HADAMARD_MUL:: the two hodlr matrices have not compatible partitioning')
end
if is_leafnode(H1)
    H.F = H1.F .* H2.F;
else
    if size(H1.U21, 1) ~= size(H2.U21, 1) || size(H1.U12, 1) ~= size(H2.U12, 1)|| ...
            size(H1.V12, 1) ~= size(H2.V12, 1) || size(H2.V21, 1) ~= size(H2.V21, 1)
        error('HODLR_HADAMARD_MUL:: the two hodlr matrices have not compatible partitioning')
    end

    H.U12 = zeros(size(H1.U12, 1), size(H1.U12, 2) * size(H2.U12, 2));
    H.V12 = zeros(size(H1.V12, 1), size(H1.V12, 2) * size(H2.V12, 2));
    k = size(H2.U12, 2);
    for j = 1 : size(H1.U12,2)
        H.U12(:,(j-1)*k+1:j*k) = H1.U12(:,j) .* H2.U12;
        H.V12(:,(j-1)*k+1:j*k) = H1.V12(:,j) .* H2.V12;
    end

    H.U21 = zeros(size(H1.U21, 1), size(H1.U21, 2) * size(H2.U21, 2));
    H.V21 = zeros(size(H1.V21, 1), size(H1.V21, 2) * size(H2.V21, 2));
    k = size(H2.U21, 2);
    for j = 1 : size(H1.U21,2)
        H.U21(:,(j-1)*k+1:j*k) = H1.U21(:,j) .* H2.U21;
        H.V21(:,(j-1)*k+1:j*k) = H1.V21(:,j) .* H2.V21;
    end

    H.A11 = hodlr_hadamard_mul_ric(H1.A11, H2.A11);
    H.A22 = hodlr_hadamard_mul_ric(H1.A22, H2.A22);
end

end
