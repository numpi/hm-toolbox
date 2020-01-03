function H = hmatrix_hadamard_mul(H1, H2, compress)
if ~exist('compress', 'var')
    compress = true;
end
H = hmatrix_hadamard_mul_ric(H1, H2);
if compress
	    H = compress_hmatrix(H);
end
end

function H = hmatrix_hadamard_mul_ric(H1, H2)
H = H1;
if xor(~is_leafnode(H1), ~is_leafnode(H2))
    error('hmatrix_HADAMARD_MUL:: the two hmatrix matrices have not compatible partitioning')
end
if is_leafnode(H1)
	if any(H1.sz ~= H2.sz)
        	error('hmatrix_HADAMARD_MUL:: the two hmatrix matrices have not compatible partitioning')
    	end
	if H1.admissible && H1.admissible
		H.U = zeros(size(H1.U, 1), size(H1.U, 2) * size(H2.U, 2));
        H.V = zeros(size(H1.V, 1), size(H1.V, 2) * size(H2.V, 2));
        k = size(H2.U, 2);
        for j = 1 : size(H1.U,2)
            H.U(:,(j-1)*k+1:j*k) = H1.U(:,j) .* H2.U;
            H.V(:,(j-1)*k+1:j*k) = H1.V(:,j) .* H2.V;
        end
	else
    		H.F = full(H1) .* full(H2);
	end
else
    H.A11 = hmatrix_hadamard_mul_ric(H1.A11, H2.A11);
    H.A22 = hmatrix_hadamard_mul_ric(H1.A22, H2.A22);
    H.A12 = hmatrix_hadamard_mul_ric(H1.A12, H2.A12);
    H.A21 = hmatrix_hadamard_mul_ric(H1.A21, H2.A21);
end

end
