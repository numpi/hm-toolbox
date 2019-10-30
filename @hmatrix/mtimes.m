function H = mtimes(H1, H2)
if isfloat(H2) % H * v
    if isscalar(H2) || all(size(H2) == 1)
	H = H1;
        if is_leafnode(H1)
	    if H1.admissible
		H.U = H1.U * H2;
	    else	
                H.F = H1.F * H2;
	    end
        else
            H.A11 = H1.A11 * H2;
            H.A22 = H1.A22 * H2;
            H.A21 = H1.A21 * H2;
            H.A12 = H1.A12 * H2;
        end
    else
        H = hmatrix_mtimes_dense(H1, H2);
    end
    
    return;
elseif isfloat(H1) % v * H
	H = (H2' * H1)';
elseif isa(H1, 'hmatrix') && isa(H2, 'hmatrix') % Multiplication of two matrices
	H = hmatrix_mtimes(H1, H2);
end

end
