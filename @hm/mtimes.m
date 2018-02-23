function H = mtimes(H1, H2)
%MTIMES Matrix multiplication

% Multiplication H * v
if isfloat(H2)
	if isscalar(H2)
		if H2 == 0
			H = hm('diagonal', zeros(size(H1,1), 1));
			return;
		end
		
		if ~isempty(H1.F)
			H = hm();
			H.F = H1.F * H2;
			H.sz = H1.sz;
		else
			H = H1;
			H.A11 = H1.A11 * H2;
			H.A22 = H1.A22 * H2;
			H.U21 = H1.U21 * H2;
			H.U12 = H1.U12 * H2;
		end
	else
		H = hmatrix_mtimes_dense(H1, H2);
	end
	
	return;
end

% Multiplication w' * H
if isfloat(H1)
	if isscalar(H1)
		H = H2 * H1;
    else
		H = dense_mtimes_hmatrix(H1, H2);
	end
	
	return;
end

% Multiplication of two H-matrices
H = hmatrix_mtimes(H1, H2);
H = compress_hmatrix(H);

end

