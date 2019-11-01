function H = hmatrix_mtimes(H1, H2)
	if size(H1, 2) ~= size(H2, 1)
		error('HMATRIX_MTIMES:: incompatible dimensions')
	end
	H = hmatrix();
	H.sz = [H1.sz(1), H2.sz(2)];
	if is_leafnode(H1) || is_leafnode(H2)
		if H1.admissible
			H.admissible = true;
			if H2.admissible
				H.U = H1.U * (H1.V' * H2.U);
				H.V = H2.V;
			else 
				H.U = H1.U;
				H.V = H2' * H1.V;
			end
		elseif H2.admissible
			H.admissible = true;
			H.U = H1 * H2.U;
			H.V = H2.V;
		else
			H.F = full(H1) * full(H2);
		end
	else
		H.A11 = H1.A11 * H2.A11 + H1.A12 * H2.A21;
		H.A12 = H1.A11 * H2.A12 + H1.A12 * H2.A22;
		H.A21 = H1.A21 * H2.A11 + H1.A22 * H2.A21;
		H.A22 = H1.A21 * H2.A12 + H1.A22 * H2.A22;	
	end
end
