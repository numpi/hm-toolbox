function x = hmatrix_mtimes_dense(H, v)
	if is_leafnode(H)
		if H.admissible
			x = H.U * (H.V' * v);
		else
			x = H.F * v;
		end
	else
		n1 = size(H.A11, 2);
		n2 = size(H.A22, 2);
		x = [H.A11 * v(1:n1, :) + H.A12 * v(n1+1:n1+n2, :);...
		H.A21 * v(1:n1, :) + H.A22 * v(n1+1:n1+n2, :)]; 
	end
end
