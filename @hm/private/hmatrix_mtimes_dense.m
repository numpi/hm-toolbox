function v = hmatrix_mtimes_dense(H1, v)
%HMATRIX_MTIMES_DENSE An H matrix times a dense matrix. 

if ~isempty(H1.F)
	v = H1.F * v;
else
	mp = H1.A11.sz(2);
	
	v = [ hmatrix_mtimes_dense(H1.A11, v(1:mp,:)) + ...
			H1.U12 * (H1.V12' * v(mp+1:end,:)) ; ...
		H1.U21 * (H1.V21' * v(1:mp,:)) + ...
			hmatrix_mtimes_dense(H1.A22, v(mp+1:end,:)) ];
end

end

