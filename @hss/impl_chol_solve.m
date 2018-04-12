function x = impl_chol_solve(A, b)
%IMPL_CHOL_SOLVE     solve the system A X = B using an implicit generalized Cholesky of A
%
	x = hss_chol_solve(A, b);
end
