function x = chol_solve(F, b)
%ULV_SOLVE     solve the system A X = B where F contains the ULV factorization of A
%
%	       X = ULV_SOLVE(F, B) computes A\B with F = ULV(A);
	if (~isstruct(F))
		error('F is not of the correct type');
	end
	x = hss_chol_fact_solve(F, b);
end
