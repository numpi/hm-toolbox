function X = mldivide(A,B)
	if isa(B,'hss')
		error('unsupported');
	else
		X = hss_solve(A, B);
	end
end
