function X = small_care_solve(A, B, C, tol, maxit)
	% Solve A' X + X A' - X B B' X + C = 0 wiht a dense method
	X = newton_care(A, B * B', C, zeros(size(A)), tol, maxit); 
end 
