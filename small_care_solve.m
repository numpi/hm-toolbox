function X = small_care_solve(A, B, C, tol, maxit)
	X = small_solve_ME(A', B, C, [], [], tol, maxit); 
end 
