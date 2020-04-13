function X = hodlr_dac_uqme(A, B, C, tol, debug, nrmtype)
% solve the UQME A X^2 + B X  + C = 0 with divide & conquer
% sparseA, sparseB additional sparse version of A and B ---> speed up in building krylov subspaces
if ~exist('tol', 'var')
	tol = hodlroption('threshold') * 1e2;
end
if ~exist('debug', 'var')
	debug = false;
end
if ~exist('nrmtype', 'var')
	nrmtype = 2;
end
[n, m] = size(A);
X = hodlr();
X.sz = [n m];
if ~isempty(A.F)	
		X.F = cr_uqme(A.F, B.F, C.F); 
		if debug
			fprintf('Base case: Dimension: %d Residual = %e\n', n, norm(A.F * X.F^2 + B.F * X.F + C.F));
		end
else
		X0 = blkdiag(hodlr_dac_uqme(A.A11, B.A11, C.A11, tol, debug, nrmtype),  hodlr_dac_uqme(A.A22, B.A22, C.A22, tol, debug, nrmtype));
		AU = [zeros(size(A.U12,1), size(A.U21, 2)), A.U12; A.U21, zeros(size(A.U21,1), size(A.U12, 2))];
		AV = blkdiag(A.V21, A.V12);
		AV = X0' * (X0' * AV); 
		BU = [zeros(size(B.U12,1), size(B.U21, 2)), B.U12; B.U21, zeros(size(B.U21,1), size(B.U12, 2))]; 
		BV = blkdiag(B.V21, B.V12);
		BV = X0' * BV;
		CU = [zeros(size(C.U12,1), size(C.U21, 2)), C.U12; C.U21, zeros(size(C.U21,1), size(C.U12, 2))]; 
		CV = blkdiag(C.V21, C.V12);
		u = [AU, BU, CU]; v = [AV, BV, CV];
		[u, v] = compress_factors(u, v);
		[dXu, dXv] = ek_uqme(A, B, X0, u, v, min(1000, n), tol, debug, nrmtype);
		if debug
			dX = hodlr('low-rank', dXu, dXv); 
			fprintf('Correction equation: Dimension %d, Rank of the sol. = %d, Residual = %e\n', n, size(dXu, 2), norm((A * dX  + A * X0 + B) * dX + A * dX * X0 + hodlr('low-rank', u, v)))
        	end
		X = X0 + hodlr('low-rank', dXu, dXv);
		if debug 
			fprintf('Complete equation: Dimension %d, HODLR rank of the sol.  = %d,  Res = %e\n', n, hodlrrank(X), norm((A * X  + B) * X + C))
		end
		
end
